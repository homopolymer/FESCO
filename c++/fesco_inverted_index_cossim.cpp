// routine to build inverted index for large-scale data
// blocked processing 

#include "Fesco.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include "LSH.h"
#include "InvertedIndex.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;

namespace po=boost::program_options;
namespace bf=boost::filesystem;
namespace bu=boost::numeric::ublas;
typedef bu::matrix<double> Matrix;


int fesco_inverted_index_cossim(int argc, char* argv[])
{

    try{
        int num_cores;
        string input;
        string output;
        bool verbose=false;

        
        po::options_description io_opts("Input and Output");
        io_opts.add_options()
            ("input,i",po::value<string>(&input),"Sketch file")
            ("output,o",po::value<string>(&output),"Similarity file");
   
        
        po::options_description pa_opts("Parameters");
        pa_opts.add_options()
            ("number-of-threads,t",po::value<int>(&num_cores)->default_value(1),"Specify the number of threads to be used");

        po::options_description misc_opts("Miscellaneous");
        misc_opts.add_options()
            ("verbose,v","Print runtime message")
            ("help,h","Print help message");

        // combine the above options
        po::options_description cmd_opts;
        cmd_opts.add(io_opts).add(pa_opts).add(misc_opts);

        po::positional_options_description p;
        p.add("input",-1);
        p.add("output",-1);

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmd_opts).positional(p).run(), vm);
        notify(vm);

        if (vm.count("help") || argc==1){
            cout << "Usage: fesco cossim [options] -i sketch -o similarity" << endl;
            cout << cmd_opts << endl;
            return 0;
        }

        if (!vm.count("input")){
            cout << "Please specify input." << endl;
            return 0;
        }else{
            bf::path abspath = bf::canonical(input);
            input = abspath.string();
        }

        if (!vm.count("output")){
            cout << "Please specify output." << endl;
            return 0;
        }

        if (vm.count("verbose")) verbose=true;

        // load inverted index
        runtime_message("Load sketch file",verbose);
        unordered_map<string,vector<string>> read_sketches;
        fesco_load_sketches(input,read_sketches);

        unsigned K = read_sketches.begin()->second.front().size();        
        unsigned L = read_sketches.begin()->second.size();
        unsigned N = read_sketches.size();
        double C = 3.141592653589793/(K*L);
        Matrix hamming(bu::zero_matrix<double>(N,N));

        runtime_message("Compute Hamming distance",verbose);

        #pragma omp parallel num_threads(num_cores) 
        {
            long i,j,l;
            #pragma omp for 
            for (i=0; i<N; ++i){
                for (j=i+1; j<N; ++j){
                    for (l=0; l<L; ++l){
                        auto iteri = next(read_sketches.begin(),i);
                        auto iterj = next(read_sketches.begin(),j);
                        double hd = HammingDistance(iteri->second[l],iterj->second[l]);
                        hamming(i,j) += hd*C;
                        hamming(j,i) += hd*C;
                    }
                }
            }
        }

        runtime_message("Compute Cosine similarity",verbose);
        #pragma omp parallel num_threads(num_cores) 
        {
            long i,j;
            #pragma omp for 
            for (i=0; i<N; ++i){
                for (j=i; j<N; ++j){
                    double cs=max(0.0,cos(hamming(i,j)));
                    if (j==i) hamming(i,j)=cs;
                    else{
                        hamming(i,j)=cs;
                        hamming(j,i)=cs;
                    }
                }
            }
        }

        runtime_message("Save results to "+output,verbose);

        unsigned W=0;
        for (auto p : read_sketches)
            W = max(W,unsigned(p.first.size()));

        ofstream ofs(output);
        #pragma omp parallel num_threads(num_cores) 
        {
            long i,j;
            #pragma omp for 
            for (i=0; i<N; ++i){
                for (j=i+1; j<N; ++j){
                    stringstream ss;
                    auto iteri = next(read_sketches.begin(),i);
                    auto iterj = next(read_sketches.begin(),j);
                    ss << iteri->first << "\t" 
                       << setw(W) << iterj->first << "\t"
                       << setw(7) << setprecision(7) << fixed << hamming(i,j) << "\n";
                    ofs << ss.str();
                    ofs.flush();
                }
            }
        }
        ofs.close();

        runtime_message("Done",verbose);
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;
}


