// routine to build inverted index for large-scale data
// blocked processing 

#include "Fesco.h"
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
using namespace std;

namespace po=boost::program_options;
namespace bf=boost::filesystem;


int fesco_inverted_index_sketch(int argc, char* argv[])
{

    try{
        int num_cores;
        int K;
        int L;
        string input;
        string output;
        bool sparse=false;
        bool compressed=false;
        bool verbose=false;

        
        po::options_description io_opts("Input and Output");
        io_opts.add_options()
            ("input,i",po::value<string>(&input),"Inverted index file")
            ("output,o",po::value<string>(&output),"Sketch file");
   
        
        po::options_description pa_opts("Parameters");
        pa_opts.add_options()
            ("projection-dimension,d",po::value<int>(&K)->default_value(16),"Specify the dimension of projected space")
            ("projection-number,n",po::value<int>(&L)->default_value(100),"Specify the number of projections")
            ("projection-sparse,s","Declare to use sparse random projection")
            ("number-of-threads,t",po::value<int>(&num_cores)->default_value(1),"Specify the number of threads to be used")
            ("compress,c","Declare inverted index file is compressed");

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
            cout << "Usage: fesco sketch [options] -i index -o sketch" << endl;
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

        if (vm.count("projection-sparse")) sparse=true;
        if (vm.count("compress")) compressed=true;
        if (vm.count("verbose")) verbose=true;

        // load inverted index
        runtime_message("Load the inverted index file",verbose);
        InvertedIndex indexer;
        indexer.Load(input,compressed);

        runtime_message("There are "+to_string(indexer.n_kmers)+" kmers and "+
                        to_string(indexer.nnz)+" non-zero values",verbose);
        
        // LSH
        runtime_message("Build "+to_string(L)+" random projections",verbose);
        vector<LSH> LSHes;
        long i;
        for (i=0; i<L; ++i)
            LSHes.push_back(LSH(indexer.n_kmers,K,sparse)); 

        // get kmer-read tf-idf matrix
        SparseMatrix csc;
        indexer.TransToCsc(csc);

        // compute sketchs for reads
        runtime_message("Compute sketches for "+to_string(indexer.n_reads)+" reads",verbose);
        unordered_map<string,vector<string>> read_sketches;
        // initialize the map
        for (auto p : indexer.read_pool)
            read_sketches[p] = vector<string>(L);

        #pragma omp parallel num_threads(num_cores) 
        {
            #pragma omp for 
            for (i=0;i<indexer.n_reads;++i){
                // read name
                string read = *(next(indexer.read_pool.begin(),i));
                // read feature vector
                SparseVector features;
                csc.get_column(i,features);
                // read sketches
                for (long j=0; j<L; ++j){
                    read_sketches[read][j] = LSHes[j].Bitcode(features);
                }
            }
        }

        // save sketches to output
        runtime_message("Save results to "+output,verbose);
        fesco_save_sketches(output,read_sketches);

        runtime_message("Done",verbose);
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;
}


void fesco_save_sketches(const string& filename, const unordered_map<string,vector<string>>& read_sketches)
{
    unsigned w = 0;
    for (auto p : read_sketches){
        w = max(w,unsigned(p.first.length()));
    }

    ofstream ofs(filename,ios::binary);

    ofs << "@N" << " " << read_sketches.begin()->second.size() << endl;
  
    for (auto p : read_sketches){
        ofs << setw(w) << p.first;
    
        for (auto q : p.second){
            ofs << " " << setw(w) << q;
        }

        ofs << endl;
    }

    ofs.close();
}

void fesco_load_sketches(const string& filename, unordered_map<string,vector<string>>& read_sketches)
{
    ifstream ifs(filename,ios::binary);
    string line;
    string firstitem;
    unsigned L;

    while(getline(ifs,line)){
        istringstream iss(line);

        string firstitem;

        iss >> firstitem;
        if (firstitem=="@N"){
            iss >> L;
        }else{
            read_sketches[firstitem] = vector<string>(L);
            for (long i=0; i<L; ++i)
                iss >> read_sketches[firstitem][i];
        }
    }

   ifs.close();
}
