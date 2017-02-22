// routine to build inverted index for large-scale data
// blocked processing 

#include "Fesco.h"
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include "InvertedIndex.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
using namespace std;

namespace po=boost::program_options;
namespace bf=boost::filesystem;


int fesco_inverted_index_stats(int argc, char* argv[])
{

    try{
        int num_cores;
        string input;
        string output;
        bool head=false;
        bool compressed=false;
        bool verbose=false;

        
        po::options_description io_opts("Input and Output");
        io_opts.add_options()
            ("input,i",po::value<string>(&input),"Inverted index file")
            ("output,o",po::value<string>(&output),"Statistics file");
   
        
        po::options_description pa_opts("Parameters");
        pa_opts.add_options()
            ("label,l","Print table head")
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
            cout << "Usage: fesco stats [options] -i index [-o stats]" << endl;
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


        if (vm.count("label")) head=true;
        if (vm.count("compress")) compressed=true;
        if (vm.count("verbose")) verbose=true;

        // load inverted index
        runtime_message("Load the inverted index file",verbose);
        InvertedIndex indexer;
        indexer.Load(input,compressed);

        runtime_message("There are "+to_string(indexer.n_kmers)+" kmers and "+
                        to_string(indexer.nnz)+" non-zero values",verbose);
        
        // internal data structure
        map<double,double,less<double>> df;
        map<double,double,less<double>> cf;

        for (auto p : indexer.inverted_index){
            double t_df = p.second.front().freq;
            double t_cf = p.second.front().cf;
            if (df.count(t_df)==0) df[t_df]=1;
            else df[t_df]+=1;
            if (cf.count(t_cf)==0) cf[t_cf]=1;
            else cf[t_cf]+=1;
        }

        // lambda function for cumulative distrubition function
        auto cdf = [](map<double,double,less<double>>& a){
            vector<double> _cdf(a.size()+1,0);
            long i=0;
            for (auto p : a){
                i++;
                _cdf[i]=_cdf[i-1]+p.second;
            }
            long z=_cdf.back();
            for (auto& p : _cdf){
                p/=z;
            }
            return _cdf;
        };

        // lambda function for quantile index
        auto quantile = [](vector<double>& a, double q){
            for (unsigned long i=1; i<a.size(); ++i)
                if (a[i]>=q && a[i-1]<=q) 
                    return i;             
            return (unsigned long)0;
        };

        // lambda function for mean value
        auto mean = [](map<double,double,less<double>>& a){
             double z = 0;
             double x = 0;
             #pragma omp parallel for shared(x,z)
             for (unsigned i=0; i<a.size(); ++i){
                 auto iter = next(a.begin(),i);
                 #pragma omp atomic
                 x += iter->first*iter->second;
                 #pragma omp atomic
                 z += iter->second;
             }
             return x/z;
        };
        
        // for df
        runtime_message("Computing statistics for df",verbose);
        vector<double> df_cdf = cdf(df);
        unsigned long df_q5i = quantile(df_cdf,0.05);
        unsigned long df_q25i = quantile(df_cdf,0.25);
        unsigned long df_q50i = quantile(df_cdf,0.5);
        unsigned long df_q75i = quantile(df_cdf,0.75);
        unsigned long df_q95i = quantile(df_cdf,0.95);
        double df_mean = mean(df);

        // for cf
        runtime_message("Computing statistics for cf",verbose);
        vector<double> cf_cdf = cdf(cf);
        unsigned long cf_q5i = quantile(cf_cdf,0.05);
        unsigned long cf_q25i = quantile(cf_cdf,0.25);
        unsigned long cf_q50i = quantile(cf_cdf,0.5);
        unsigned long cf_q75i = quantile(cf_cdf,0.75);
        unsigned long cf_q95i = quantile(cf_cdf,0.95);
        double cf_mean = mean(cf);

        // report statistics
        cout << setw(10) << "Statistics" << " " << setw(10) << "Mean" <<  " " 
             << setw(10) << "q5" <<  " " 
             << setw(10) << "q25" <<  " " 
             << setw(10) << "q50" << " " 
             << setw(10) << "q75" <<  " " 
             << setw(10) << "q95" << endl;
        cout << setw(10) << "df" << " " << setw(10) << df_mean << " "
             << setw(10) << next(df.begin(),df_q5i-1)->first << " "
             << setw(10) << next(df.begin(),df_q25i-1)->first << " "
             << setw(10) << next(df.begin(),df_q50i-1)->first << " "
             << setw(10) << next(df.begin(),df_q75i-1)->first << " "
             << setw(10) << next(df.begin(),df_q95i-1)->first << endl;
        cout << setw(10) << "cf" << " " << setw(10) << cf_mean << " "
             << setw(10) << next(cf.begin(),cf_q5i-1)->first << " "
             << setw(10) << next(cf.begin(),cf_q25i-1)->first << " "
             << setw(10) << next(cf.begin(),cf_q50i-1)->first << " "
             << setw(10) << next(cf.begin(),cf_q75i-1)->first << " "
             << setw(10) << next(cf.begin(),cf_q95i-1)->first << endl;

        // output term statistics
        if (vm.count("output"))
        {
            runtime_message("Savint k-mer info to "+output,verbose);
        
            ofstream ofs(output);

            if (head){
                stringstream ss;
                ss << "K-MER" << "\t" << setw(15) << "DF" << "\t"
                   << setw(15) << "IDF" << "\t"
                   << setw(15) << "CF" << "\n";
                ofs << ss.str();
                ofs.flush();
            }

            #pragma omp parallel num_threads(num_cores)
            {
                #pragma omp for 
                for (unsigned long i=0; i<indexer.inverted_index.size(); ++i){
                    auto iter=next(indexer.inverted_index.begin(),i);
                    stringstream ss;
                    ss << iter->first << "\t" << setw(15) << iter->second.front().freq << "\t"
                       << setw(15) << iter->second.front().idf << "\t"
                       << setw(15) << iter->second.front().cf << "\n";
                    ofs << ss.str();
                    ofs.flush();
                }
            }

            ofs.close();
        }
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;
}
