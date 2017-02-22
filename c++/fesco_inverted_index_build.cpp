// routine to build inverted index for large-scale data
// blocked processing 

#include "Fesco.h"
#include <ctime>
#include <iomanip>
#include <zlib.h>
#include "kseq.h"
#include <list>
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include "InvertedIndex.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
using namespace std;


KSEQ_INIT(gzFile, gzread)


namespace po=boost::program_options;
namespace bf=boost::filesystem;

struct Block
{
    Block(){}
    Block(long _id, long _pos, long _len)
        : id(_id),pos(_pos),len(_len) 
    {}
    Block(const Block& other)
        : id(other.id),pos(other.pos),len(other.len)
    {}
    ~Block(){}
    long id;
    long pos;
    long len;
};

struct Kmer
{
    Kmer():pos(-1),str(""){}
    Kmer(const Kmer& other):pos(other.pos),str(other.str){}
    Kmer(long& _pos, string& _str):pos(_pos),str(_str){}
    ~Kmer(){}

    friend bool operator==(const Kmer& left, const Kmer& right)
    {
        return (left.pos==right.pos) && (left.str==right.str);
    }
    friend bool operator!=(const Kmer& left, const Kmer& right)
    {
        return (left.pos!=right.pos) || (left.str!=right.str);
    }

    long pos;
    string str;
};


void print_inverted_index(const InvertedIndex& indexer)
{
    for (auto p : indexer.kmers){
        auto piter = indexer.inverted_index.find(p);
        auto iter = piter->second.begin();
        auto end  = piter->second.cend();

        cout << "KMER [" << p << "," << iter->freq << "," << iter->cf << "," << iter->idf << "]: ";
        iter++;

        cout << "<" << iter->abs_id << "," << iter->freq << "," << iter->tfidf << ">";
        for (++iter;iter!=end;++iter){
            cout << "->" << "<" << iter->abs_id << "," 
                                << iter->freq << ","
                                << iter->tfidf
                         << ">";
        }
        cout << endl;
    }
}

/*
void runtime_message(const string& msg, bool verbose)
{
    time_t result = time(nullptr);
    char mbstr[128];
    strftime(mbstr,sizeof(mbstr),"%c %Z", localtime(&result));

    string msg2 = "[" + string(mbstr) + "]: " + msg + "\n";
    
    if (verbose) cerr << msg2;
}//*/


void blocked_inverted_index(const unordered_map<string,string>& reads, 
                            const vector<int>& kmers, 
                            const Block& block, 
                            const long minimizer_w,
                            InvertedIndex& indexer, 
                            bool verbose)
{
    runtime_message("Process block "+to_string(block.id)+", size="+to_string(block.len),verbose);

    unsigned long  w = minimizer_w;
    auto iter = next(reads.begin(),block.pos);

    auto cmp = [](const Kmer& left, const Kmer& right) {return left.str < right.str;};

    // modify read initial number
    indexer.read_id_0 = block.pos;

    // scan read-kmer pair
    for (long rc=0; rc<block.len; ++rc,++iter){
        for (auto k : kmers){
            long rl = iter->second.length();
            vector<Kmer> queue;
            Kmer minimizer;
            for (long i=0; i<rl-k+1; ++i){
                string kmer = iter->second.substr(i,k);
                if (queue.size()<w){
                    queue.push_back(Kmer(i,kmer));
                }else if (queue.size()==w){
                    queue.erase(queue.begin());
                    queue.push_back(Kmer(i,kmer));
                }else{
                    throw overflow_error("Error in blocked processing: overflow in queue");
                }

                if ((i+1)%w==0){
                    auto t_minimizer = min_element(queue.begin(),queue.end(),cmp);
                    if (*t_minimizer!=minimizer){
                        minimizer = *t_minimizer;
                        indexer.Insert(minimizer.str,iter->first);
                    }
                }
            }
        }
    }

    // sort
    indexer.Sort();

    runtime_message("Done block "+to_string(block.id),verbose);
}

int fesco_inverted_index_build(int argc, char* argv[])
{

    try{
        int num_blocks;
        int num_cores;
        int minimizer_w;
        string input;
        string output;
        bool compressed=false;
        bool verbose=false;
        vector<int> kmers;

        // to be developed in the future
        //
        // float tf_thres
        // float df_thres
        // float idf_thres
        // float tfidf_thres
        
        
        po::options_description io_opts("Input and Output");
        io_opts.add_options()
            ("input,i",po::value<string>(&input),"Sequencing reads (fasta or fastq)")
            ("output,o",po::value<string>(&output),"Inverted index");
   
        
        po::options_description pa_opts("Parameters");
        pa_opts.add_options()
            ("kmer-set,k",po::value<vector<int>>(&kmers)->multitoken()->default_value(vector<int>{11,21,31},"11,21,31"),
                          "Specify a set of k-mers for indexing")
            ("minimizer-w,w",po::value<int>(&minimizer_w)->default_value(1),"Specify the parameter w for minimizer")
            ("number-of-threads,t",po::value<int>(&num_cores)->default_value(1),"Specify the number of threads to be used")
            ("number-of-blocks,b",po::value<int>(&num_blocks)->default_value(1),"Specify the number of blocks for data division")
            ("compress,c","Save result to compressed file (gzip)");

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
            cout << "Usage: fesco build [options] -i sequences -o index" << endl;
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

        if (vm.count("compress")) compressed=true;
        if (vm.count("verbose")) verbose=true;

        // load sequencing reads
        runtime_message("Loading sequences",verbose);
        unordered_map<string,string> reads;
        long num_reads;
      
        {
            gzFile fp;
            kseq_t *seq;
            int l;
            
            fp = gzopen(input.c_str(),"r");
            seq = kseq_init(fp);
            while ((l=kseq_read(seq)) >=0){
                reads[seq->name.s] = seq->seq.s;
            }

            kseq_destroy(seq);
            gzclose(fp);
        }

        num_reads = reads.size();
        runtime_message("Loaded "+to_string(num_reads)+" sequences",verbose);

        // block division
        num_blocks = num_blocks<num_reads?num_blocks:num_reads;
        vector<Block> blocks(num_blocks);
        long block_size = floor(float(num_reads)/float(num_blocks));
        for (long i=0; i<num_blocks; ++i){
             long pos=i*block_size;
             long len=(i+1==num_blocks)?(num_reads-i*block_size):block_size;            
             blocks[i]=Block(i,pos,len);
        }

        // build blocked index
        runtime_message("Building blocked inverted index",verbose);
        vector<InvertedIndex> blocked_indices(num_blocks);

        num_cores = num_cores<num_blocks?num_cores:num_blocks;
        #pragma omp parallel num_threads(num_cores)
        {
            #pragma omp for 
            for (long i=0; i<num_blocks; ++i){
                blocked_inverted_index(reads, kmers, 
                                       blocks[i],
                                       minimizer_w,
                                       blocked_indices[i],
                                       verbose);
            
            }
        }

        // some statistics
        long nnz=0;
        for (auto p : blocked_indices){
            nnz+=p.nnz;
        }
        runtime_message(to_string(num_blocks)+" blocks contain "
                        +to_string(nnz)+" non-zero values",verbose);

        // release memory
        reads.erase(reads.begin(),reads.end());

        // merge blocked index files
        runtime_message("Merging blocked inverted indices",verbose);

        long ss=1,mi=0;
        while (ss<num_blocks){
             
            #pragma omp parallel num_threads(num_cores)
            {
                #pragma omp for 
                for (long i=0; i<num_blocks; i+=2*ss){
                    if (i+ss>=num_blocks) continue;
                    blocked_indices[i].Merge(blocked_indices[i+ss]);
                    blocked_indices[i+ss].Clean();
                }
            }

            ss*=2;
        }

        runtime_message("There are "+to_string(blocked_indices[mi].n_kmers)+" kmers",verbose);

        // calculate TF-IDF
        runtime_message("Computing TF-IDF",verbose);
        blocked_indices[mi].CalculateTfIdf();

        // output
        runtime_message("Save results to "+output,verbose);
        blocked_indices[mi].Save(output,compressed);

        runtime_message("Done!",verbose);
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;
}
