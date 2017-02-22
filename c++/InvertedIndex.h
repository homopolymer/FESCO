#ifndef MY_INVERTED_INDEX
#define MY_INVERTED_INDEX

#include <unordered_map>
#include <list>
#include <string>
#include "SparseMatrix.h"
using namespace std;

// using boost serialization for read and save
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/list.hpp>



class InvertedIndex
{
public:
    InvertedIndex(){nnz=0;n_kmers=0;n_reads=0;read_id_0=0;sorted=false;}
    InvertedIndex(long _read_id_0)
        : read_id_0(_read_id_0)
    {nnz=0;n_kmers=0;n_reads=0;sorted=false;}
    InvertedIndex(const InvertedIndex& other)
        : nnz(other.nnz),
          n_kmers(other.n_kmers),
          n_reads(other.n_reads),
          read_id_0(other.read_id_0),
          sorted(other.sorted)
    {
         kmers=other.kmers;
         read_pool=other.read_pool;
         read_kmers=other.read_kmers;
         read_max_kmer=other.read_max_kmer;
         for (auto p: other.inverted_index) {inverted_index[p.first]=p.second;}
    }
    ~InvertedIndex(){}

public:
    // insert a <kmer, read> pair to the inverted index
    void Insert(const string& kmer, const string& read);

    // sort kmers
    void Sort();

    // clean data
    void Clean();

    // calculate tf-idf
    void CalculateTfIdf();

    // remove some kmers, of which frequency is lower than threshold
    void FilterLow(int thresh=2);
    // remove some kmers, of which frequency is higher than threshold
    void FilterHigh(int thresh);

    // convert the inverted index to a compressed column storage sparse matrix
    void TransToCsc(SparseMatrix& csc);

    // save to a local file
    void Save(const string& filename, bool compressed=true);
    // load from a local file
    void Load(const string& filename, bool compressed=true);

    // merge other inverted index
    void Merge(const InvertedIndex& other);

public:
    // an item of a posting
    struct PostingItem
    {
        PostingItem():id(-1),abs_id(-1),freq(0.0),cf(0.0),idf(0.0),tfidf(0.0){}
        PostingItem(long _id, long double _freq):id(_id),abs_id(-1),freq(_freq),cf(0.0),idf(0.0),tfidf(0.0){}
        PostingItem(long _id, long _abs_id, long double _freq):id(_id),abs_id(_abs_id),freq(_freq),cf(0.0),idf(0.0),tfidf(0.0){}
        long id;     // relative id in local data
        long abs_id; // absolute id in whole data
        long double freq;
        long double cf;
        long double idf;
        long double tfidf;


        // serialize
        template<class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & id;
            ar & abs_id;
            ar & freq;
            ar & cf;
            ar & idf;
            ar & tfidf;
        }
    };

    typedef list<PostingItem> Posting;

public:
    unordered_map<string,Posting> inverted_index; // inverted indices
    list<string> read_pool; // read counter
    unordered_map<long,long> read_kmers; // read kmer counter
    unordered_map<long,long> read_max_kmer; // read most-frequent kmer count
    list<string> kmers; // all kmers
    long nnz; // number of non-zeros
    long n_kmers;
    long n_reads;
    long read_id_0; // begin number of reads
    bool sorted;

    typedef unordered_map<string,Posting>::const_iterator inverted_index_iterator;
public:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & nnz;
        ar & n_kmers;
        ar & n_reads;
        ar & read_id_0;
        ar & kmers;
        ar & read_pool;
        ar & read_kmers;
        ar & read_max_kmer;
        ar & inverted_index;
    }
};

#endif
