#include "InvertedIndex.h"
#include <cmath>
#include <iterator>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace bar=boost::archive;
namespace bio=boost::iostreams;

void InvertedIndex::Insert(const string& kmer, const string& read)
{
    long read_id;
    // assign id to read
    if (this->read_pool.empty()){
        read_id = 0;
        read_pool.push_back(read);
        this->n_reads+=1;
    }else if (this->read_pool.back()==read){
        read_id = this->read_pool.size()-1;
    }else{
        read_id = this->read_pool.size();
        read_pool.push_back(read);
        this->n_reads+=1;
    }

    // find a posting indicated by kmer
    auto iter = this->inverted_index.find(kmer);
    
    // a new kmer
    if (iter==this->inverted_index.end()){
        this->inverted_index[kmer] = Posting(1,PostingItem());
        this->n_kmers+=1;

        // record new kmer
        this->kmers.push_back(kmer);
    }   
   
    // add the pair <kmer, read> to the posting
    iter = this->inverted_index.find(kmer);

    if (iter->second.front().id==read_id){
        // update tf
        iter->second.back().freq += 1.0;
        // update cf
        iter->second.front().cf += 1.0;

        // update most-frequent kmer count
        if (this->read_max_kmer[read_id]<iter->second.back().freq)
            this->read_max_kmer[read_id] = iter->second.back().freq;
    }else{
        iter->second.push_back(PostingItem(read_id,read_id+this->read_id_0,1.0));
        iter->second.front().id=read_id;
        iter->second.front().freq+=1.0;
        iter->second.front().cf+=1.0;
        // count non-zero elements
        this->nnz+=1;

        // update read-kmer counter
        auto riter=this->read_kmers.find(read_id);
        if (riter==this->read_kmers.end()){
            this->read_kmers[read_id]=1;

            // update most-frequent kmer count
            this->read_max_kmer[read_id]=1;
        }else{
            riter->second+=1;

            // update most-frequent kmer count
            if (this->read_max_kmer[read_id]<iter->second.back().freq)
                this->read_max_kmer[read_id] = iter->second.back().freq;
        }
    }

}


void InvertedIndex::CalculateTfIdf()
{
    long double eps=0;
    long double N = this->read_pool.size();

    for (auto kiter=this->inverted_index.begin(); kiter!=this->inverted_index.end(); ++kiter){
        auto piter=kiter->second.begin();

        // idf
        piter->idf=log(N/piter->freq)+eps;       

        // tf-idf
        for(++piter; piter!=kiter->second.end(); ++piter){
            piter->tfidf=piter->freq*kiter->second.front().idf;
        }
    }
}


void InvertedIndex::FilterLow(int thresh)
{
    long new_id=0;
    for (auto iter=inverted_index.begin(); iter!=inverted_index.end(); ){
        if (iter->second.front().freq<thresh){
            iter = inverted_index.erase(iter);
        }else{
            if (iter->second.front().id!=new_id){
                iter->second.begin()->id=new_id;
            }
            ++iter;
            ++new_id;
        }
    }

    n_kmers = inverted_index.size();
}


void InvertedIndex::FilterHigh(int thresh)
{
    long new_id=0;
    for (auto iter=inverted_index.begin(); iter!=inverted_index.end(); ){
        if (iter->second.front().freq>thresh){
            iter = inverted_index.erase(iter);
        }else{
            if (iter->second.front().id!=new_id){
                iter->second.begin()->id=new_id;
            }
            ++iter;
            ++new_id;
        }
    }

    n_kmers = inverted_index.size();
}

void InvertedIndex::TransToCsc(SparseMatrix& csc)
{
    // resize
    csc.resize(n_kmers,n_reads,nnz);
    csc.cj[0] = 0;

    long i,j,nnzj,nni;
    for (j=0,nni=0; j<n_reads; ++j){
        i = 0;
        nnzj = 0;
        for (auto iter=inverted_index.begin(); iter!=inverted_index.end(); ++iter,++i){
            auto riter = iter->second.begin();
            for (++riter; riter!=iter->second.end(); ++riter){

                if (riter->id>j) break;

                if (riter->id==j){
                    csc.data[nni]=riter->tfidf;
                    csc.ri[nni]=i;
                    nni++;
                    nnzj++;
                }
            }
        }
        csc.cj[j+1]=csc.cj[j]+nnzj;
    }
}

void InvertedIndex::Sort()
{
    this->kmers.sort();
    this->sorted = true;
}


void InvertedIndex::Save(const string& filename, bool compressed)
{
    ofstream ofs(filename);
    
    bio::filtering_stream<bio::output> f;
    if (compressed) f.push(bio::gzip_compressor());
    f.push(ofs);
    bar::binary_oarchive oa(f);
    //*/

    oa << (*this);
}


void InvertedIndex::Load(const string& filename, bool compressed)
{
    ifstream ifs(filename);

    bio::filtering_stream<bio::input> f;
    if (compressed) f.push(bio::gzip_decompressor());
    f.push(ifs);
    bar::binary_iarchive ia(f);
    //*/

    ia >> (*this);
}

void InvertedIndex::Merge(const InvertedIndex& other)
{
    long read_shift = this->n_reads;

    // update number of non-zeros
    this->nnz += other.nnz;
    // update number of reads
    this->n_reads += other.n_reads;
    // update read pool
    copy(other.read_pool.begin(), other.read_pool.end(), back_inserter(this->read_pool));
    // update read kmer count
    for (auto p : other.read_kmers){
        this->read_kmers[p.first+read_shift] = p.second;
    }
    // update read max kmer count
    for (auto p : other.read_max_kmer){
        this->read_max_kmer[p.first+read_shift] = p.second;
    }

    // update inverted index and kmers

    auto this_kmer_iter = this->kmers.begin();
    auto other_kmer_iter = other.kmers.begin();
    
    while(this_kmer_iter!=this->kmers.end() &&
          other_kmer_iter!=other.kmers.end()){

        string this_kmer = *this_kmer_iter;
        string other_kmer = *other_kmer_iter;

        if (this_kmer==other_kmer){
            auto this_post_iter = this->inverted_index[this_kmer].begin();

            auto other_kmer_post_iter = other.inverted_index.find(other_kmer);
            auto other_post_iter = other_kmer_post_iter->second.begin();
            auto other_post_end = other_kmer_post_iter->second.end();

            long this_post_size = this->inverted_index[this_kmer].size();

            // update kmer frequency
            this_post_iter->freq += other_post_iter->freq;
            // update collection frequency
            this_post_iter->cf += other_post_iter->cf;
            
            // insert alien posting to the end of current posting
            copy(++other_post_iter,
                 other_post_end,
                 back_inserter(this->inverted_index[this_kmer]));

            // modify read information
            advance(this_post_iter, this_post_size);
            for (;this_post_iter!=this->inverted_index[this_kmer].end(); ++this_post_iter){
                this_post_iter->id += read_shift;
            }

            // udpate iterator
            this_kmer_iter++;
            other_kmer_iter++;

        }else if (this_kmer>other_kmer){
            // insert alien kmer
            this->kmers.insert(this_kmer_iter,other_kmer);
            
            // update inverted index
            this->inverted_index[other_kmer] = other.inverted_index.find(other_kmer)->second;

            // modify read information
            auto post_iter = this->inverted_index[other_kmer].begin();
            for (++post_iter; post_iter!=this->inverted_index[other_kmer].end(); ++post_iter){
                post_iter->id += read_shift;
            }

            // update iterator
            other_kmer_iter++;

        }else{
            this_kmer_iter++;

        }
     
    }

    for (; other_kmer_iter!=other.kmers.end(); ++other_kmer_iter){
        string other_kmer = *other_kmer_iter;
        // insert alien kmer
        this->kmers.insert(this_kmer_iter,other_kmer);
        
        // update inverted index
        this->inverted_index[other_kmer] = other.inverted_index.find(other_kmer)->second;
        
        // modify read information
        auto post_iter = this->inverted_index[other_kmer].begin();
        for (++post_iter; post_iter!=this->inverted_index[other_kmer].end(); ++post_iter){
            post_iter->id += read_shift;
        }
    }

    this->n_kmers = this->kmers.size();
}

void InvertedIndex::Clean()
{
    if (!kmers.empty()) kmers.erase(kmers.begin(),kmers.end());
    if (!read_pool.empty()) read_pool.erase(read_pool.begin(),read_pool.end());
    if (!read_kmers.empty()) read_kmers.erase(read_kmers.begin(),read_kmers.end());
    if (!read_max_kmer.empty()) read_max_kmer.erase(read_max_kmer.begin(),read_max_kmer.end());
    if (!inverted_index.empty()){
        for (auto p : inverted_index)
            if (!p.second.empty())
                p.second.erase(p.second.begin(), p.second.end());

        inverted_index.erase(inverted_index.begin(),inverted_index.end());
    }
}
