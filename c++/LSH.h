#ifndef MY_LSH
#define MY_LSH

#include <map>
#include <list>
#include <vector>
#include <string>
#include <random>
#include "SparseVector.h"
using namespace std;

class LSH
{
public:
    LSH(long _d,long _k,bool _sparse)
        :d(_d),k(_k),sparse(_sparse)
    {
        P=new SparseVector[k];
        // random setting
        seed_seq seed{rd(),rd(),rd(),rd(),rd(),rd()};
        rng.seed(seed);
        if (sparse) GenerateSparseRandomProjection();
        else GenerateRandomProjection();
        GramSchmidtProcess();
    }
    LSH(const LSH& other)
        :d(other.d),k(other.k),sparse(other.sparse)
    {
        P=new SparseVector[k];
        for (long i=0; i<k; ++i)
            P[i]=other.P[i];
    }
    ~LSH(){clean();}

private:
    void clean()
    {
        if (P!=nullptr) delete[] P;
        P=nullptr;
    }

    void GramSchmidtProcess();
    void GenerateRandomProjection();
    void GenerateSparseRandomProjection();
    void Fingerprint(SparseVector& x, string& c);

public:
    void Bitcode(SparseVector& x, string& b);
    string Bitcode(SparseVector& x);

public:
    long d; // dimension of original space
    long k; // dimension of projected space
    bool sparse;
    SparseVector *P; // projection matrix

    random_device rd;
    mt19937 rng; 
};

unsigned HammingDistance(const string& bitcode1, const string& bitcode2);

#endif
