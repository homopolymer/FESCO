#include "LSH.h"
#include <cmath>
#include <random>
#include <sstream>
#include <iostream>
#include <numeric>
#include <functional>
#include <boost/dynamic_bitset.hpp>
using namespace std;

void LSH::GenerateRandomProjection()
{
    normal_distribution<> dice(0,1);

    long i,j;
    double r;
    for (i=0; i<k; ++i){
        long _nnz = 0;

        P[i].d=d;
        P[i].data = new double[d];
        P[i].ri=new long[d];

        for (j=0; j<d; ++j){
            r=dice(rng);

            P[i].data[_nnz]=r;
            P[i].ri[_nnz]=j;
            _nnz++;

        }

        P[i].nnz = _nnz;
    }
}

void LSH::GenerateSparseRandomProjection()
{
    uniform_real_distribution<> dice(0,1);

    long i,j;
    double r,s=3,l=1/s,ro=0.5,t=sqrt(k);
    double a=sqrt(s/t);
    for (i=0; i<k; ++i){
        long _nnz = 0;

        P[i].d=d;
        P[i].data = new double[d];
        P[i].ri=new long[d];

        for (j=0; j<d; ++j){
            r=dice(rng);

            if (r<ro*l){
                P[i].data[_nnz]=-a;
                P[i].ri[_nnz]=j;
                _nnz++;
            }else if (r>1-l+ro*l){
                P[i].data[_nnz]=a;
                P[i].ri[_nnz]=j;
                _nnz++;
            }//*/

        }

        P[i].nnz = _nnz;
    }
}

void LSH::GramSchmidtProcess()
{
    long i,j,n=min(k,d),l=ceil(k/n),t;
    SparseVector *Q = new SparseVector[k];
    

    for (t=0; t<l; t++){
        if (t+1==l) n = min(n,k-t*n);

        Q[t*n] = P[t*n]/P[t*n].norm();

        for (i=1; i<n; ++i){
            Q[t*n+i] = P[t*n+i];
            for (j=0; j<i; ++j){
                Q[t*n+i] = Q[t*n+i] - inner_dot(Q[t*n+i],Q[t*n+j])*Q[t*n+j];
            }       
            Q[t*n+i] /= Q[t*n+i].norm();
        }
    }

    for (i=0; i<k; ++i)
        P[i] = Q[i];

    delete [] Q;
    Q = nullptr;
}

void LSH::Fingerprint(SparseVector& x, string& c)
{
    stringstream buff;
    long double z=0;

    for (long i=0; i<k; ++i){
        // projection
        z=inner_dot(P[i],x);                 
        buff << ((z>0)?1:0);
    }

    c = buff.str();
}

void LSH::Bitcode(SparseVector& x, string& c)
{
    Fingerprint(x,c);
}

string LSH::Bitcode(SparseVector& x)
{
    string c;
    Bitcode(x,c);
    return c;
}

unsigned HammingDistance(const string& bitcode1, const string& bitcode2)
{
    boost::dynamic_bitset<> b1(bitcode1);
    boost::dynamic_bitset<> b2(bitcode2);   

    b1^=b2;

    return b1.count();
}
