#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <vector>
#include "SparseVector.h"
using namespace std;

// compressed column storage
class SparseMatrix
{
public:
    SparseMatrix() {data=nullptr;ri=nullptr;cj=nullptr;}
    SparseMatrix(long _m, long _n, long _nnz)
    {
        resize(_m,_n,_nnz);
    }

    ~SparseMatrix()
    {
        clean();
    }


public:
    void clean()
    {
        if (data!=nullptr) delete[] data;
        if (ri!=nullptr) delete[] ri;
        if (cj!=nullptr) delete[] cj;
        data=nullptr;
        ri=nullptr;
        cj=nullptr;
    }

    void resize(long _m, long _n, long _nnz)
    {
        m = _m;
        n = _n;
        nnz = _nnz;

        clean();
        allocate();
    }

    void get_column(long j, SparseVector& mj);
    SparseVector operator[](long i)
    {
        SparseVector column;
        get_column(i,column);
        return column;
    }

private:
    void allocate()
    {
        data=new double[nnz];
        cj=new long[n+1];
        ri=new long[nnz];
    }

public:
    double *data;
    long *cj; // compressed column index 
    long *ri; // row index
    long m; // row dimension
    long n; // column dimension
    long nnz; // number of non-zeros
};

#endif
