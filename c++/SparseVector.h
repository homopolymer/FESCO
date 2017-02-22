#ifndef SPARSE_VECTOR
#define SPARSE_VECTOR

#include <vector>
#include <sstream>
#include <ostream>
#include <algorithm>
using namespace std;

class SparseVector
{
public:
    SparseVector()
        :d(0),nnz(0)
    {
        allocate();
    }
    SparseVector(long _d, long _nnz)
        :d(_d),nnz(_nnz)
    {
        allocate();
    }

    SparseVector(const SparseVector& other)
        :d(other.d),nnz(other.nnz)
    {
        allocate();
        copy(other.data,other.data+other.nnz,data);
        copy(other.ri,other.ri+other.nnz,ri);
    }

    ~SparseVector()
    {
        clean();
    }

    
public:
    void clean()
    {
        if (data!=nullptr) delete[] data;
        if (ri!=nullptr) delete[] ri;
        data=nullptr;
        ri=nullptr;
    }
    void resize(long _d, long _nnz)
    {
        d = _d;
        nnz = _nnz;

        clean();
        allocate();
    }

    
public:
    double norm();
    SparseVector& operator-=(const SparseVector& other);
    SparseVector& operator*=(double scale);
    SparseVector& operator/=(double scale);
    SparseVector& operator=(const SparseVector& other);

private:
    void allocate()
    {
        data=new double[nnz];
        ri=new long[nnz];
    }

public:
    double *data;
    long   *ri;
    long    d;
    long    nnz;
};

ostream& operator<<(ostream& stream, const SparseVector& v);
SparseVector operator*(const SparseVector& v, double s);
SparseVector operator*(double s, const SparseVector& v);
SparseVector operator/(const SparseVector& v, double s);
SparseVector operator-(const SparseVector& u, const SparseVector& v);
double inner_dot(const SparseVector& a, const SparseVector& b);

#endif
