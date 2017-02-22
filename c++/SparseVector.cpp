#include <cmath>
#include "SparseVector.h"
using namespace std;

double SparseVector::norm()
{
    return sqrt(inner_dot(*this,*this));
}

SparseVector& SparseVector::operator-=(const SparseVector& other)
{
   *this = *this - other; 
   return *this;
}

SparseVector& SparseVector::operator*=(double scale)
{
    *this = (*this)*scale;
    return *this;
}

SparseVector& SparseVector::operator/=(double scale)
{
    *this = (*this)/scale;
    return *this;
}

SparseVector& SparseVector::operator=(const SparseVector& other)
{
    if (this!=&other){
        this->nnz=other.nnz;
        this->d=other.d;
        this->resize(d,nnz);
        copy(other.ri,other.ri+other.nnz,this->ri);
        copy(other.data,other.data+other.nnz,this->data);
    }
    return *this;
}

ostream& operator<<(ostream& out, const SparseVector& x)
{
    stringstream ss;
    ss << "(" << x.d << "," << x.nnz << "): ";
    ss << "[";
    ss << "(" << *x.ri << "," << *x.data << ")";
    for (long i=1; i<x.nnz; ++i){
        ss << "," << "(" << *(x.ri+i) << "," << *(x.data+i) << ")";
    }
    ss << "]"; 
    out << ss.str();

    return out;
}

SparseVector operator*(const SparseVector& x, double s)
{
    SparseVector y(x.d,x.nnz);
    for (long i=0; i<x.nnz; i++){
        y.ri[i] = x.ri[i];
        y.data[i] = s*x.data[i];
    }
    return y;
}

SparseVector operator*(double s, const SparseVector& x)
{
    SparseVector y(x.d,x.nnz);
    for (long i=0; i<x.nnz; i++){
        y.ri[i] = x.ri[i];
        y.data[i] = s*x.data[i];
    }
    return y;
}

SparseVector operator/(const SparseVector& x, double s)
{
    SparseVector y(x.d,x.nnz);
    for (long i=0; i<x.nnz; i++){
        y.ri[i] = x.ri[i];
        y.data[i] = x.data[i]/s;
    }
    return y;
}

SparseVector operator-(const SparseVector& u, const SparseVector& v)
{
    long *ri = new long[u.d];
    double *data = new double[u.d];
    long nnz = 0;
    long d = u.d;

    long i,j,k=0,x=0,y=0;
    while (x<u.nnz && y<v.nnz){
        i = u.ri[x];
        j = v.ri[y];
        if (i==j){
            ri[k] = i;
            data[k] = u.data[x] - v.data[y];
            k++;
            x++;
            y++;
        }else if (i<j){
            ri[k] = i;
            data[k] = u.data[x];
            k++;
            x++;
        }else if (i>j){
            ri[k] = j;
            data[k] = -v.data[y];
            k++;
            y++;
        }
    }

    for(;x<u.nnz;++x,++k){
        ri[k] = u.ri[x];
        data[k] = u.data[x];
    }

    for(;y<v.nnz;++y,++k){
        ri[k] = v.ri[y];
        data[k] = -v.data[y];
    }

    nnz = k;
    SparseVector z(d,nnz);
    copy(data,data+nnz,z.data);
    copy(ri,ri+nnz,z.ri);

    delete [] data;
    delete [] ri;
    data = nullptr;
    ri = nullptr;

    return z;
}

double inner_dot(const SparseVector& a, const SparseVector& b)
{
    long i,j,x=0,y=0;
    double val=0;

    while (x<a.nnz && y<b.nnz){
        i=a.ri[x];
        j=b.ri[y];
        if (i==j){
            val+=a.data[x]*b.data[y];
            x++;
            y++;
        }else if (i<j){
            x++;
        }else if (i>j){
            y++;
        }
    }

    return val;
}

