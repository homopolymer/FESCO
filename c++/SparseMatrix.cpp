#include "SparseMatrix.h"

void SparseMatrix::get_column(long j, SparseVector& mj)
{
    long dj = m;
    long nnzj = cj[j+1]-cj[j];

    // resize 
    mj.resize(dj,nnzj);

    // retrieval
    for (long i=cj[j]; i<cj[j+1]; ++i){
        mj.data[i-cj[j]] = data[i];
        mj.ri[i-cj[j]] = ri[i];
    }
}
