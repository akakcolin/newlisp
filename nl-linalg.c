#include "newlisp.h"
#include "linalq-protos.h"
#include "float.h"

/*
  dot
  diag
  diagonal
  tri
  eye
  vdot
  inner
  outer(a, b[, out])
  matmul
  tensordot
  einsum
  matrix_power
  kron
  lstsq
  lll_reduce
  niggli_reduce
  
 */
void lstsq(double *A, double *b, int M, int N) {
  // numpy.linalg.lstsq
  // A_dim1 == b_dim0
  int m = M, n = N, nrhs = 1, lda = M, ldb = N, rank = 0, info = 1, lwork = -1;
  int iwork[11 * m]; // 3*min(m,n)*n1v1 + 11*min(m,n)
  double rcond = -1;
  double wkopt = 0.0;
  double *work;
  double s[M];
  dgelsd_(&m, &n, &nrhs, A, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork,
          iwork, &info);
  lwork = (int)wkopt;
  work = (Real *)malloc(lwork * sizeof(Real));
  // solve A*x = b
  dgelsd_(&m, &n, &nrhs, A, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
          iwork, &info);
  free(work);
}
