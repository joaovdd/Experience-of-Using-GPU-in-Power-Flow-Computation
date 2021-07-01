#pragma once
#include <cuComplex.h>
#include <cusparse.h>
#include "cusolverDn.h"
#include "cusolverSp.h"
#include "opcoesDeCompilacao.h"

//// defina como verdadeiro para que seja usado double,
//// caso contrário será utilizado float 
//// não esquecer de fazer mudança no outro arquivo também!!!!!!!!!!
//#define DOUBLE_MODE false 
//
//// *****************************************************************************************
//// Não alterar a partir daqui***************************************************************
//// *****************************************************************************************
//
//// definição do tipo de ponto flutuante
//#if DOUBLE_MODE
//typedef double float_type;
//#else
//typedef float float_type;
//#endif
//
//// definição do tipo de número complexo
//#if DOUBLE_MODE
//typedef cuDoubleComplex complex_type;
//#else
//typedef cuFloatComplex complex_type;
//#endif
//
//// definição do constritor complexo
//#if DOUBLE_MODE
//#define MK_COMPLEX_FUNCTION(x, y) make_cuDoubleComplex(x, y) 
//#else
//#define MK_COMPLEX_FUNCTION(x, y) make_cuFloatComplex(x, y)
//#endif

__host__ __device__ static inline cuDoubleComplex _cuAdd(cuDoubleComplex x, cuDoubleComplex y) { return cuCadd(x, y); }
__host__ __device__ static inline cuDoubleComplex _cuSub(cuDoubleComplex x, cuDoubleComplex y) { return cuCsub(x, y); }
__host__ __device__ static inline cuDoubleComplex _cuDiv(cuDoubleComplex x, cuDoubleComplex y) { return cuCdiv(x, y); }
__host__ __device__ static inline cuDoubleComplex _cuMul(cuDoubleComplex x, cuDoubleComplex y) { return cuCmul(x, y); }
__host__ __device__ static inline cuDoubleComplex _cuCon(cuDoubleComplex x) { return cuConj(x); }
__host__ __device__ static inline double _cuAbs(cuDoubleComplex x) { return cuCabs(x); }
__host__ __device__ static inline double _cuReal(cuDoubleComplex x) { return cuCreal(x); }
__host__ __device__ static inline double _cuImag(cuDoubleComplex x) { return cuCimag(x); }

// static inline cuDoubleComplex _mkComplex(double x, double y) { return make_cuDoubleComplex(x, y); }

__host__ __device__ static inline cuFloatComplex _cuAdd(cuFloatComplex x, cuFloatComplex y) { return cuCaddf(x, y); }
__host__ __device__ static inline cuFloatComplex _cuSub(cuFloatComplex x, cuFloatComplex y) { return cuCsubf(x, y); }
__host__ __device__ static inline cuFloatComplex _cuDiv(cuFloatComplex x, cuFloatComplex y) { return cuCdivf(x, y); }
__host__ __device__ static inline cuFloatComplex _cuMul(cuFloatComplex x, cuFloatComplex y) { return cuCmulf(x, y); }
__host__ __device__ static inline cuFloatComplex _cuCon(cuFloatComplex x) { return cuConjf(x); }
__host__ __device__ static inline float _cuAbs(cuFloatComplex x) { return cuCabsf(x); }
__host__ __device__ static inline float _cuReal(cuFloatComplex x) { return cuCrealf(x); }
__host__ __device__ static inline float _cuImag(cuFloatComplex x) { return cuCimagf(x); }

__host__ __device__ static inline complex_type _mkComplex(float_type x, float_type y) { return MK_COMPLEX_FUNCTION(x, y); }


// funções das bibliotecas CUDA


#if DOUBLE_MODE
__host__ __device__ static inline cusparseStatus_t _cusparseNnz(cusparseHandle_t          handle,
                                                                  cusparseDirection_t      dirA,
                                                                  int                      m,
                                                                  int                      n,
                                                                  const cusparseMatDescr_t descrA,
                                                                  const cuDoubleComplex*   A,
                                                                  int                      lda,
                                                                  int*                     nnzPerRowColumn,
                                                                  int*                     nnzTotalDevHostPtr) { 
                                                                      return cusparseZnnz(handle, dirA, m, n, descrA, A, lda, nnzPerRowColumn, nnzTotalDevHostPtr);
                                                                  }
#else
__host__ __device__ static inline cusparseStatus_t _cusparseNnz(cusparseHandle_t          handle,
                                                                  cusparseDirection_t      dirA,
                                                                  int                      m,
                                                                  int                      n,
                                                                  const cusparseMatDescr_t descrA,
                                                                  const cuFloatComplex*    A,
                                                                  int                      lda,
                                                                  int*                     nnzPerRowColumn,
                                                                  int*                     nnzTotalDevHostPtr) { 
                                                                      return cusparseCnnz(handle, dirA, m, n, descrA, A, lda, nnzPerRowColumn, nnzTotalDevHostPtr);
                                                                  }

#endif

#if DOUBLE_MODE
static inline cusparseStatus_t
_cusparseDense2csr(
    cusparseHandle_t         handle,
    int                      m,
    int                      n,
    const cusparseMatDescr_t descrA,
    const cuDoubleComplex* A,
    int                      lda,
    const int* nnzPerRow,
    cuDoubleComplex* csrValA,
    int* csrRowPtrA,
    int* csrColIndA) {
    return cusparseZdense2csr(handle, m, n, descrA, A, lda, nnzPerRow, csrValA, csrRowPtrA, csrColIndA);
}
#else
static inline cusparseStatus_t
_cusparseDense2csr(
    cusparseHandle_t         handle,
    int                      m,
    int                      n,
    const cusparseMatDescr_t descrA,
    const cuFloatComplex* A,
    int                      lda,
    const int* nnzPerRow,
    cuFloatComplex* csrValA,
    int* csrRowPtrA,
    int* csrColIndA) {
    return cusparseCdense2csr(handle, m, n, descrA, A, lda, nnzPerRow, csrValA, csrRowPtrA, csrColIndA);
}
#endif

#if DOUBLE_MODE
static inline cusolverStatus_t
_cusolverDnGetrf_bufferSize(
    cusolverDnHandle_t handle,
    int m,
    int n,
    double* A,
    int lda,
    int* Lwork) {
    return cusolverDnDgetrf_bufferSize(handle, m, n, A, lda, Lwork);
}
#else
static inline cusolverStatus_t
_cusolverDnGetrf_bufferSize(
    cusolverDnHandle_t handle,
    int m,
    int n,
    float* A,
    int lda,
    int* Lwork) {
    return cusolverDnSgetrf_bufferSize(handle, m, n, A, lda, Lwork);
}
#endif


#if DOUBLE_MODE
static inline cusolverStatus_t
_cusolverSpCsrlsvqr(cusolverSpHandle_t handle,
                    int m,
                    int nnz,
                    const cusparseMatDescr_t descrA,
                    const double *csrValA,
                    const int *csrRowPtrA,
                    const int *csrColIndA,
                    const double *b,
                    double tol,
                    int reorder,
                    double *x,
                    int *singularity) {
                        return cusolverSpDcsrlsvqr(handle, m, nnz, descrA, csrValA, csrRowPtrA, csrColIndA, b, tol, reorder, x, singularity);
                    }
#else
static inline cusolverStatus_t
_cusolverSpCsrlsvqr(cusolverSpHandle_t handle,
                    int m,
                    int nnz,
                    const cusparseMatDescr_t descrA,
                    const float *csrValA,
                    const int *csrRowPtrA,
                    const int *csrColIndA,
                    const float *b,
                    float tol,
                    int reorder,
                    float *x,
                    int *singularity) {
                        return cusolverSpScsrlsvqr(handle, m, nnz, descrA, csrValA, csrRowPtrA, csrColIndA, b, tol, reorder, x, singularity);
                    }
#endif

#if DOUBLE_MODE
cusparseStatus_t
_cusparseCsr2dense(cusparseHandle_t         handle,
                   int                      m,
                   int                      n,
                   const cusparseMatDescr_t descrA,
                   const double*            csrValA,
                   const int*               csrRowPtrA,
                   const int*               csrColIndA,
                   double*                  A,
                   int                      lda) {
                       return cusparseDcsr2dense(handle, m, n, descrA, csrValA, csrRowPtrA, csrColIndA, A, lda);
                   }

#else
cusparseStatus_t
_cusparseCsr2dense(cusparseHandle_t         handle,
                   int                      m,
                   int                      n,
                   const cusparseMatDescr_t descrA,
                   const float*            csrValA,
                   const int*               csrRowPtrA,
                   const int*               csrColIndA,
                   float*                  A,
                   int                      lda) {
                       return cusparseScsr2dense(handle, m, n, descrA, csrValA, csrRowPtrA, csrColIndA, A, lda);
                   }
#endif

#if DOUBLE_MODE
static inline cusolverStatus_t
_cusolverDnGetrf(
    cusolverDnHandle_t handle,
    int m,
    int n,
    double* A,
    int lda,
    double* Workspace,
    int* devIpiv,
    int* devInfo) {
    return cusolverDnDgetrf(handle, m, n, A, lda, Workspace, devIpiv, devInfo);
}
#else
static inline cusolverStatus_t
_cusolverDnGetrf(
    cusolverDnHandle_t handle,
    int m,
    int n,
    float* A,
    int lda,
    float* Workspace,
    int* devIpiv,
    int* devInfo) {
    return cusolverDnSgetrf(handle, m, n, A, lda, Workspace, devIpiv, devInfo);
}
#endif

#if DOUBLE_MODE
static inline cusolverStatus_t
_cusolverDnGetrs(cusolverDnHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const double *A,
                 int lda,
                 const int *devIpiv,
                 double *B,
                 int ldb,
                 int *devInfo ) {
                     return cusolverDnDgetrs(handle, trans, n, nrhs, A, lda, devIpiv, B, ldb, devInfo);
                 }
#else
static inline cusolverStatus_t
_cusolverDnGetrs(cusolverDnHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const float *A,
                 int lda,
                 const int *devIpiv,
                 float *B,
                 int ldb,
                 int *devInfo ) {
                     return cusolverDnSgetrs(handle, trans, n, nrhs, A, lda, devIpiv, B, ldb, devInfo);
                 }
#endif

#if DOUBLE_MODE

#else

#endif