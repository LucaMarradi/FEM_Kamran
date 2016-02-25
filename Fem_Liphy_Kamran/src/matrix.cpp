#include "matrix.h"
#include <cassert>
//#include <iostream>

//using std::cout;

/**
 * Matrix operations should be performed using the
 * PBLAS library (parallel version)
 * @link http://www.netlib.org/scalapack/pblas_qref.html
 */

Matrix::Matrix(){}
Matrix::~Matrix(){}

/**
 * [Matrix::Transpose description]
 * @param[in] A matrix to Transpose
 * @param[in] m number of columns
 * @param[in] n number of rows
 * @param[out] B matrix Transposed
 */
void Matrix::Transpose(  double **A,  unsigned int m,  unsigned int n, double **B )
{
	for( int i = 0; i < n; i++ ){
		for( int j = 0; j < m; j++ ){
		B[ i ][ j ] = A[ j ][ i ];
		}
	}
}

/**
 * [Matrix::MulMatVec description]
 * @param[in] A input matrix
 * @param[in] m number of rows
 * @param[in] n number of columns
 * @param u [description]
 * @param o [description]
 * @param[out] v output vector
 */
void Matrix::MulMatVec(  double **A,  unsigned int m,  unsigned int n,  double *u,  unsigned int o, double *v )
{
	assert( n == o );
	double sum = 0.0;
	for( int i = 0; i < m; i++ ){
		sum = 0.0;
		for( int j = 0; j < n; j++ ){
			sum += A[ i ][ j ] * u[ j ];
		}
		v[ i ] = sum;
	}
}


/**
 * [Matrix::MulMatMat description]
 * @author Kamran Karimi
 * @param[in] A input matrix
 * @param[in] m number of C-rows
 * @param[in] n number of B-rows and A-columns
 * @param[in] B input matrix
 * @param o [description]
 * @param p [description]
 * @param[out] C output matrix
 *
 * @brief The function multiply couple of Matrix. Noticing
 *        that the the number of columns of A matrix must
 *        be equal to the B one.
 */
void Matrix::MulMatMat(  double **A,  unsigned int m,  unsigned int n,
			 double **B,  unsigned int o,  unsigned int p, double **C )
{
	assert( n == o );
	double sum = 0.0;
	for( int i = 0; i < m; i++ ){
		for( int k = 0; k < p; k++ ){
			sum = 0.0;
			for( int j = 0; j < n; j++ )
				sum += A[ i ][ j ] * B[ j ][ k ];
			C[ i ][ k ] = sum;
		}
	}
}

/*
int main()
{
	Matrix *obj;

	double **A = new double*[ 2 ];
	A[ 0 ] = new double[ 2 ];
	A[ 1 ] = new double[ 2 ];
	A[ 0 ][ 0 ] = 0.0;
	A[ 0 ][ 1 ] = 1.0;
	A[ 1 ][ 0 ] = 2.0;
	A[ 1 ][ 1 ] = 3.0;

	double *u = new double[ 2 ];
	u[ 0 ] = 1.0;
	u[ 1 ] = 2.0;

	double *v = new double[ 2 ];

	obj->MulMatVec( A, 2, 2, u, 2, v );
	cout << v[ 0 ] << '\n' << v[ 1 ] << '\n';

	delete [] A;
	delete [] u;
	delete [] v;
}
*/
