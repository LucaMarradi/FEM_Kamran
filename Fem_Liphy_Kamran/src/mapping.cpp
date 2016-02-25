//#include "Python.h"
#include "mapping.h"
#include "matrix.h"
#include "math.h"
#include <cassert>
#include <iostream>
using std::cout;

Wrap::Wrap( double *lohi ):
lohiPtr( lohi ), ndime( 2 )
{
	set_h();
	memPtr = new Memory; // allocate memory
	matPtr = new Matrix;
}

Wrap::~Wrap()
{
//	cout << "destructor called!\n";
	delete memPtr;
	delete matPtr;
	delete [] h[ 0 ];
	delete [] h[ 1 ];
	delete [] h;
	delete [] h_inv[ 0 ];
	delete [] h_inv[ 1 ];
	delete [] h_inv;
}
//----------------------------------------------
void Wrap::CoordWrapper(  unsigned int npoin, double **coord )
{
	double **ksiETA;
	double **tmp;
	memPtr->doubleMat( ksiETA, npoin, ndime );
	memPtr->doubleMat( tmp, npoin, ndime );

	GetDimensionlessCords( npoin, coord, ksiETA ); //--- output ksiETA
	// --- assumption: the origin for the xy frame and lattice frmae is the same! 
	double s0, s1;
	for( int ipoin = 0; ipoin < npoin; ipoin++ ){ 
		s1 = ksiETA[ ipoin ][ 1 ];
		s1 = fmod( s1 , 1.0 ) + ( s1 < 0.0 ) * 1.0; //s1 % 1.0
		assert ( 0.0 <= s1 && s1 < 1.0 );

		s0 = ksiETA[ ipoin ][ 0 ];
		s0 = fmod( s0 , 1.0 ) + ( s0 < 0.0 ) * 1.0; //s0 % 1.0
		assert ( 0.0 <= s0 && s0 < 1.0 );

		tmp[ ipoin ][ 0 ] = s0;
		tmp[ ipoin ][ 1 ] = s1;
	}
	GetCords( npoin, tmp, coord ); // --- output coord
	
	for( int ipoin = 0; ipoin < npoin; ipoin++ ){
		delete [] ksiETA[ ipoin ];
		delete [] tmp[ ipoin ];
	}
	delete [] ksiETA;
	delete [] tmp;
}
//----------------------------------------------
void Wrap::GetCords(  unsigned int npoin,  double **dimensioless_cords, double **coord ) //output coord
{
	double coord_mapped[ ndime ];
	for( int ipoin = 0; ipoin < npoin; ipoin++ ){ 
		matPtr->MulMatVec( h, ndime, ndime, dimensioless_cords[ ipoin ], ndime, coord_mapped ); // h * ksi
		coord[ ipoin ][ 0 ] = coord_mapped[ 0 ] + lohiPtr[ 0 ];
		coord[ ipoin ][ 1 ] = coord_mapped[ 1 ] + lohiPtr[ 2 ];
	}
}
//----------------------------------------------
void Wrap::GetDimensionlessCords(  unsigned int npoin,  double **coord, double **dimensioless_cords ) // output dimensioless_cords
{
	double x[ 2 ];
	for( int ipoin = 0; ipoin < npoin; ipoin++ ){ 
		x[ 0 ] = coord[ ipoin ][ 0 ] - lohiPtr[ 0 ];
		x[ 1 ] = coord[ ipoin ][ 1 ] - lohiPtr[ 2 ];
		matPtr->MulMatVec( h_inv, ndime, ndime, x, ndime, dimensioless_cords[ ipoin ] ); // h_inv * x
	}	
}
//----------------------------------------------
void Wrap::set_h()
{
	double xlo = lohiPtr[ 0 ];
	double xhi = lohiPtr[ 1 ];
	double ylo = lohiPtr[ 2 ];
	double yhi = lohiPtr[ 3 ];
	double xy  = lohiPtr[ 4 ];
	memPtr->doubleMat( h, ndime, ndime );
	h[ 0 ][ 0 ] = xhi - xlo;
	h[ 0 ][ 1 ] = xy;
	h[ 1 ][ 0 ] = 0.0;
	h[ 1 ][ 1 ] = yhi - ylo;
	double det = ( yhi - ylo ) * ( xhi - xlo );
	memPtr->doubleMat( h_inv, ndime, ndime );
	h_inv[ 0 ][ 0 ] = h[ 1 ][ 1 ];
	h_inv[ 0 ][ 1 ] = -h[ 0 ][ 1 ];
	h_inv[ 1 ][ 0 ] = -h[ 1 ][ 0 ];
	h_inv[ 1 ][ 1 ] =  h[ 0 ][ 0 ];
	for( int i = 0; i < ndime; i++ ){
		for( int j = 0; j < ndime; j++ ){
			h_inv[ i ][ j ] /=  det;
		}
	}
}

/*
int main()
{
	double lohi[ 5 ] = { 0.0, 1.0, 0.0, 1.0, 0.0 };
	Wrap *obj = new Wrap( lohi );
	Memory *memPtr;
	double **coord;
	memPtr->doubleMat( coord, 1, 2 ); 
	coord[ 0 ][ 0 ] = -1.5;
	coord[ 0 ][ 1 ] = 0.5;
	obj->CoordWrapper( 1, coord );
	cout << coord[ 0 ][ 0 ] << '\t' << coord[ 0 ][ 1 ] << '\n';
}
*/
