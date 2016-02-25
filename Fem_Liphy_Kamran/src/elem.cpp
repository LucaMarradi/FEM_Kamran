#include "elem.h"
#include "fem.h"
#include "gauss_coord.h"
#include "stdlib.h"
#include "matrix.h"
#include <cassert>
#include <iostream>

using std::cout;

Elem::Elem( Fem *femObj ):
femPtr( femObj ), 
memPtr( new Memory ), randOBJ( new RandNumbers( lamda, rand() ) ),
domain( femObj->domain )
{}

Elem::~Elem()
{
	delete [] ltype;
	memPtr->Delete2ndMat( lnods, nelem );
	delete [] kmodu;
	delete [] gmodu;
	delete [] uniax;
	delete [] strny;
	memPtr->Delete2ndMat( stran, nelem * ngaus );
	memPtr->Delete2ndMat( strni, nelem * ngaus );
	memPtr->Delete2ndMat( strst, nelem * ngaus );
	memPtr->Delete2ndMat( strsi, nelem * ngaus );
	delete [] epstn;
	memPtr->Delete2ndMat( stranRate, nelem * ngaus );
	memPtr->Delete2ndMat( plasticStrain, nelem * ngaus );
	memPtr->Delete2ndMat( plasticStrain0, nelem * ngaus );
	delete [] vol;
	delete [] actField;

	memPtr->Delete2ndMat( deriv, domain->ndime );
	memPtr->Delete2ndMat( XJACM, domain->ndime );
	memPtr->Delete2ndMat( XJACI, domain->ndime );
	delete [] shape;
	delete memPtr;
	delete randOBJ;

	memPtr->Delete2ndMat( cordg, nelem * ngaus );
}

void Elem::Init()
{
//	printf( "ndime=%d\n", femPtr->domain->ndime );
//	printf( "nnode=%d\n", nnode );
	memPtr->doubleVec( shape, nnode );
	memPtr->doubleMat( deriv, domain->ndime, nnode );
	memPtr->doubleMat( XJACM, domain->ndime, domain->ndime );
	memPtr->doubleMat( XJACI, domain->ndime, domain->ndime );
	memPtr->doubleMat( cordg, nelem * ngaus, domain->ndime ); // allocate cordg	
	gcObj = new GaussCoord( femPtr );
	gcObj->Init();
	gcObj->GetCoord( cordg );
	minE  = 0.25 * meanGmodu * 0.035 * 0.035; //--- hard-coded!
}
//------------------------------------------------------------------
void Elem::ksiTOxy( double *output,  double *ksi,  double **nodalXY )
{
//	double *shape, **deriv;
	Shape2( shape, deriv, ksi[ 0 ], ksi[ 1 ] ); // output shape & deriv functions 
	for( int idime = 0; idime < domain->ndime; idime++ )
		for( int inode = 0; inode < nnode; inode++ )
			output[ idime ] += shape[ inode ] * nodalXY[ inode ][ idime ];
//	delete [] shape;
//	for( int idime = 0; idime < ndime; idime++ ){
//		delete [] deriv[ idime ];
//	}
//	delete [] deriv;

}
//------------------------------------------------------------------
void Elem::set(  unsigned int ielem )
{
	// set gmodu and kmodu
	gmodu[ ielem ] = meanGmodu;//femPtr->elem->meanGmodu;
	kmodu[ ielem ] = meanKmodu; //femPtr->elem->meanKmodu;
//	double minE  = 0.25 * femPtr->elem->meanGmodu * 0.07 * 0.07; //kam
	// set yield strain	
//	randOBJ->Expon( &( uniax[ ielem ] ), lamda );
	randOBJ->Expon( &( uniax[ ielem ] ) );
	uniax[ ielem ] += minE;
	assert( uniax[ ielem ] >= minE );
	uniax[ ielem ] = 2.0 * sqrt( gmodu[ ielem ] * uniax[ ielem ] ); // energy to stress
	strny[ ielem ] = meanRecoveryStrain; //femPtr->elem->meanRecoveryStrain;
}
//------------------------------------------------------------------
void Elem::Shape2( double *shape,  /*output*/
		   double **deriv, /*output*/
		    double S, 
		    double T )
{
// --- ( S, T ) material coordinates ranging between -1 and 1
// --- SHAPE FUNCTIONS AND DERIVATIVES FOR LINEAR TRIANGLE

	if( nnode == 3 )
	{
		p = 1.0 - S - T;
		shape[ 0 ] = p; 
		shape[ 1 ] = S;
		shape[ 2 ] = T;

		deriv[ 0 ][ 0 ] = -1.0;
		deriv[ 0 ][ 1 ] = 1.0;
		deriv[ 0 ][ 2 ] = 0.0;
		deriv[ 1 ][ 0 ] = -1.0;
		deriv[ 1 ][ 1 ] = 0.0;
		deriv[ 1 ][ 2 ] = 1.0;
	}
	else if( nnode == 4 ){
		S2 = S * 2.0;
		T2 = T * 2.0;
		SS = S * S;
		TT = T * T;
		ST = S * T;
		SST = S * S * T;
		STT = S * T * T;
		ST2 = S * T * 2.0;
//--- SHAPE FUNCTIONS FOR 4 NODED ELEMENT
		shape[ 0 ] = ( 1.0 - S - T + ST ) * 0.25;
		shape[ 1 ] = ( 1.0 + S - T - ST ) * 0.25;
		shape[ 2 ] = ( 1.0 + S + T + ST ) * 0.25;
		shape[ 3 ] = ( 1.0 - S + T -ST ) * 0.25;
//--- AND DERIVATIVES
		deriv[ 0 ][ 0 ] = ( -1.0 + T ) * 0.25;
		deriv[ 0 ][ 1 ] = -deriv[ 0 ][ 0 ];
		deriv[ 0 ][ 2 ] = ( 1.0 + T ) * 0.25;
		deriv[ 0 ][ 3 ] = -deriv[ 0 ][ 2 ];
		deriv[ 1 ][ 0 ] = ( -1.0 + S ) * 0.25;
		deriv[ 1 ][ 1 ] = ( -1.0 - S ) * 0.25;
		deriv[ 1 ][ 2 ] = -deriv[ 1 ][ 1 ];
		deriv[ 1 ][ 3 ] = -deriv[ 1 ][ 0 ];
	}
}
//------------------------------------------------------------------------
void Elem::Blarge( double **bmatx, /*output*/
		    double **cartd )
{
	ngash = -1;
	for( int inode = 0; inode < nnode; inode++ )
	{
		mgash = ngash + 1;
		ngash = mgash + 1;
		bmatx[ 0 ][ mgash ] = cartd[ 0 ][ inode ];
		bmatx[ 0 ][ ngash ] = 0.0;
		bmatx[ 1 ][ mgash ] = 0.0;
		bmatx[ 1 ][ ngash ] = cartd[ 1 ][ inode ];
		bmatx[ 2 ][ mgash ] = cartd[ 1 ][ inode ]; //* 0.5
		bmatx[ 2 ][ ngash ] = cartd[ 0 ][ inode ]; //* 0.5
	}	
}
//------------------------------------------------------------------------
void Elem::Jacob2( double **cartd, /*output*/
		   double *djacb,  /*output*/
		    double S, 
		    double T, 
		    double **elcod, 
		    unsigned int ielem )
{
//	double *shape, **deriv;

//	int ndime = femPtr->domain->ndime;
//	memPtr->doubleVec( shape, nnode );
//	memPtr->doubleMat( deriv, ndime, nnode );
//	memPtr->doubleMat( XJACM, ndime, ndime );
//	memPtr->doubleMat( XJACI, ndime, ndime );
		
	Shape2( shape, deriv, S, T ); // output shape, deriv
	matPtr->MulMatMat( deriv, domain->ndime, nnode, elcod, nnode, domain->ndime, XJACM ); //XJACM = deriv * elcod
//	for( int i = 0; i < 4; i++ ){
//		for( int j = 0; j < 2; j++ )
//			cout << elcod[ i ][ j ] << '\t';
//		cout << '\n';
//	}
	*djacb = det( XJACM ); 
//	cout << *djacb << '\n';
//	exit( EXIT_FAILURE );
	if( *djacb <= 0.0 )
	{
		printf( "STOP IN JACOB2: element %d\t%f\n", ielem, *djacb );
		exit( EXIT_FAILURE );
	}
	inv( XJACI, XJACM ); //inverse
	matPtr->MulMatMat( XJACI, domain->ndime, domain->ndime, deriv, domain->ndime, nnode, cartd ); //output cartd = XJACI * deriv

//	delete [] shape;
//	for( int idime = 0; idime < ndime; idime++ ){
//		delete [] deriv[ idime ];
//		delete [] XJACM[ idime ];
//		delete [] XJACI[ idime ];
//	}
//	delete [] deriv;
//	delete [] XJACM;
//	delete [] XJACI;
}

//-----------------------------------------------------------------------------
inline void Elem::inv( double **XJACI,  double **XJACM )
{
	determinant = XJACM[ 0 ][ 0 ] * XJACM[ 1 ][ 1 ] - XJACM[ 0 ][ 1 ] * XJACM[ 1 ][ 0 ];
	assert( determinant != 0.0 );
	XJACI[ 0 ][ 0 ] = XJACM[ 1 ][ 1 ] / determinant;
	XJACI[ 1 ][ 1 ] = XJACM[ 0 ][ 0 ] / determinant;
	XJACI[ 0 ][ 1 ] = -XJACM[ 0 ][ 1 ] / determinant;
	XJACI[ 1 ][ 0 ] = -XJACM[ 1 ][ 0 ] / determinant;
}
//-----------------------------------------------------------------------------
inline double Elem::det(  double ** XJACM )
{
	determinant = XJACM[ 0 ][ 0 ] * XJACM[ 1 ][ 1 ] - XJACM[ 0 ][ 1 ] * XJACM[ 1 ][ 0 ];
	return determinant;
}
//-----------------------------------------------------------------------------
void Elem::GAUSSQ( double **posgp, double *weigp )
{
//	nnode = femPtr->elem->nnode;
//	ngaus = femPtr->elem->ngaus;
//	int ndime = femPtr->domain->ndime;
	if( nnode == 3 )
	{
		posgp[ 0 ][ 0 ] = 1.0 / 3.0;
		posgp[ 0 ][ 1 ] = 1.0 / 3.0;
		weigp[ 0 ] = 0.5;
	}
	else if( nnode == 4 && ngaus == 4 )
	{
		double G = 0.577350269189626;
		posgp[ 0 ][ 0 ] = -1.0;
		posgp[ 0 ][ 1 ] = -1.0;
		posgp[ 1 ][ 0 ] = -1.0;
		posgp[ 1 ][ 1 ] = 1.0;
		posgp[ 2 ][ 0 ] = 1.0;
		posgp[ 2 ][ 1 ] = -1.0;
		posgp[ 3 ][ 0 ] = 1.0;
		posgp[ 3 ][ 1 ] = 1.0;
		for( int igaus = 0; igaus < ngaus; igaus++ )
		{
			weigp[ igaus ] = 1.0;
			for( int idime = 0; idime < domain->ndime; idime++ )
				posgp[ igaus ][ idime ] *= G;
		}
	}
	else if( nnode == 4 && ngaus == 1 )
	{
	        posgp[ 0 ][ 0 ] = 0.0;
	        posgp[ 0 ][ 1 ] = 0.0;
	        weigp[ 0 ] = 4.0;
	}

	
}
