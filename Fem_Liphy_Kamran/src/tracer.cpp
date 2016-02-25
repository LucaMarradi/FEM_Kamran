#include "tracer.h"
#include "fem.h"
#include "stdlib.h"
#include <cassert>
#include "matrix.h" 

Tracer::Tracer( Fem *femObj ):
femPtr( femObj ), memObj( new Memory ),//, nodei( 0 ), ievab( 0 ), 
node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ),
ts( femObj->ts ), fd( femObj->fd ),
msd( 0 )
{}
//------------------------------------------------------
Tracer::~Tracer()
{

//	memObj->Delete2ndMat( posgp, elem->ngaus );
//	memObj->Delete2ndMat( elcod, elem->nnode );
//	memObj->Delete2ndMat( elcod0, elem->nnode );
//	memObj->Delete2ndMat( veloc, elem->nnode );
//	memObj->Delete2ndMat( Cords, 2 );
//	memObj->Delete2ndMat( dimensioless_cords, 2 );
//	memObj->Delete2ndMat( lohi, 2 );
//	memObj->Delete2ndMat( deriv, domain->ndime );
//	memObj->Delete2ndMat( coord, elem->nelem * elem->ngaus );
	
	delete [] weigp;
	delete [] indx;
	delete [] shape;
	delete [] velocity;
	delete mpObj[ 0 ];
	delete mpObj[ 1 ];
	delete r;
	delete theta;
	delete cross_pbc;
	delete [] veloi;
	delete [] velot;
	delete [] dispt;
	delete [] tdisp;		
}
//------------------------------------------------------
void Tracer::Init() //--- interpolate velocities onto tracers
{
	//--- allocate memory
	memObj->doubleMat( posgp, elem->ngaus, domain->ndime );
	memObj->doubleMat( elcod, elem->nnode, domain->ndime ); 
	memObj->doubleMat( elcod0, elem->nnode, domain->ndime ); 
	memObj->doubleMat( veloc, elem->nnode, domain->ndime ); 
	memObj->doubleMat( Cords, 2, domain->ndime );
	memObj->doubleMat( lohi, 2, 5 );
	memObj->doubleMat( dimensioless_cords, 2, domain->ndime ); // allocate memory
	memObj->doubleMat( deriv, domain->ndime, elem->nnode );
	memObj->doubleVec( weigp, elem->ngaus );
	indx = new unsigned int[ domain->ndime ];
	memObj->doubleVec( shape, elem->nnode );
	memObj->doubleVec( velocity, domain->ndime );
	
	// --- some initialization
	// --- box bounds
	lohi[ 0 ][ 0 ] = domain->xlo0;
	lohi[ 0 ][ 1 ] = domain->xhi0;
	lohi[ 0 ][ 2 ] = domain->ylo0;
	lohi[ 0 ][ 3 ] = domain->yhi0;
	lohi[ 0 ][ 4 ] = domain->xy0;
	lohi[ 1 ][ 0 ] = domain->xlo0; //--- need it to get proper disp. 
	lohi[ 1 ][ 1 ] = domain->xhi0;
	lohi[ 1 ][ 2 ] = domain->ylo0;
	lohi[ 1 ][ 3 ] = domain->yhi0;
	lohi[ 1 ][ 4 ] = domain->xy0 + fd->shearRATE * ts->dt * ( domain->yhi0 - domain->ylo0 );
	cordPtr[ 0 ] = node->cordi;
	cordPtr[ 1 ] = node->coord;
	mpObj[ 0 ] = new Wrap( lohi[ 0 ] );
	mpObj[ 1 ] = new Wrap( lohi[ 1 ] );
	r = new double;
	theta = new double;
	cross_pbc = new int;

	elem->GAUSSQ( posgp, weigp ); // gauss points: pos and weight
	kgaus = 0;
			
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			exisp = posgp[ igaus ][ 0 ]; 
			etasp = posgp[ igaus ][ 1 ];
			GetNodalCords( elcod, ielem, 1, 0 ); // xyz with pbc( 1 ) and initial frame( 0 )
			//--- update coord
			Nod2Gas( coord[ kgaus ], elcod, exisp, etasp ); //--- update coord		
			kgaus++;
		} //--- end of gauss loop
	} //--- end of elem. loop
//	printf("Ntracer=%g\t%g\n", coord[ 0 ][ 0 ], dispt[ 0 ] );
	 //--- wrap coord
	Wrap mpObjLoc = lohi[ 0 ];
 	mpObjLoc.CoordWrapper( elem->nelem * elem->ngaus, coord );

}
//------------------------------------------------------
void Tracer::TracerUpdate() //--- interpolate velocities onto tracers
{
	kgaus = 0;
	isvab = 0;
	unsigned int *elemi; //--- elem. id
	elemi = new unsigned int;
	double *exisp, *etasp;
	exisp = new double; 
	etasp = new double; 
	double dldis[ elem->nevab ];

//	double strag[ elem->nstre ];
//	Matrix matPtr;

	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			FindElem( kgaus, elemi ); //--- update elemi
//			printf( "indx = %d\t%d\n", indx[ 0 ], indx[ 1 ] );
			GetNodalCords( elcod, *elemi, 1, 0 ); // xyz with pbc( 1 ) and initial frame( 0 )
			GetExispEtasp( coord[ kgaus ], elcod, exisp, etasp ); //--- update exisp, etasp
			//--- update veloi and dispt
			GetDisp( dldis, *elemi ); // --- output node-based dldis ( from tr->dispt )
//			matPtr.MulMatVec( femPtr->force->bmatx0, elem->nstre, elem->nevab, 
//					   dldis, elem->nevab, strag ); // strain = b * u
			Nod2Gas( velocity, dldis, *exisp, *etasp ); //--- update tracer velocity	
			for( unsigned int idime = 0; idime < domain->ndime; idime++ )
			{
				dispt[ isvab ] += velocity[ idime ]; //--- update dispt
				isvab++;
			}
			kgaus++;
		}
	}
//	printf("exy[ %d ] = %g\n", kgaus, strag[ 2 ] );
	delete elemi;
	delete exisp;
	delete etasp;
}
//-----------------------------
void Tracer::FindElem( unsigned int kgaus, unsigned int *elemPtr )
{
	indx[ 0 ] = ( int )floor( coord[ kgaus ][ 0 ] ) % ( int )( domain->xhi0 - domain->xlo0 );
	indx[ 0 ] += ( int )( domain->xhi0 - domain->xlo0 ) * ( indx[ 0 ] < 0 );
	indx[ 1 ] = ( int )floor( coord[ kgaus ][ 1 ] ) % ( int )( domain->yhi0 - domain->ylo0 );
	indx[ 1 ] += ( int )( domain->xhi0 - domain->xlo0 ) * ( indx[ 1 ] < 0 );
	assert(  0 <= indx[ 0 ] && indx[ 0 ] < ( int )( domain->xhi0 - domain->xlo0 ) );
	assert(  0 <= indx[ 1 ] && indx[ 1 ] < ( int )( domain->yhi0 - domain->ylo0 ) );
	*elemPtr = indx[ 1 ] * ( int )( domain->xhi0 - domain->xlo0 ) + indx[ 0 ]; //--- host element!
	assert( 0 <= *elemPtr && *elemPtr < elem->nelem );
}
//-----------------------------
void Tracer::GetExispEtasp( double *cords, double **elcod, double *exisp, double *etasp )
{
	//--- output exisp & etasp	
	double **cordLoc, **dimensioless_cords;
	cordLoc = new double*[ 1 ];
	cordLoc[ 0 ] = cords;
	dimensioless_cords = new double*[ 1 ];
	dimensioless_cords[ 0 ] = new double[ 2 ];
	
	double lohiLoc[ 5 ] = { elcod[ 0 ][ 0 ], elcod[ 1 ][ 0 ], elcod[ 0 ][ 1 ], elcod[ 2 ][ 1 ], 0.0 }; //--- local lohi (element-based)
	assert( lohiLoc[ 0 ] < lohiLoc[ 1 ] && lohiLoc[ 2 ] < lohiLoc[ 3 ] );
	Wrap mpObjLoc = lohiLoc;
	mpObjLoc.GetDimensionlessCords( 1, cordLoc, dimensioless_cords ); // output dimension...
//	printf( "exisp, etasp = %g\t%g\n", dimensioless_cords[ 0 ][ 0 ], dimensioless_cords[ 0 ][ 1 ] );
	
	*exisp = dimensioless_cords[ 0 ][ 0 ];
	*etasp = dimensioless_cords[ 0 ][ 1 ];
	if( !( 0.0 <= *exisp && *exisp < 1.0 ) ) //--- if not inside
	{
		if( *exisp < 0.0 )
			cords[ 0 ] += ( domain->xhi0 - domain->xlo0 ); //--- shift to right
		else
			cords[ 0 ] -= ( domain->xhi0 - domain->xlo0 );
		mpObjLoc.GetDimensionlessCords( 1, cordLoc, dimensioless_cords ); // output dimension...
		*exisp = dimensioless_cords[ 0 ][ 0 ];
	}
	if( !( 0.0 <= *etasp && *etasp < 1.0 ) )
	{
		if( *etasp < 0.0 )
			cords[ 1 ] += ( domain->yhi0 - domain->ylo0 );
		else
			cords[ 1 ] -= ( domain->yhi0 - domain->ylo0 );
		mpObjLoc.GetDimensionlessCords( 1, cordLoc, dimensioless_cords ); // output dimension...
		*etasp = dimensioless_cords[ 0 ][ 1 ];
	}
	assert( 0.0 <= *exisp && *exisp < 1.0 ); 
	assert( 0.0 <= *etasp && *etasp < 1.0 );
	
	*exisp = 2.0 * ( *exisp - 0.5 ); //--- square: -1 <= ksi < 1
	*etasp = 2.0 * ( *etasp - 0.5 );

	delete [] cordLoc;
	delete [] dimensioless_cords[ 0 ];
	delete [] dimensioless_cords;
}
//-----------------------------
void Tracer::Nod2Gas( double *gausBased, double *nodeBased, double S, double T ) //--- update gausBased
{
	elem->Shape2( shape, deriv, S, T ); //--- output shape, deriv
//	printf( "exisp, etasp = %g\t%g\n", deriv[ 0 ][ 0 ], deriv[ 1 ][ 0 ] );
	for( int idime = 0; idime < domain->ndime; idime++ )
	{
		gausBased[ idime ] = 0.0;
		for( int inode = 0; inode < elem->nnode; inode++ )
			gausBased[ idime ] += shape[ inode ] * nodeBased[ inode * elem->nnode + idime ];
	}
}
//-----------------------------------------------------------------------
//-----------------------------
void Tracer::Nod2Gas( double *gausBased, double **nodeBased, double S, double T ) //--- update gausBased
{
	elem->Shape2( shape, deriv, S, T ); //--- output shape, deriv
//	printf( "exisp, etasp = %g\t%g\n", deriv[ 0 ][ 0 ], deriv[ 1 ][ 0 ] );
	for( int idime = 0; idime < domain->ndime; idime++ )
	{
		gausBased[ idime ] = 0.0;
		for( int inode = 0; inode < elem->nnode; inode++ )
			gausBased[ idime ] += shape[ inode ] * nodeBased[ inode ][ idime ];
	}
}
//-----------------------------------------------------------------------
void Tracer::GetNodalCords( double **elcod, unsigned int ielem, int pbc, int initial )
{
	for( int inode = 0; inode < elem->nnode; inode++ ) // loop over nodes
	{
		nodei = elem->lnods[ ielem ][ inode ]; // node global id
		for( int idime = 0; idime < domain->ndime; idime++ ) // assign coordantes
			elcod[ inode ][ idime ] = cordPtr[ initial ][ nodei ][ idime ];
		if( pbc && inode != 0 ) // --- returns xyz relative to node0 
		{
			Xj = elcod[ 0 ][ 0 ]; 
			Yj = elcod[ 0 ][ 1 ];
			for( int idime = 0; idime < domain->ndime; idime++ )
			{
				Cords[ 0 ][ idime ] = elcod[ inode ][ idime ];
				Cords[ 1 ][ idime ] = elcod[ 0 ][ idime ];
			}
//			printf("hello!\t%d\n", inode );
			Distance( Cords, initial, r, theta, cross_pbc ); // find the true distance
			Xj += ( *r ) * cos( *theta );
			Yj += ( *r ) * sin( *theta );
			elcod[ inode ][ 0 ] = Xj;
			elcod[ inode ][ 1 ] = Yj;
		}
	}
}
//-----------------------------------------------------------------------
void Tracer::GetVeloc( double *output, double **elcod, unsigned int ielem )
{
	for( int inode = 0; inode < elem->nnode; inode++ ) // loop over nodes
	{
		nodei = elem->lnods[ ielem ][ inode ]; // node global id
		for( int idime = 0; idime < domain->ndime; idime++ ) // assign coordantes
			elcod[ inode ][ idime ] = output[ nodei * domain->ndime + idime];
	}
}
//-----------------------------------------------------------------------
void Tracer::Distance(  double **cords,  int initial, double *r, double *theta, int *cross_pbc )
{
	mpObj[ initial ]->GetDimensionlessCords( 2, cords, dimensioless_cords ); // output dimension...

	s0i = dimensioless_cords[ 0 ][ 0 ];
	s1i = dimensioless_cords[ 0 ][ 1 ];
	s0j = dimensioless_cords[ 1 ][ 0 ];
	s1j = dimensioless_cords[ 1 ][ 1 ];
	s0_ij = s0j - s0i; // dimensionless distance
	s1_ij = s1j - s1i;

	// --- translate point j to find the proper distance rij
	*cross_pbc = 0;
        if( s0_ij > 0.5 )
	{ 
		s0j -= 1.0;
		*cross_pbc = 1;
	}
        else if ( s0_ij < - 0.5 )
	{ 
		s0j += 1.0;
		*cross_pbc = 1;
        }
	if( s1_ij > 0.5 )
	{
		s1j -= 1.0;
		*cross_pbc = 1;
        }
	else if( s1_ij < - 0.5 )
	{
		s1j += 1.0;
		*cross_pbc = 1;
	}
	// --- find new coordinates
	dimensioless_cords[ 0 ][ 0 ] = s0i;
	dimensioless_cords[ 0 ][ 1 ] = s1i;
	dimensioless_cords[ 1 ][ 0 ] = s0j;
	dimensioless_cords[ 1 ][ 1 ] = s1j;
	mpObj[ initial ]->GetCords( 2, dimensioless_cords, cords ); // --- update coord

	xi = cords[ 0 ][ 0 ];
	yi = cords[ 0 ][ 1 ];
        xj = cords[ 1 ][ 0 ];
        yj = cords[ 1 ][ 1 ];
        xij = xi - xj;
        yij = yi - yj;
        rsq = xij * xij + yij * yij;
        *r = sqrt( rsq );
	*theta = atan2( yij, xij );
}

//------------------------------------------------------
/*
void Tracer::ApplyHOMO()
{
	kgaus = 0;
	// --- update displacements
	double dGamma = ts->dt * fd->shearRATE;
	// --- initialize dispt
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			dldis[ kgaus * domain->ndime ] = dGamma * ( coord[ kgaus ][ 1 ] - domain->ylo );
			dldis[ kgaus * domain->ndime + 1 ] = 0.0;
			kgaus++;
		}
	}
}
*/
//------------------------------------------------------
void Tracer::UpdateCords()
{
	kgaus = 0;
	isvab = 0;
	// --- update displacements
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{
		for( unsigned int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			for( unsigned int idime = 0; idime < domain->ndime; idime++ )
			{
				coord[ kgaus ][ idime ] += dispt[ isvab ];
				isvab++;
			}
			kgaus++;
		}
		
	}
	//---	map coord
	Wrap mpObjLoc = lohi[ 0 ];
 	mpObjLoc.CoordWrapper( elem->nelem * elem->ngaus, coord );
}
//-----------------------------------------------------------------------
void Tracer::GetDisp( double *dldis,  unsigned int ielem )
{
// --- compute proper dldis
	GetNodalCords( elcod0, ielem, 0, 0 ); // xyz in initial frame( 0 ): no pbc( 0 )
	unsigned int ievab = 0, isvab = 0;	
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		nodei = elem->lnods[ ielem ][ inode ];
		isvab = domain->ndime * nodei;
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			dldis[ ievab ] = tdisp[ isvab ]; // displacements
			elcod0[ inode ][ idofn ] += dldis[ ievab ];  // update original coords
			isvab++;
			ievab++;

		}
	}

	for( int jnode = 1; jnode <  elem->nnode; jnode++ )
	{
		Xj = elcod0[ 0 ][ 0 ];
		Yj = elcod0[ 0 ][ 1 ];
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			Cords[ 0 ][ idime ] = elcod0[ jnode ][ idime ];
			Cords[ 1 ][ idime ] = elcod0[ 0 ][ idime ];
		}
		Distance( Cords, 1, r, theta, cross_pbc ); // find coords in the current frame
		Xj += ( *r ) * cos( *theta );
		Yj += ( *r ) * sin( *theta );
		elcod0[ jnode ][ 0 ] = Xj;
		elcod0[ jnode ][ 1 ] = Yj;
	}	
	// --- proper xyz respecting pbc
	GetNodalCords( elcod, ielem, 1, 0 );
	ievab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			dldis[ ievab ] = elcod0[ inode ][ idime ] - elcod[ inode ][ idime ];
			ievab++;
		}
	}
}
