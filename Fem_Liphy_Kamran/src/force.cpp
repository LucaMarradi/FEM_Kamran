#include <iostream>
#include "force.h"
#include "fem.h"
#include "domain.h"
#include "math.h"
#include "mapping.h"
#include "stdlib.h" // exit, rand,srand
#include "matrix.h"
#include "time.h" // time
#include "cassert" // assert

using std::cout;

Force::Force( Fem *femObj ):
femPtr( femObj ), memObj( new Memory ), nodei( 0 ), ievab( 0 ), 
isvab( 0 ), kgaus( 0 ), itype( 0 ),
node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ), 
fv( femObj->fv ), ts( femObj->ts ), fd( femObj->fd )
{}
//--------------------------------------------------------------------
Force::~Force()
{
	memObj->Delete2ndMat( posgp, elem->ngaus );
	memObj->Delete2ndMat( elcod, elem->nnode );
	memObj->Delete2ndMat( elcod0, elem->nnode );
	memObj->Delete2ndMat( Dmatx, elem->nstre );
	memObj->Delete2ndMat( Vmatx, elem->nstre );
	memObj->Delete2ndMat( cartd, domain->ndime );
	memObj->Delete3rdMat( bmatx, elem->nelem * elem->ngaus, elem->nstre  );
	memObj->Delete3rdMat( bt, elem->nelem * elem->ngaus, elem->nevab );
	memObj->Delete2ndMat( Cords, 2 );
	memObj->Delete2ndMat( dev, 2 );
	memObj->Delete2ndMat( dimensioless_cords, 2 );
	memObj->Delete2ndMat( lohi, 2 );
	delete [] weigp;
	delete [] dldis;
	delete [] veloc;
	delete [] strag;
	delete [] dsigm;
	delete [] sigmV;
	delete [] sgtri;
	delete [] bt_sigma;
	delete [] eload;
	delete [] deviatoric_stress;
	delete [] dev_trial;
	delete [] dev_strsi;
	delete [] dev_strag;
	delete djacb;
	delete r;
	delete theta;
	delete cross_pbc;
	delete mpObj[ 0 ];
	delete mpObj[ 1 ];
	delete memObj;
}
//--------------------------------------------------------------------
void Force::Init()
{
	//--- allocate memory
	memObj->doubleMat( posgp, elem->ngaus, domain->ndime );
	memObj->doubleMat( elcod, elem->nnode, domain->ndime ); 
	memObj->doubleMat( elcod0, elem->nnode, domain->ndime ); 
	memObj->doubleMat( Dmatx, elem->nstre, elem->nstre ); 
	memObj->doubleMat( Vmatx, elem->nstre, elem->nstre ); 
	memObj->doubleMat( cartd, domain->ndime, elem->nnode ); 
	bmatx = new double**[ elem->nelem * elem->ngaus ];
	bt    = new double**[ elem->nelem * elem->ngaus ];
	for( int igaus = 0; igaus < elem->nelem * elem->ngaus; igaus++ )
	{
		memObj->doubleMat( bmatx[ igaus ], elem->nstre, elem->nevab ); 
		memObj->doubleMat( bt[ igaus ], elem->nevab, elem->nstre ); 
	}	
	memObj->doubleMat( Cords, 2, domain->ndime );
	memObj->doubleMat( dev, 2, 2 );
	memObj->doubleMat( dimensioless_cords, 2, domain->ndime ); // allocate memory
	memObj->doubleVec( weigp, elem->ngaus );
	memObj->doubleVec( dldis, elem->nevab );
	memObj->doubleVec( veloc, elem->nevab );
	memObj->doubleVec( strag, elem->nstre );
	memObj->doubleVec( dsigm, elem->nstre );
	memObj->doubleVec( sigmV, elem->nstre );
	memObj->doubleVec( sgtri, elem->nstre );
	memObj->doubleVec( bt_sigma, elem->nevab );
	memObj->doubleVec( eload, elem->nevab );
	memObj->doubleMat( lohi, 2, 5 );
	memObj->doubleVec( deviatoric_stress, elem->nstre );
	memObj->doubleVec( dev_trial, elem->nstre );
	memObj->doubleVec( dev_strsi, elem->nstre );
	memObj->doubleVec( dev_strag, elem->nstre );
	djacb = new double;
	r = new double;
	theta = new double;
	cross_pbc = new int;

	// --- box bounds
	lohi[ 0 ][ 0 ] = domain->xlo0;
	lohi[ 0 ][ 1 ] = domain->xhi0;
	lohi[ 0 ][ 2 ] = domain->ylo0;
	lohi[ 0 ][ 3 ] = domain->yhi0;
	lohi[ 0 ][ 4 ] = domain->xy0;

	lohi[ 1 ][ 0 ] = domain->xlo0;
	lohi[ 1 ][ 1 ] = domain->xhi0;
	lohi[ 1 ][ 2 ] = domain->ylo0;
	lohi[ 1 ][ 3 ] = domain->yhi0;
	lohi[ 1 ][ 4 ] = domain->xy0 + fd->shearRATE * ts->dt * ( domain->yhi0 - domain->ylo0 );
	// --- some initialization
	cordPtr[ 0 ] = node->cordi;
	cordPtr[ 1 ] = node->coord;
	mpObj[ 0 ] = new Wrap( lohi[ 0 ] );
	mpObj[ 1 ] = new Wrap( lohi[ 1 ] );
		
	// --- set up element obj
	elem->Init();
	elem->GAUSSQ( posgp, weigp ); // gauss points: pos and weight
	visc_coeff = fv->cdrag0;
	second_visc_coeff = fv->cdrag1;

	// --- element loop starts: update vol[ kgaus ], bmat[ kgaus ] , bt[ kgaus ]
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{ 
		GetNodalCords( elcod, ielem, 1, 0 ); // xyz with pbc( 1 ) and initial frame( 0 ) 
/*		if( ielem == 33 )
		{
			for( int i = 0; i < elem->nnode; i++ )
			{
				nodei = elem->lnods[ ielem ][ i ];
				printf( "x[ %d ]\t", i );
				for( int j = 0; j < domain->ndime; j++ )
					printf( "%g\t", elcod[ i ][ j ] );
				for( int j = 0; j < domain->ndime; j++ )
					printf( "%g\t", node->coord[ nodei ][ j ] );
				printf( "\n" );
			}
		}
*/		// loop over quad points starts
		for( int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			exisp = posgp[ igaus ][ 0 ]; 
			etasp = posgp[ igaus ][ 1 ];
			elem->Jacob2( cartd, djacb, exisp, etasp, elcod, ielem ); // output cartd, djacb
			dvolu = ( *djacb ) * weigp[ igaus ];
			elem->vol[ kgaus ] = dvolu;
			elem->Blarge( bmatx[ kgaus ], cartd ); // output bmatx
			matPtr->Transpose( bmatx[ kgaus ], elem->nstre, elem->nevab, bt[ kgaus ] ); // output bt
			kgaus++;
		} //--- end of quad. loop

	} // --- end of element loop
}
//--------------------------------------------------------------------
void Force::ComputeForce( int counter )
{
	int deform = 0; //--- q-st
	if( counter == 1 )
        	deform = 1;

//	double tmp[ nstre ];
//	double dmass = 0.0, imass = 0.0, emass[ nevab ];
//	double rho = 1.0 * node->npoin;
//	double shape[ nnode ], **deriv;
//	memObj->doubleMat( deriv, ndime, nnode );

	// --- initialize fintl
	for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ ) 
		node->fintl[ isvab ] = 0.0;
	// --- quad. counter
	kgaus = 0;
//	for( unsigned int ipoin = 0; ipoin < node->npoin; ipoin++ ) 
//		node->mass[ ipoin ] = 0.0;

//	double vol = 0.0;
//	double sum = 0.0;		
//	cout << node->veloi[ 0 ] << '\n';	
//			if( kgaus == 0 )
	// --- element loop starts
	for( unsigned int ielem = 0; ielem < elem->nelem; ielem++ )
	{ 
		// --- initialize eload
		for( int iEvab = 0; iEvab < elem->nevab; iEvab++ )
		{
			eload[ iEvab ] = 0.0;
			dldis[ iEvab ] = 0.0;
			veloc[ iEvab ] = 0.0;
//			emass[ iEvab ] = 0.0;
		}
		GetDisp( dldis, ielem, deform ); // --- output dldis
		GetVeloc( veloc, ielem ); // --- output veloc
		
/*		if( ielem == 33 )
		{
			for( int i = 0; i < elem->nevab; i++ )
				printf( "x[ %d ] = %g\n", i, dldis[ i ] );
		}
*/		
		// --- computes elasticity matrix
		itype = elem->ltype[ ielem ];
		strny = elem->strny[ itype ];
		kmodu = elem->kmodu[ itype ];
		gmodu = elem->gmodu[ itype ];
		uniax = elem->uniax[ itype ];
		Modulus( Dmatx, kmodu, gmodu, elem->nstre );
		Modulus( Vmatx, second_visc_coeff, visc_coeff, elem->nstre ); //--- visc. matrix 
		// loop over quad points starts
		for( int igaus = 0; igaus < elem->ngaus; igaus++ )
		{
			dvolu = elem->vol[ kgaus ]; 
//			vol += dvolu;
			// --- stress and plastic strain update
			matPtr->MulMatVec( bmatx[ kgaus ], elem->nstre, elem->nevab, 
					   dldis, elem->nevab, strag ); // strain = b * u
			matPtr->MulMatVec( bmatx[ kgaus ], elem->nstre, elem->nevab, 
					   veloc, elem->nevab, elem->stranRate[ kgaus ] ); // strain_rate = b * v
			matPtr->MulMatVec( Dmatx, elem->nstre, elem->nstre, 
					   strag, elem->nstre, dsigm ); // dsigma = D * strain
			for( int istre = 0; istre < elem->nstre; istre++ ) 
				sigmV[ istre ] = 0.0;
			matPtr->MulMatVec( Vmatx, elem->nstre, elem->nstre, 
					   elem->stranRate[ kgaus ], elem->nstre, sigmV ); //--- viscous sigma
			for( int istre = 0; istre < elem->nstre; istre++ )
			{
				elem->strst[ kgaus ][ istre ] = elem->strsi[ kgaus ][ istre ]; // initialize stress
				sgtri[ istre ] = dsigm[ istre ] + elem->strsi[ kgaus ][ istre ]; // trial stress
				elem->stran[ kgaus ][ istre ] = elem->strni[ kgaus ][ istre ]; // total strain
				elem->plasticStrain[ kgaus ][ istre ] = elem->plasticStrain0[ kgaus ][ istre ]; // initialize plastic strain
			}
/*			if( ielem == 33 && igaus == 0 )
			{
				for( int i = 0; i < elem->nstre; i++ )
					printf( "e[ %d ] = %g\n", i, strag[ i ] );
			} */
			ReturnStress( elem->strst[ kgaus ], elem->stran[ kgaus ], 
				      elem->plasticStrain[ kgaus ], &( elem->actField[ kgaus ] ), ielem ); // update strst, stran, plasticStrain 

		//--- CALCULATE THE EQUIVALENT NODAL FORCES
			for( int istre = 0; istre < elem->nstre; istre++ ) 
				sigmV[ istre ] = sigmV[ istre ] * ( 1 - elem->actField[ kgaus ] ) + elem->strst[ kgaus ][ istre ]; //--- only elastic elemnts (kelvin type)
//			matPtr->MulMatVec( bt[ kgaus ], elem->nevab, elem->nstre, 
//					   elem->strst[ kgaus ], elem->nstre, bt_sigma ); //bt * sigma
			matPtr->MulMatVec( bt[ kgaus ], elem->nevab, elem->nstre, 
					   sigmV, elem->nstre, bt_sigma ); // ---bt * sigma
			for( int ievab = 0; ievab < elem->nevab; ievab++ )
				eload[ ievab ] += bt_sigma[ ievab ] * dvolu;
			kgaus++;
		} // gauss points lopp ends
		ievab = 0;
		for( int inode = 0; inode < elem->nnode; inode++ )
		{
			nodei = elem->lnods[ ielem ][ inode ];
			isvab = domain->ndime * nodei;
//			node->mass[ nodei ] += emass[ ievab ]; // update mass 
			for( int idofn = 0; idofn < domain->ndime; idofn++ )
			{
//				if( isvab == 0 )
//					printf( "f[ %d ] = %15.14e\n", ievab, eload[ ievab ] );
				node->fintl[ isvab ] += eload[ ievab ];
				isvab++;
				ievab++;
			}
		}
	} // element loop ends
//	exit( EXIT_FAILURE );
}
//-----------------------------------------------------------------------
inline double Force::Smean(  double *stress )
{
	return 0.5 * ( stress[ 0 ] + stress[ 1 ] );
}
//-----------------------------------------------------------------------
inline void Force::DeviaStrain( double *dev /*output*/,  double *T )
{
	dev[ 0 ] = T[ 0 ] - T[ 1 ];
	dev[ 1 ] = - dev[ 0 ];
	dev[ 2 ] = T[ 2 ]; // [ ( exx - eyy ), ( eyy - exx ), 2.0 * exy ]
}
//-----------------------------------------------------------------------
inline double Force::EffectiveStrain(  double *strain )
{
	exx_yy = strain[ 0 ];
	exy = 0.5 * strain[ 2 ];
	return sqrt( exx_yy * exx_yy + 4.0 * exy * exy ); // gamma^2 = 2 e_dev : e_dev
}
//-----------------------------------------------------------------------
inline void Force::Devia( double *dev /*output*/,  double *T )
{
	smean = Smean( T );
	dev[ 0 ] = T[ 0 ] - smean;
	dev[ 1 ] = T[ 1 ] - smean;
	dev[ 2 ] = T[ 2 ];
}
//-----------------------------------------------------------------------
double Force::Invart(  double *u )
{
// --- invariants of tensor T

	smean = Smean( u );
	dev[ 0 ][ 0 ] = u[ 0 ] - smean;
	dev[ 0 ][ 1 ] = u[ 2 ];
	dev[ 1 ][ 0 ] = u[ 2 ];
	dev[ 1 ][ 1 ] = u[ 1 ] - smean;
	varj2 = dev[ 0 ][ 1 ] * dev[ 1 ][ 0 ] + \
	0.5 * ( dev[ 0 ][ 0 ] * dev[ 0 ][ 0 ] + dev[ 1 ][ 1 ] * dev[ 1 ][ 1 ] );
	return sqrt( varj2 );
	
}
//-----------------------------------------------------------------------
void Force::ReturnStress( double *stress, double *stran, double *plasticStrain, int *actField /*output*/,  unsigned int ielem )
{

	effectiveStress = Invart( sgtri );	
	effectivePlasStrain = EffectiveStrain( plasticStrain );	
	yfunc = effectiveStress - uniax;

//	if( kgaus == 0 )
//		printf( "yf=%g\n", yfunc );
//	if( yfunc >= 0.0 || effectivePlasStrain > 0.0 ) // --- elastic to plastic or plastic-plastic!
	if( yfunc >= 0.0 || *actField || effectivePlasStrain > 0.0 ) // --- elastic to plastic or plastic-plastic!
	{
		*actField = 1;
		Devia( dev_trial, sgtri );
		Devia( dev_strsi, stress );
		DeviaStrain( dev_strag, strag );
		for( int istre = 0; istre < elem->nstre; istre++ )
		{
			deviatoric_stress[ istre ] = dev_trial[ istre ] - dev_strsi[ istre ] * ts->dt * gmodu / visc_coeff;
			plasticStrain[ istre ] += fabs( ts->dt * dev_strsi[ istre ] / visc_coeff );
			stran[ istre ] += fabs( dev_strag[ istre ] ); //dt * dev_strsi[ istre ] / visc_coeff;
		}
		pstre = Smean( sgtri );
		stress[ 0 ] = deviatoric_stress[ 0 ] + pstre;
		stress[ 1 ] = deviatoric_stress[ 1 ] + pstre;
		stress[ 2 ] = deviatoric_stress[ 2 ];
	}
	else 
	{
		for( int istre = 0; istre < elem->nstre; istre++ )
			stress[ istre ] = sgtri[ istre ]; //--- elastic-elastic
	}
	
	if ( ( *actField ) && EffectiveStrain( stran ) >= strny )
//	if ( ( *actField ) && EffectiveStrain( plasticStrain ) - 0.35 * uniax / gmodu >= 0.0  )
//	if ( ( *actField ) && EffectiveStrain( plasticStrain ) - uniax / gmodu >= 0.0  )
	{
		for( int istre = 0; istre < elem->nstre; istre++ )
		{
			plasticStrain[ istre ] = 0.0;
			stran[ istre ] = 0.0;
			*actField = 0;
		}
		elem->set( ielem ); // new material parameters
	}
			
}

//-----------------------------------------------------------------------
inline void Force::Modulus( double **Dmatx, 
		      double kmod, 
		      double gmod, 
		      int    nstre )
{	
	Dmatx[ 0 ][ 0 ] = kmod + gmod;
	Dmatx[ 0 ][ 1 ] = kmod - gmod;
	Dmatx[ 0 ][ 2 ] = 0.0;

	Dmatx[ 1 ][ 0 ] = kmod - gmod;
	Dmatx[ 1 ][ 1 ] = kmod + gmod;
	Dmatx[ 1 ][ 2 ] = 0.0;
		
	Dmatx[ 2 ][ 0 ] = 0.0;
	Dmatx[ 2 ][ 1 ] = 0.0;
	Dmatx[ 2 ][ 2 ] = 1.0 * gmod;
}
//-----------------------------------------------------------------------
void Force::GetVeloc( double *veloc,  unsigned int ielem )
{
//	unsigned int nodei, isvab;
	ievab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		nodei = elem->lnods[ ielem ][ inode ];
		isvab = domain->ndime * nodei;
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			veloc[ ievab ] = node->velot[ isvab ]; // displacements
			isvab++;
			ievab++;
		}
	}
}
//-----------------------------------------------------------------------
void Force::GetDisp( double *dldis,  unsigned int ielem, int deform )
{
// --- find proper displacements for stress update
//	double **lcord0; // xyz in the current frame
//	double **lcord;  // xyz in the initial frame
//	int ievab = 0;   // counter for element dof's
//	unsigned int isvab; // counter for global dof's
//	unsigned int nodei; // node id
//	double **cords; // xyz for pairs of nodes 
//	double *r = new double; // bond length
//	double *theta = new double; // bond angle
//	int *cross_pbc = new int; // cross PBC?


//	memObj->doubleMat( lcord0, elem->nnode, domain->ndime ); // allocation
//	memObj->doubleMat( lcord, elem->nnode, domain->ndime ); // allocation
//	memObj->doubleMat( cords, 2, domain->ndime );

// --- compute proper dldis
	GetNodalCords( elcod0, ielem, 0, 0 ); // xyz in initial frame( 0 ): no pbc( 0 )
/*	if( ielem == 33 )
	{
		printf( "getDisp\n" );
		for( int i = 0; i < elem->nnode; i++ )
		{
			printf( "x[ %d ]\t", i );
			for( int j = 0; j < domain->ndime; j++ )
				printf( "%g\t", elcod0[ i ][ j ] );
			printf( "\n" );
		}
		printf( "\n" );
	} */
//	if( ielem == 28 ){
//		cout << "initial\n";
//		cout << lcord0[ 2 ][ 0 ] << '\t' << lcord0[ 2 ][ 1 ] << '\n';
//	}
	ievab = 0;	
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		nodei = elem->lnods[ ielem ][ inode ];
		isvab = domain->ndime * nodei;
		for( int idofn = 0; idofn < domain->ndime; idofn++ )
		{
			dldis[ ievab ] = node->dldis[ isvab ]; // displacements
			elcod0[ inode ][ idofn ] += dldis[ ievab ];  // update original coords
			isvab++;
			ievab++;

		}
	}

/*	if( ielem == 33 )
	{
		printf( "getDispUpdated\n" );
		for( int i = 0; i < elem->nnode; i++ )
		{
			printf( "x[ %d ]\t", i );
			for( int j = 0; j < domain->ndime; j++ )
				printf( "%g\t", elcod0[ i ][ j ] );
			printf( "\n" );
		}
		printf( "\n" );
	} */
//	lohi[ 0 ] = domain->xlo;
//	lohi[ 1 ] = domain->xhi;
//	lohi[ 2 ] = domain->ylo;
//	lohi[ 3 ] = domain->yhi;
//	lohi[ 4 ] = domain->xy;
//	for( int idime = 0; idime < domain->ndime; idime++ )
//		Cords[ 1 ][ idime ] = elcod0[ 0 ][ idime ];
	for( int jnode = 1; jnode <  elem->nnode; jnode++ )
	{
		Xj = elcod0[ 0 ][ 0 ];
		Yj = elcod0[ 0 ][ 1 ];
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			Cords[ 0 ][ idime ] = elcod0[ jnode ][ idime ];
			Cords[ 1 ][ idime ] = elcod0[ 0 ][ idime ];
		}
//		cout << "here" << '\n';	
		Distance( Cords, deform, r, theta, cross_pbc ); // find coords in the current frame
		Xj += ( *r ) * cos( *theta );
		Yj += ( *r ) * sin( *theta );
		elcod0[ jnode ][ 0 ] = Xj;
		elcod0[ jnode ][ 1 ] = Yj;
	}	
/*	if( ielem == 33 )
	{
		printf( "getDispElem\n" );
		printf( "xy = %g\n", lohi[ 1 ][ 4 ] );
		for( int i = 0; i < elem->nnode; i++ )
		{
			printf( "x[ %d ]\t", i );
			for( int j = 0; j < domain->ndime; j++ )
				printf( "%g\t", elcod0[ i ][ j ] );
			printf( "\n" );
		}
		printf( "\n" );
	} */
//	if( ielem == 28 ){
//		cout << "xy new frame\n";
//		cout << lcord0[ 2 ][ 0 ] << '\t' << lcord0[ 2 ][ 1 ] << '\n';
//	}
	// --- proper xyz respecting pbc
	GetNodalCords( elcod, ielem, 1, 0 );
/*	if( ielem == 33 )
	{
		printf( "getDispInit\n" );
		for( int i = 0; i < elem->nnode; i++ )
		{
			printf( "x[ %d ]\t", i );
			for( int j = 0; j < domain->ndime; j++ )
				printf( "%g\t", elcod[ i ][ j ] );
			printf( "\n" );
		}
		printf( "\n" );
	} */
//	if( ielem == 28 ){
//		cout << "xy init frame proper\n";
//		cout << domain->xy0 << '\n';
//		cout << lcord[ 2 ][ 0 ] << '\t' << lcord[ 2 ][ 1 ] << '\n';
//	}
	ievab = 0;
	for( int inode = 0; inode < elem->nnode; inode++ )
	{
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			dldis[ ievab ] = elcod0[ inode ][ idime ] - elcod[ inode ][ idime ];
			ievab++;
		}
	}
//	if( ielem == 28 )
//		cout << "jnode = " << dldis[ 4 ] << '\n';	
//	exit( EXIT_FAILURE );
}
//-----------------------------------------------------------------------
void Force::GetNodalCords( double **elcod,  unsigned int ielem, int pbc, int initial )
{
// --- returns elcod: nodal coordinates of the element

//	unsigned int nodei;
//	double **Cords; // array to store xyz of each pair of particles
//	double *r = new double; //distance between two nodes and angle of the connecting line
//	double *theta = new double;
//	int *cross_pbc = new int; // cross PBC?
//	double lohi[ 5 ] = { domain->xlo, domain->xhi,\
			      domain->ylo, domain->yhi, domain->xy };
//	lohi[ 0 ] = domain->xlo;
//	lohi[ 1 ] = domain->xhi;
//	lohi[ 2 ] = domain->ylo;
//	lohi[ 3 ] = domain->yhi;
//	lohi[ 4 ] = domain->xy;
//	coord = node->coord; // ptr to coord
//	double xj, yj; // xyz
//	memObj->doubleMat( Cords, 2, domain->ndime );

	// initial or current frame?
//	if( initial )
//	{
//		lohi[ 0 ] = domain->xlo0;
//		lohi[ 1 ] = domain->xhi0;
//		lohi[ 2 ] = domain->ylo0;
//		lohi[ 3 ] = domain->yhi0;
//		lohi[ 4 ] = domain->xy0;
//		coord = node->cordi;
//	}
//	for( int idime = 0; idime < domain->ndime; idime++ )
//		Cords[ 1 ][ idime ] = elcod[ 0 ][ idime ];
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
			Distance( Cords, initial, r, theta, cross_pbc ); // find the true distance
			Xj += ( *r ) * cos( *theta );
			Yj += ( *r ) * sin( *theta );
			elcod[ inode ][ 0 ] = Xj;
			elcod[ inode ][ 1 ] = Yj;
		}
	}
}

void Force::Distance(  double **coord,  int initial, double *r, double *theta, int *cross_pbc )
{
	// --- get s0, s1
//	double **dimensioless_cords; // declare an array to hold s0, s1 
//	double s0i, s1i, s0j, s1j, s0_ij, s1_ij;
//	memObj->doubleMat( dimensioless_cords, 2, domain->ndime ); // allocate memory

//	mpObj = new Wrap( box_dim ); // create Wrap obj
	mpObj[ initial ]->GetDimensionlessCords( 2, coord, dimensioless_cords ); // output dimension...

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
	mpObj[ initial ]->GetCords( 2, dimensioless_cords, coord ); // --- update coord

	xi = coord[ 0 ][ 0 ];
	yi = coord[ 0 ][ 1 ];
        xj = coord[ 1 ][ 0 ];
        yj = coord[ 1 ][ 1 ];
        xij = xi - xj;
        yij = yi - yj;
        rsq = xij * xij + yij * yij;
        *r = sqrt( rsq );
	*theta = atan2( yij, xij );

//	delete mpObj;	
}

