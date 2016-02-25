#include <iostream>
#include "fem.h"
#include "stdlib.h"
#include <cmath>
#include <cassert>

using std::cout;

Fem::Fem()
{
	infile = stdin;
	logFile = stdout; //fopen( "log.fem", "w" );	
	domain = new Domain();
	node   = new Node();
	srand( time( NULL ) );
	elem   = new Elem( this );
	ts	   = new TimeStep();
	fd	   = new FixDeform();
	fvl	   = new FixVeloc();
	ff	   = new FixForce();
	qs	   = new QuasiStatic();
	fv	   = new FixViscous();
	run    = new Run();
	ths	   = new ThermoStyle( this );
	dp	   = new Dump( this );
	input  = new Input( this );
	force  = new Force( this );
	rs     = new Restart( this );
	tr     = new Tracer( this );
	cds     = new Compute_ds( this, fopen( "ds.xyz", "w" ) );
}

Fem::~Fem()
{
	delete force;
	delete node;
	delete elem;
	delete domain;
	delete ts;
	delete fd;
	delete fv;
	delete ff;
	delete fvl;
	delete qs;
	delete run;
	delete ths;
	delete dp;
	delete input;
	delete rs;
	delete tr;
	delete cds;
//	fclose( logFile );
}
//------------------------------------------------------
//--- perfom a loop for numerical integration
//------------------------------------------------------
void Fem::Loop()
{
	unsigned int ntime = run->nstep;
	unsigned int *counter = new unsigned int; //--- min. steps 
	*counter = 0;
	//--- computes nodal internal forces at x( n )
	Init(); // initialize state variables			
	if( fvl->fix_vel )
		SetVeloc(); // set velocities for "wall" nodes
	if( ff->fix_force )
		SetForce(); // set velocities for "wall" nodes
	// --- std output!
	ths->Process( ts->itime0 );
	// --- output data
	for( int iDump = 0; iDump < dp->nDump; iDump++ )
		dp->Process( 0, iDump );
	for( int itime = 0; itime < ntime; itime++ ){ // time integration starts
		if( ! fvl->fix_vel && ! ff->fix_force ) 
			ApplyHOMO(); // --- apply homogeneous deformation
		while( Iterate( counter ) ) //--- q-st
		{
			Verlet( node->veloi, node->velot, ts->dt, ts->dt ); //output node->veloi, node->dldis: nonaffine piece
		// --- computes nodal internal forces at x( n + 1 )
			force->ComputeForce( *counter ); //--- if counter == 1-> update boundaries 
		// --- calculates v( n + 1 )
			Verlet( node->velot, node->veloi, ts->dt, 0.0 ); //output node->velot
			Switch(); //--- switch state variables
		}
		if( tr->msd ) //--- tracer particle
		{
			tr->TracerUpdate();// --- update tr->dispt ( use tr->tdisp ) 
			tr->UpdateCords(); //--- output tr->coord
		}

		if(  rs->n != 0 ) //--- write restart
		{
			if( ( itime + 1 ) % rs->n == 0 ) 
				rs->output( ts->itime0 + itime + 1 );
		}
		// --- std output!
		if ( ( itime + 1 ) % ths->n == 0 ){
			ths->Process( ts->itime0 + itime + 1 );
		}
		// --- output data
		for( int iDump = 0; iDump < dp->nDump; iDump++ ){
			if( ( itime + 1 ) % dp->n[ iDump ] == 0 )
				dp->Process( 0 + itime + 1, iDump );
		}	
		cds->getAvalanch();
	} // end of time loop
	delete counter;
}
//------------------------------------------------------
bool Fem::Iterate( unsigned int *counter )
{
        ( *counter )++;
        //--- residual force
        double fnorm = 0.0;
        for( unsigned int isvab = 0; isvab < node->nsvab; isvab++ )
                fnorm += node->fintl[ isvab ] * node->fintl[ isvab ];
//      fprintf( logFile, "%d\t%e\n", counter, fnorm / node->nsvab );
        if( *counter == 1 ) //--- run at least once!
                return 1;
	if( *counter == 2 && ( ! qs->min ) ) //--- no minimization
	{
                *counter = 0;
		return 0;
	}
        //--- meet the min. criterion 
	assert( qs->min );
        if( ( fnorm / node->nsvab < qs->ftol ) ||
                ( *counter == 8 && fnorm / node->nsvab < 10.0 * qs->ftol ) ) //--- avoid further min. in the elastic regime (underdamped) 
        {
                *counter = 0;
                return 0; //--- success
        }
        else
                return 1;
}
//------------------------------------------------------ 
void Fem::SetForce()
{
        double y = 0.0;
        unsigned int isvab = 0, count = 0;
	for( int isvab = 0; isvab < node->nsvab; isvab++ )
		node->extForce[ isvab ] = 0.0; //--- initialize external force

	//--- initialize fixID
        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{
                for( int idime = 0; idime < domain->ndime; idime++ )
                        node->fixID[ isvab ] = 0;
                        node->rigidID[ isvab ] = 0;
                        isvab++;
        }

	// --- assign fixID and extForce
        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{
		isvab = ipoin * domain->ndime;
                y = node->coord[ ipoin ][ 1 ];
                if( fabs( y - domain->yhi ) < 1.0e-08 )
		{
			node->effMass += node->mass[ ipoin ]; //--- effective mass
                        node->rigidID[ isvab ] = 1; // --- nodes with fext
			node->extForce[ isvab ] = node->fx;
                        node->fixID[ isvab + 1 ] = 1; //--- constrained in y
                        node->presVeloc[ isvab + 1 ] = 0.0; // --- vy
//                        node->rigidID[ isvab + 1 ] = 1; // --- nodes with fext
//			node->extForce[ isvab + 1 ] = node->fy;
                }
                if( fabs( y - domain->ylo ) < 1.0e-08 )
		{
                        node->fixID[ isvab ] = 1;
                        node->fixID[ isvab + 1 ] = 1;
                        node->presVeloc[ isvab ] = 0.0;
                        node->presVeloc[ isvab + 1 ] = 0.0;
                }
        }
}
//------------------------------------------------------
void Fem::SetVeloc() 
{                       
        double y = 0.0;        
        unsigned int isvab = 0, count = 0;
        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
        {
                for( int idime = 0; idime < domain->ndime; idime++ )
                        node->fixID[ isvab ] = 0;
                        isvab++;
        }

        for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
        {
                isvab = ipoin * domain->ndime;
                y = node->coord[ ipoin ][ 1 ];
                if( fabs( y - domain->yhi ) < 1.0e-08 )
                {
                        node->fixID[ isvab ] = 1;
                        node->fixID[ isvab + 1 ] = 1;
                        node->presVeloc[ isvab ] = fd->shearRATE * ( domain->yhi - domain->ylo ); //--- unit length
                        node->presVeloc[ isvab + 1 ] = 0.0;
                }
                if( fabs( y - domain->ylo ) < 1.0e-08 )
                {
                        isvab = ipoin * domain->ndime;
                        node->fixID[ isvab ] = 1;
                        node->fixID[ isvab + 1 ] = 1;
                        node->presVeloc[ isvab ] = 0.0;
                        node->presVeloc[ isvab + 1 ] = 0.0;
                }
        }
}
//------------------------------------------------------
void Fem::Init()
{
	domain->xlo = domain->xlo0;
	domain->xhi = domain->xhi0;
	domain->ylo = domain->ylo0;
	domain->yhi = domain->yhi0;
	domain->xy = domain->xy0;
	for( int igaus = 0; igaus < elem->nelem * elem->ngaus; igaus++ ){ //strsi = strst
		for( int istre = 0; istre < elem->nstre; istre++ ){
			elem->strst[ igaus ][ istre ] = elem->strsi[ igaus ][ istre ];
			elem->stran[ igaus ][ istre ] = elem->strni[ igaus ][ istre ];
			elem->plasticStrain[ igaus ][ istre ] = elem->plasticStrain0[ igaus ][ istre ];
		}
	}
	unsigned int isvab = 0;
	for( int ipoin = 0; ipoin < node->npoin; ipoin++ ){ // cordi = coord
		for( int idime = 0; idime < domain->ndime; idime++ ){
	                node->coord[ ipoin ][ idime ] = node->cordi[ ipoin ][ idime ];
	                node->coordUnWrapped[ ipoin ][ idime ] = node->cordiUnWrapped[ ipoin ][ idime ];
			node->dldis[ isvab ] = 0.0;
			isvab++;
		}
	}
	
	force->Init(); // --- initialize "force" member
	if( tr->msd )
		tr->Init(); //--- init "tracer"
}
//------------------------------------------------------
void Fem::Switch()
{
	for( int igaus = 0; igaus < elem->nelem * elem->ngaus; igaus++ ) //strsi = strst
	{
		for( int istre = 0; istre < elem->nstre; istre++ )
		{
			elem->strsi[ igaus ][ istre ] = elem->strst[ igaus ][ istre ];
			elem->strni[ igaus ][ istre ] = elem->stran[ igaus ][ istre ];
			elem->plasticStrain0[ igaus ][ istre ] = elem->plasticStrain[ igaus ][ istre ];
		}
	}
	unsigned int isvab = 0;
	for( int ipoin = 0; ipoin < node->npoin; ipoin++ ) // cordi = coord
	{
		for( int idime = 0; idime < domain->ndime; idime++ )
		{
			node->dispt[ isvab ] += node->dldis[ isvab ];
			tr->tdisp[ isvab ] += node->dldis[ isvab ]; //--- set zero in ApplyHomo
			node->dldis[ isvab ] = 0.0;
	                node->cordi[ ipoin ][ idime ] = node->coord[ ipoin ][ idime ];
	                node->cordiUnWrapped[ ipoin ][ idime ] = node->coordUnWrapped[ ipoin ][ idime ];
			isvab++;
		}
	}
}
//------------------------------------------------------
void Fem::ApplyHOMO()
{
	domain->xy = domain->xy0 + fd->shearRATE * ts->dt * ( domain->yhi - domain->ylo );
	// --- update displacements
	dGamma = ts->dt * fd->shearRATE;
	// --- initialize dispt
	for( int isvab = 0; isvab < node->nsvab; isvab++ )
		tr->tdisp[ isvab ] = 0.0;
	for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
		node->dldis[ ipoin * domain->ndime ] = dGamma * ( node->coord[ ipoin ][ 1 ] - domain->ylo );
	
}
//------------------------------------------------------
void Fem::Verlet( double *output, double *v, double dtf, double dtv )
{
	double effForce[ domain->ndime ];
	unsigned int isvab = 0;
	if( ff->fix_force )
	{
	//--- effective force
	for( int ipoin = 0; ipoin < node->npoin; ipoin++ )
	{
        	for( int idofn = 0; idofn < domain->ndime; idofn++ )
        	{
        	        if( node->rigidID[ isvab ] == 1 )
        	                effForce[ idofn ] += node->fintl[ isvab ];
        	        isvab++;
        	}
	}	
	}
	for( int isvab = 0; isvab < node->nsvab; isvab++ ){
		dtfm = dtf / node->mass[ isvab / domain->ndime ];
		output[ isvab ] = v[ isvab ] + 0.5 * dtfm * ( - node->fintl[ isvab ] );
		if( ( fvl->fix_vel || ff->fix_force ) && node->fixID[ isvab ] == 1 ) //--- set veloc
			output[ isvab ] = node->presVeloc[ isvab ];
		if( ff->fix_force && node->rigidID[ isvab ] == 1 ) // --- nodes with external force
			output[ isvab ] = v[ isvab ] + 0.5 * dtf * ( - effForce[ isvab % domain->ndime ] + node->extForce[ isvab ] ) / node->effMass;
		if( dtv != 0.0 )
			node->dldis[ isvab ] += dtv * output[ isvab ];
	}
}
