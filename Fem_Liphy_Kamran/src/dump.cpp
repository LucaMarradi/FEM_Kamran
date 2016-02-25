//#include "Python.h"
#include "dump.h"
#include "fem.h"
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <string>

using std::cout;
using std::string;

Dump::Dump( Fem *femObj ):
node( femObj->node ), elem( femObj->elem ), domain( femObj->domain ), 
ts( femObj->ts ), fd( femObj->fd ), femPtr( femObj ),
//gcObj( new GaussCoord( femObj ) ),
nDump( 0 ), memPtr( new Memory ) 
{}

Dump::~Dump()
{
	size_t found;
	for( int iDump = 0; iDump < nDump; iDump++ )
	{
		string stdStr( output[ iDump ] ); // --- standard str
		found = stdStr.find( '*' ); 
		if( found == std::string::npos ) // --- only close single file
			fclose( outptFile[ iDump ] );
	}
	if( nDump != 0 )
	{
		delete [] outptFile;
//		memPtr->Delete2ndMat( cordg, elem->nelem * elem->ngaus );
	}
	delete memPtr;
//	delete gcObj;
}

void insert( char *str1, size_t pos, char *str2 )
{
	string stdStr1( str1 );
	stdStr1.erase( pos, 1 );
	string stdStr2( str2 );
	stdStr1.insert( pos, stdStr2 );
	strcpy( str1, stdStr1.c_str() );
	
}
 void Dump::Init()
{
//	memPtr->doubleMat( cordg, elem->nelem * elem->ngaus, domain->ndime ); // allocate cordg	
//	gcObj->Init();
	outptFile = new FILE*[ nDump ];
}
void Dump::Process(  unsigned int itime,  int iDump )
{
	if( itime == 0 && iDump == 0 ) //--- define file array 
		Init();
	char tmp0[ 64 ], tmp1[ 64 ];
	sprintf( tmp0, "%d", itime ); // --- int -> string
	strcpy( tmp1, output[ iDump ] ); // --- assign tmp1
	string stdStr( tmp1 ); // --- standard str
	size_t found = stdStr.find( '*' ); // --- dump.*.xyz

	if( found != std::string::npos ) //--- separate files need to be opened at every call!
	{
		insert( tmp1, found, tmp0 ); //--- format tmp1
		outptFile[ iDump ] = fopen( tmp1, "w" );
	}
	else // --- single file
	{
		if( itime == 0 ) // open at first call
			outptFile[ iDump ] = fopen( output[ iDump ], "w" );
	}


	data = new double*[ narg[ iDump ] ];
	for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
	{
		if ( !strcmp( args[ iDump ][ iarg ], "vi_x" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->veloi[ ipoin * domain->ndime ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "vi_y" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->veloi[ ipoin * domain->ndime + 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "dx" ) )
		{ 
			N = node->npoin;
			if( femPtr->tr->msd )
				N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			if( femPtr->tr->msd )
			{
				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = femPtr->tr->dispt[ ipoin * domain->ndime ];	
			}
			else
			{
				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->dispt[ ipoin * domain->ndime ];
			}
		}
		if ( !strcmp( args[ iDump ][ iarg ], "dy" ) )
		{ 
			N = node->npoin;
			if( femPtr->tr->msd )
				N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			if( femPtr->tr->msd )
			{
				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = femPtr->tr->dispt[ ipoin * domain->ndime + 1 ];	
			}
			else
			{
				for( unsigned int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->dispt[ ipoin * domain->ndime + 1 ];
			}
		}
		if ( !strcmp( args[ iDump ][ iarg ], "vx" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->velot[ ipoin * domain->ndime ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "vy" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->velot[ ipoin * domain->ndime + 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "fx" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->fintl[ ipoin * domain->ndime ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "fy" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->fintl[ ipoin * domain->ndime + 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "sxy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strst[ igaus ][ 2 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "sxx" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strst[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "syy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->strst[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "exx" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->stran[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eyy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->stran[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "exy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = 0.5 * elem->stran[ igaus ][ 2 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eDoTxy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = 0.5 * elem->stranRate[ igaus ][ 2 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eDoTxx" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->stranRate[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "eDoTyy" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->stranRate[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "vol" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->vol[ igaus ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "epX-Y" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->plasticStrain[ igaus ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "epY-X" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->plasticStrain[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "epXY" ) )
		{ 
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = 0.5 * elem->plasticStrain[ igaus ][ 2 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "x" ) )
		{ 
			N = node->npoin;
			if( femPtr->tr->msd )
				N = elem->nelem * elem->ngaus ;
			data[ iarg ] = new double[ N ];
			if( femPtr->tr->msd )
			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = femPtr->tr->coord[ ipoin ][ 0 ];
			}
			else
			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->coord[ ipoin ][ 0 ];
			}
		}
		if( !strcmp( args[ iDump ][ iarg ], "xu" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->coordUnWrapped[ ipoin ][ 0 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "yu" ) )
		{ 
			N = node->npoin;
			data[ iarg ] = new double[ N ];
			for( int ipoin = 0; ipoin < N; ipoin++ )
				data[ iarg ][ ipoin ] = node->coordUnWrapped[ ipoin ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "y" ) )
		{ 
			N = node->npoin;
			if( femPtr->tr->msd )
				N = elem->nelem * elem->ngaus ;
			data[ iarg ] = new double[ N ];
			if( femPtr->tr->msd )
			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = femPtr->tr->coord[ ipoin ][ 1 ];
			}
			else
			{
				for( int ipoin = 0; ipoin < N; ipoin++ )
					data[ iarg ][ ipoin ] = node->coord[ ipoin ][ 1 ];
			}
		}
		if ( !strcmp( args[ iDump ][ iarg ], "n0" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( unsigned int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->lnods[ ielem ][ 0 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "n1" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( unsigned int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->lnods[ ielem ][ 1 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "n2" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( unsigned int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->lnods[ ielem ][ 2 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "xg" ) )
		{ 
			// computes cordg
//			gcObj->GetCoord( cordg );
			// assign data
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( unsigned int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->cordg[ igaus ][ 0 ];
		}
		if ( !strcmp( args[ iDump ][ iarg ], "yg" ) )
		{ 
//			gcObj->GetCoord( cordg );
			// assign data
			N = elem->nelem * elem->ngaus;
			data[ iarg ] = new double[ N ];
			for( unsigned int igaus = 0; igaus < N; igaus++ )
				data[ iarg ][ igaus ] = elem->cordg[ igaus ][ 1 ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "K" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->kmodu[ ielem ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "G" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->gmodu[ ielem ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "sy" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->uniax[ ielem ];
		}
		if( !strcmp( args[ iDump ][ iarg ], "Ey" ) )
		{ 
			N = elem->nelem;
			data[ iarg ] = new double[ N ];
			for( int ielem = 0; ielem < N; ielem++ )
				data[ iarg ][ ielem ] = elem->strny[ ielem ];
		}
	}
	// --- write output
	fprintf( outptFile[ iDump ], "Timestep\n%d\n", itime );
	fprintf( outptFile[ iDump ], "Number of Nodes\n%d\n", N );
	fprintf( outptFile[ iDump ], "Box Bounds\n%15.10g\t%15.10g\n%15.10g\t%15.10g\n%15.10g\n",
		domain->xlo, domain->xhi, domain->ylo, 
		domain->yhi, domain->xy );
// --- print items specified in args dictionary
	fprintf( outptFile[ iDump ], "Atoms\t" );
	for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
		fprintf( outptFile[ iDump ], "%s\t", args[ iDump ][ iarg ] );
	fprintf( outptFile[ iDump ], "\n" );
	for( int ipoin = 0; ipoin < N; ipoin++ )
	{
		fprintf( outptFile[ iDump ], "%d\t", ipoin );
		for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
			fprintf( outptFile[ iDump ], "%15.10g\t", data[ iarg ][ ipoin ] );
		fprintf( outptFile[ iDump ], "\n" );
	}
	for( int iarg = 0; iarg < narg[ iDump ]; iarg++ )
		delete [] data[ iarg ];
	delete [] data;
	if( found != string::npos ) //--- separate files need to be closed at every call!
	{
//		printf( "bug!\n" );
		fclose( outptFile[ iDump ] );
	}
}

