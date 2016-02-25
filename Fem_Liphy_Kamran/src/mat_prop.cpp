#include <iostream>
#include <stdio.h>
#include <cstring>
#include "stdlib.h"
#include "mat_prop.h"

//---------------------------------------------------------
MatProp::MatProp( char* dataFile, Fem *femObj ):
MAXSTRLEN( 100 ), MAXARG( 10 )
{
	infile = fopen( dataFile, "r" );
	line = new char[ MAXSTRLEN ]; // command string
	args = new char*[ MAXARG ]; // string for a parsed command
	for( int i = 0; i < MAXARG; i++ ) 
		args[ i ] = new char[ MAXSTRLEN ];
	femPtr = femObj;	
}
//---------------------------------------------------------
MatProp::~MatProp()
{
	delete [] line;
	delete [] args;
//	delete [] femPtr->elem->kmodu;
//	delete [] femPtr->elem->gmodu;
//	delete [] femPtr->elem->uniax;
//	delete [] femPtr->elem->strny;
}
//---------------------------------------------------------
void MatProp::process()
{
	//junk
	fgets( line, MAXSTRLEN, infile );
	// loop starts: read material properties
	unsigned int nelem = femPtr->elem->nelem;
	femPtr->elem->kmodu = new double[ nelem ];
	femPtr->elem->gmodu = new double[ nelem ];
	femPtr->elem->uniax = new double[ nelem ];
	femPtr->elem->strny = new double[ nelem ];
	for( int ielem = 0; ielem < nelem; ielem++ ){
		fgets( line, MAXSTRLEN, infile );
		parse();
		femPtr->elem->kmodu[ ielem ] = atof( args[ 1 ] );
		femPtr->elem->gmodu[ ielem ] = atof( args[ 2 ] );
		femPtr->elem->uniax[ ielem ] = atof( args[ 3 ] );
		femPtr->elem->strny[ ielem ] = atof( args[ 4 ] );
		}
	fclose( infile );
	femPtr->elem->meanGmodu = femPtr->elem->gmodu[ 0 ];
	femPtr->elem->meanKmodu = femPtr->elem->kmodu[ 0 ];	
	femPtr->elem->meanRecoveryStrain = femPtr->elem->strny[ 0 ];
}

//---------------------------------------------------------
//--- parse the string
void MatProp::parse()
{
	char *ptr;
	ptr = strtok( line, " " );
	int narg = 0;
	while( ptr != NULL ){
		strcpy( args[ narg ], ptr );
		ptr = strtok( NULL, " " );
		narg++;
	}
}
