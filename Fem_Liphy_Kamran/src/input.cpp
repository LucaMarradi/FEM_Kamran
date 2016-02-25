//#include "Python.h"
#include <cstring>
#include <cassert>
#include <iostream>
#include "stdlib.h"
#include "input.h"
#include "fem.h"
#include "read_data.h"
#include "mat_prop.h"

using std::cout;

#define MAXLINE	2048
#define MAXARGS 10
#define MAXARGLEN 100
#define MAXDUMP 100

// constructor
Input::Input( Fem *fem ):
narg( 0 )
{
	infile = fem->infile;	
	line = new char[ MAXLINE ]; 
	copy = new char[ MAXLINE ];
	allocate_memory();	
	femObj = fem;
	
}
// destructor
Input::~Input()
{
	delete line;
	delete [] arg;
}

/*-----------------------------------------------------------------
read lammps script
-------------------------------------------------------------------*/
void Input::file()
{
	int n;
	while( 1 ) //--- reads the input file
	{
		if ( fgets( line, MAXLINE, infile ) == NULL ) n = 0;
		else n = strlen( line ) + 1;
		if ( n == 0 ) break;
		//--- echo the command
		fprintf( femObj->logFile, line, "\n" );
		parse();
		if ( command == NULL ) continue;
		execute_command();
	}
}

//thermo_style
void Input::dump()
{
	if( femObj->dp->nDump == 0 ){
		femObj->dp->n = new int[ MAXDUMP ];
		femObj->dp->narg = new int[ MAXDUMP ];
		femObj->dp->output = new char*[ MAXDUMP ];
		femObj->dp->args = new char**[ MAXDUMP ];
	}	 
	femObj->dp->nDump++;
	int dump_id = atoi( arg[ 0 ] );	

	femObj->dp->n[ dump_id ] = atoi( arg[ 1 ] ); //nevery

	femObj->dp->output[ dump_id ] = new char[ MAXARGLEN ];
	strcpy( femObj->dp->output[ dump_id ], arg[ 2 ] ); //output-file name
 
	// copy args
	assert ( narg - 3 > 0 );
	femObj->dp->narg[ dump_id ] = narg - 3;
	femObj->dp->args[ dump_id ] = new char*[ narg - 3 ];
	for( int iarg = 0; iarg < narg - 3; iarg++ ){
		femObj->dp->args[ dump_id ][ iarg ] = new char[ MAXARGLEN ];
		strcpy( femObj->dp->args[ dump_id ][ iarg ], arg[ iarg + 3 ] );
	}
}

//thermo_style
void Input::thermo_style()
{
	femObj->ths->n = atoi( arg[ 0 ] );
	assert ( narg - 1 > 0 );
	femObj->ths->narg = narg - 1;
	// copy args
	femObj->ths->args = new char*[ narg - 1 ];
	for( int iarg = 0; iarg < narg - 1; iarg++ ){
		femObj->ths->args[ iarg ] = new char[ MAXARGLEN ];
		strcpy( femObj->ths->args[ iarg ], arg[ iarg + 1 ] );
	}
}

//run
void Input::msd()
{
	femObj->tr->msd = 1;
}
//run
void Input::min()
{
	femObj->qs->min = 1;
	femObj->qs->ftol = atof( arg[ 0 ] );
}
//run
void Input::fix_force(){
	femObj->ff->fix_force = 1;
	femObj->node->fx = atof( arg[ 0 ] );
	femObj->node->fy = atof( arg[ 1 ] );
}
//run
void Input::fix_veloc(){
	femObj->fvl->fix_vel = atoi( arg[ 0 ] );
}
//run
void Input::run(){
	femObj->run->nstep = atoi( arg[ 0 ] );
}

//fix viscous
void Input::fix_viscous(){
	femObj->fv->cdrag0 = atof( arg[ 0 ] );
	femObj->fv->cdrag1 = atof( arg[ 1 ] );
}

//fix deform
void Input::fix_deform(){
	femObj->fd->shearRATE = atof( arg[ 0 ] );
}

//time step
void Input::time_step(){
	femObj->ts->dt = atof( arg[ 0 ] );
}

//restart
void Input::restart()
{
	femObj->rs->n = atoi( arg[ 0 ] );
}

//material properties
void Input::mat_prop(){
	fprintf( femObj->logFile,  "Reading material properties ...\n") ;
	MatProp *mpObj = new MatProp( arg[ 0 ], femObj );
	mpObj->process();
	delete mpObj;
}

//data file
void Input::read_data(){
	fprintf( femObj->logFile,  "Reading data file ...\n") ;
	ReadData *rdObj = new ReadData( arg[ 0 ], femObj );
	rdObj->process();
	delete rdObj;
}
/*-----------------------------------------------------------------
process a single parsed command
-------------------------------------------------------------------*/
void Input::execute_command()
{
	if ( !strcmp( command, "read_data" ) ) read_data();
	else if ( !strcmp( command, "mat_prop" ) ) mat_prop();
	else if ( !strcmp( command, "timestep" ) ) time_step();
	else if ( !strcmp( command, "restart" ) ) restart();
	else if ( !strcmp( command, "fix_deform" ) ) fix_deform();
	else if ( !strcmp( command, "fix_force" ) ) fix_force();
	else if ( !strcmp( command, "fix_veloc" ) ) fix_veloc();
	else if ( !strcmp( command, "min" ) ) min();
	else if ( !strcmp( command, "msd" ) ) msd();
	else if ( !strcmp( command, "fix_viscous" ) ) fix_viscous();
	else if ( !strcmp( command, "run" ) ) run();
	else if ( !strcmp( command, "thermo_style" ) ) thermo_style();
	else if ( !strcmp( command, "dump" ) ) dump();
}
/*-----------------------------------------------------------------
parse copy of the command line
-------------------------------------------------------------------*/
void Input::parse()
{
	strcpy( copy, line );
	char *ptr = copy;
	while( *ptr )
	{
		if( *ptr == '#' )
		{
			*ptr = '\0';
			break;
		} 
		ptr++;
	}
	// command = 1st arg
	command = strtok( copy, " \t\n\r\f" );
	if( command == NULL )
		return;
	// next args
	
	narg = 0;
	while( 1 )
	{
		arg[ narg ] = strtok( NULL, " \t\n\r\f");
		if ( arg[ narg ] ) narg++;
		else break;
	}

}

// utility function for memory allocation
void Input::allocate_memory()
{
	arg = new char*[ MAXARGS ];
	for( int i = 0; i < MAXARGS; i++ )
	{
		arg[ i ] = new char[ 20 ];
	}
}


//int main()
//{
//	Fem *femobj = new Fem;	
//	Input *inputObj = new Input;//( femobj );
//	inputObj->file();
//	delete femobj;
//	delete inputObj;

//}
