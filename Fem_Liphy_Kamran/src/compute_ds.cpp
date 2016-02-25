#include "fem.h"
#include "compute_ds.h"


Compute_ds::Compute_ds( Fem *femObj, FILE *dump_ds ):
kount( 0 ), init( 0 ), ds( 0.0 ), duration( 0 ),
ths( femObj->ths ), outFile( dump_ds )
{}
//-----------------------------------------------------------
Compute_ds::~Compute_ds()
{
	fclose( outFile );
}
//-----------------------------------------------------------
void Compute_ds::getAvalanch()
{
	double sdr;
	if( kount == 0 )
		a = ths->cSxyObj->GetSxy();
	else if( kount == 1 )
		b = ths->cSxyObj->GetSxy();
	else
	{
		c = ths->cSxyObj->GetSxy();
		sdr = 0.5 * ( c - a ); //--- drivative
		if( sdr > 0.0 ) //--- elastic loading
			 init = kount + 1; //--- initialize init
		else //--- avalanche startes!		
		{
			ds += sdr;
			duration++;
		}
		if( init - kount == 1 && duration != 0 ) //--- avalanche ends!
		{
			//--- output t, s
			fprintf( outFile, "%d\t%8.7e\n", duration, - ds );
			fflush( outFile );
//			printf( "avalanch: ds = %g\n", -ds );
			ds = 0.0;
			duration = 0;
		}
		a = b;
		b = c;
	}
	kount++;
}
