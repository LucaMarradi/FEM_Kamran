#include "compute_sxy.h"
#include <stdio.h>
#include <math.h>

ComputeSxy::ComputeSxy()
{}
ComputeSxy::~ComputeSxy()
{}

void ComputeSxy::Init(  double **strain,  double **stress,  double *vol,
		        unsigned int nelem0,  unsigned int ngaus0 )
{
	s = stress;
	eps =  strain; 
	area = vol; 
	nelem = nelem0;
	ngaus = ngaus0;
}
double ComputeSxy::GetSxy()
{
	sxyAvg = 0.0;
	sum = 0.0;
	totArea = 0.0;
	kgaus = 0;

	for( unsigned int ielem = 0; ielem < nelem; ielem++ )
	{
		for( int igaus = 0; igaus < ngaus; igaus++ )
		{
			sum += s[ kgaus ][ 2 ] * area[ kgaus ];
			totArea += area[ kgaus ];
			kgaus++;
		}
	}
	sxyAvg = sum;
	if( totArea != 0.0 )
		sxyAvg /= totArea;
	return sxyAvg;
}

double ComputeSxy::GetP()
{
	sxyAvg = 0.0;
	sum = 0.0;
	totArea = 0.0;
	kgaus = 0;

	for( unsigned int ielem = 0; ielem < nelem; ielem++ )
	{
		for( int igaus = 0; igaus < ngaus; igaus++ )
		{
				sum += ( s[ kgaus ][ 0 ] + s[ kgaus ][ 1 ] ) * area[ kgaus ] * ( - 0.5 );
				totArea += area[ kgaus ];
			kgaus++;
		}
	}
	sxyAvg = sum;
	if( totArea != 0.0 )
		sxyAvg /= totArea;
	return sxyAvg;
}

double ComputeSxy::GetEv()
{

	sxyAvg = 0.0;
	sum = 0.0;
	totArea = 0.0;
	kgaus = 0;
	for( unsigned int ielem = 0; ielem < nelem; ielem++ )
	{
		for( int igaus = 0; igaus < ngaus; igaus++ )
		{
				sum += ( eps[ kgaus ][ 0 ] + eps[ kgaus ][ 1 ] ) * area[ kgaus ];
				totArea += area[ kgaus ];
			kgaus++;
		}
	}
	sxyAvg = sum;
	if( totArea != 0.0 )
		sxyAvg /= totArea;
	return sxyAvg;
}
