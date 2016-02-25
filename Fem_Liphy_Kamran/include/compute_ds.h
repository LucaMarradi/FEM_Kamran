#ifndef COMPUTE_DS_H
#define COMPUTE_DS_H
#include "thermo_style.h"
class Fem;
class Compute_ds
{
	public:
		Compute_ds( Fem *, FILE * );
		~Compute_ds();
		void getAvalanch();


	private:

		unsigned int kount, init, duration;
		double a, b, c, ds;
		FILE *outFile;
		ThermoStyle *ths;

};

#endif
