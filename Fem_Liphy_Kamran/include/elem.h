#ifndef ELEM_H
#define ELEM_H
#include "memory.h"
#include "moduli.h" // for material parameters
#include "rand_numbers.h"
#include "domain.h"
#include "matrix.h"


class Fem;
class GaussCoord;
class Elem
{
	public:
		Elem( Fem * );
		~Elem();

		int nnode, ngaus, nstre;
		unsigned int nelem, ntype, nevab;
		double meanGmodu, meanKmodu, meanRecoveryStrain;
		unsigned int *ltype;
//		double *mass;
		unsigned int **lnods;
		double *kmodu, *gmodu, *uniax, *strny;
		double **stran, **strni, **strst, **strsi, *epstn, **stranRate, **plasticStrain, **plasticStrain0, *vol;
		double **cordg;
		int *actField;

		void Init();
		void set(  unsigned int ); // sets shear modulus etc.		
		void GAUSSQ( double **, double * ); // quadrature points
		void Jacob2( double **, double *,  double,  double,  double **,  unsigned int );
		void Shape2( double *, double **,  double,  double );
		void Blarge( double **,  double ** );
		void ksiTOxy( double *,  double *,  double ** );
	private:
		Fem *femPtr;
		Domain *domain;
		Memory *memPtr;
		Matrix *matPtr;
		inline double det(  double ** );
		inline void inv( double **,  double ** );
		RandNumbers *randOBJ; // for random numbers	
		double *shape, **deriv;
		double **XJACM, **XJACI;
		double minE;
		double p;
		double S2, T2, SS, TT, ST, SST, STT, ST2;
		int ngash, mgash;
		double determinant;
		GaussCoord *gcObj;
		
};

#endif
