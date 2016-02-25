#ifndef TRACER_H
#define TRACER_H
#include "memory.h"
#include "node.h"
#include "elem.h"
#include "domain.h"
#include "mapping.h"
#include "time_step.h"
#include "fix_deform.h"

class Fem;
class Tracer
{
	public:
		Tracer( Fem * );
		~Tracer();	
		void TracerUpdate();
		void Init();
		void ApplyHOMO();
		void UpdateCords();
		double **coord, *veloi, *velot, *dispt, *tdisp;
		bool msd;
	private:
		void GetDisp( double *,  unsigned int );
		void GetNodalCords( double **,  unsigned int, int pbc = 1, int initial = 0 ); // nodal coordinates
		void GetVeloc( double *, double **,  unsigned int );
		void Distance(  double **,  int, double *, double *, int * ); // find true distance
		void FindElem( unsigned int, unsigned int * );
		void GetExispEtasp( double *, double **, double *, double * );
		void Nod2Gas( double *, double **, double, double );
		void Nod2Gas( double *, double *, double, double );
		Fem *femPtr;
		Memory *memObj;
		Node *node;
		Elem *elem;
		Domain *domain;
		TimeStep *ts;   // time step
		FixDeform *fd;  // shear 
		

		Wrap *mpObj[ 2 ];	

		double **cordPtr[ 2 ];
		unsigned int kgaus, // gauss point counter
			     elemi, //--- elem. id	
			     nodei,
			     isvab,
			     *indx; //--- elem. indx
		double **posgp, // position of gauss points
		       **elcod, **elcod0, **veloc, // nodal xyz
		       **Cords,
		       **lohi,	
		       **dimensioless_cords, // declare an array to hold s0, s1	
		       **deriv,
		       **cartd; //shape funcs derivatives
		double  exisp, etasp, // gauss points positions
			Xj, Yj,
			xj, yj, xi, yi, xij, yij,
			s0i, s1i, s0j, s1j, s0_ij, s1_ij,
			rsq;
		 double *weigp, // quadrature points: weight
			*shape, *velocity;
		 double /**djacb, //determinant of jacobian*/
			*r, //distance between two nodes and angle of the connecting line
			*theta;
		 int *cross_pbc; // cross PBC?

};

#endif
