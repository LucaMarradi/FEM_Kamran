#ifndef FORCE_H
#define FORCE_H
#include "memory.h"
#include "matrix.h"
#include "node.h"
#include "elem.h"
#include "domain.h"
#include "fix_viscous.h"
#include "fix_deform.h"
#include "time_step.h"
#include "mapping.h"

class Fem;
class Force
{
	public:
		Force( Fem * );
		~Force();
		void Init();
		void ComputeForce( int );
	private:
		void GetDisp( double *,  unsigned int, int );
		void GetVeloc( double *,  unsigned int );
		void GetNodalCords( double **,  unsigned int, int pbc = 1, int initial = 0 ); // nodal coordinates
		void Distance(  double **,  int, double *, double *, int * ); // find true distance
		inline void Modulus( double **,  double,  double,  int );
		void ReturnStress( double *, double *, double *, int *, unsigned int );
		inline double Invart(  double * );
		inline void Devia( double *,  double * );
		inline void DeviaStrain( double *,  double * );
		inline double EffectiveStrain(  double * );
		inline double Smean(  double * );
		Fem *femPtr;
		 Memory *memObj;
		 Matrix *matPtr;
		Node *node;
		Elem *elem;
		Domain *domain;
		 FixViscous *fv;  // viscosity 
		 FixDeform *fd;
		 TimeStep *ts;   // time step
		Wrap *mpObj[ 2 ];	
		double **cordPtr[ 2 ];
		unsigned int nodei, // node id
			     ievab, isvab, // dof's
			     kgaus, // gauss point counter
			     itype; // material id
		double **posgp, // position of gauss points
		       **elcod, **elcod0, // nodal xyz
		       **Dmatx, // elasticity matrix
		       **Vmatx, // elasticity matrix
		       **cartd, //shape funcs derivatives
		       **Cords,
		       **dimensioless_cords, // declare an array to hold s0, s1	
		       **dev,
		       **lohi;	
//		       **bmatx, //b-matrix
//		       **bt; // transpose of b-matrix
		double ***bmatx, ***bt; //--- b-matrix, transpose of b-matrix
		double  uniax, // yield stress/strain
			strny, // healing strain
			kmodu, gmodu, // K, G
			visc_coeff, second_visc_coeff,
			exisp, etasp, // gauss points positions
			dvolu, // area for each element
			xj, yj, xi, yi, xij, yij,
			s0i, s1i, s0j, s1j, s0_ij, s1_ij,
			Xj, Yj,
			rsq,
			exx_yy, exy, smean,
			varj2, effectiveStress, effectivePlasStrain, yfunc,
			pstre,
			randVal;
			  
		 double *weigp, // quadrature points: weight
			*dldis, *veloc, // element dipslacements
			*strag, // strain
			*dsigm, // delta sigma
			*sigmV, // delta sigma
			*sgtri, // trial stress = s_n+D * strain
			*bt_sigma, *eload, // b^t * sigma & element load vector
			*deviatoric_stress,
			*dev_trial, *dev_strsi, *dev_strag;
		 double *djacb, //determinant of jacobian
			*r, //distance between two nodes and angle of the connecting line
			*theta;
		 int *cross_pbc; // cross PBC?

};

#endif
