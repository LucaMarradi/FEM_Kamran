#ifndef FEM_H
#define FEM_H
#include "node.h"
#include "stdio.h"
#include "elem.h"
#include "domain.h"
#include "time_step.h"
#include "fix_deform.h"
#include "fix_force.h"
#include "fix_veloc.h"
#include "quasi_static.h"
#include "fix_viscous.h"
#include "run.h"
#include "thermo_style.h"
#include "dump.h"
#include "input.h"
#include "force.h"
#include "mapping.h"
#include "restart.h"
#include "compute_ds.h"
#include "tracer.h"

class Fem
{
	public:
		Fem();
		~Fem();

		FILE *infile; // standard input
		FILE *logFile; // standard input
		Node *node; // node-based quantities
		Elem *elem; // element-based variables
		Domain *domain; // simulation box
		TimeStep *ts;   // time step
		FixDeform *fd;  // shear
		FixVeloc *fvl;  // shear
		FixForce *ff;
		QuasiStatic *qs;
		FixViscous *fv;  // viscosity
		Run *run;
		ThermoStyle *ths;
		Dump *dp;
		Input *input; // loading input from hard disk
		Force *force;  // inter-particle forces
		Restart *rs;
		Compute_ds *cds;
		Tracer *tr;
		void Loop();

	private:
		void Verlet( double *, double *, double, double );
		void ApplyHOMO();
		void Switch();
		void Init();
		void SetVeloc();
		void SetForce();
		bool Iterate( unsigned int * );
		double dGamma, dtfm;
};

#endif
