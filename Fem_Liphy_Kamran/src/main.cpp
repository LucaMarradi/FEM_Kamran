/**
 * @file    example_action.h
 * @author  Kamran Karimi (kamran.karimi1362@gmail.com)
 * @author  Alexandre Nicolas (email)
 * @date    September, 2008
 * @version 0.1
 * @brief   This is the main file.
 * @note 		- The standard adopted to name the variables is the
 *         	  "Camel Case"
 *         	- Maybe could be interesting using the FEM parallel library
 *         	  deal_II (https://www.dealii.org)
 */

//#include "Python.h"
#include "fem.h"
#include <ctime>
#include <iostream.h>

using std::cout;

int main()
{
	clock_t begin = clock();
	Fem *femPtr = new Fem;
	femPtr->input->file();
	femPtr->Loop();
	clock_t end = clock();
	double elapsedSecs = double ( end - begin ) / CLOCKS_PER_SEC;
	fprintf( femPtr->logFile, "Loop time of %g s for %d steps with %d nodes\n", elapsedSecs, femPtr->run->nstep, femPtr->node->npoin );
	delete femPtr;
	return 0;
}
