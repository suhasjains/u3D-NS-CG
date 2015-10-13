//my Navier Stokes solution


#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#include<math.h>
#include <algorithm>
#include <cstring>
#include"source/variables/variables.h"
#include"preProcessing/preProcessing.c"
#include"preProcessing/default_Input.c"
#include"preProcessing/dynamic_Allocate.c"
#include"preProcessing/initial_Time_Values.c"
#include"preProcessing/boundaryConditions/pressure_BC.c"
#include"preProcessing/boundaryConditions/velocity_BC.c"
#include"source/pressure/p_Correction_Coefficients.c"
#include"source/pressure/p_Correction_Equation.c"
#include"source/pv_Coupling/simple.c"
#include"source/pv_Coupling/schemes.c"
#include"source/solver/ADI_TDMA.c"
#include"source/update/correct_Face_Flux.c"
#include"source/update/update_Pressure_Velocity.c"
#include"source/velocity/u_Equation.c"
#include"source/velocity/v_Equation.c"
#include"source/velocity/w_Equation.c"
#include"source/velocity/velocity_Equation_Coefficients.c"
#include"source/controls/solver_Controls.c"
#include"postProcessing/write_Output.c"
#include"preProcessing/grid_Generation/create_Mesh.c"
#include"preProcessing/grid_Generation/boundary_Coordinates.c"
#include"preProcessing/grid_Generation/uniform_Grid.c"
#include"preProcessing/grid_Generation/write_Mesh.c"

int main()
{

	defaultInput(); 	   	//Inputting the values required

	dynamicAllocate();		//dynamically allocates the arrays

	zeroTimeValues();   	//Gives initial values to all the variables and fields at time 0

	stepSize();

	createMesh();

	startCPUClock();

	setting_ulid();


	while(simulationTime<=endTime)         	//time loop
 	{
	  initialize();

	  simple(); 		    				//use simple or simpler algorithm

	  if(t % writeOutputSteps == 0)
	  writeOutput();           		 		//Outputs values to files

	  shift();		    					//shifts the values to l=0

	  courantNumber();	    				//calculates courant number

	  incrementTime();	    				//increments time variables

	  finalize();

	}


}
