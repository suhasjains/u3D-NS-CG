#include"pv_Coupling.h"
#include"../velocity/velocity.h"
#include"../update/update.h"


//SIMPLE algorithm loop
void simple()
{

	printf("Starting SIMPLE loop\n");
	
	
	for (simplerVar = 1;simplerVar <= nSimpleLoops;simplerVar++)
 {

		nLoops();	//Initializes no of velocity  and pressure inner loops
		
		boundary_pressure(p);
		
		uEqnCoeff();
		
		velocityBoundaryConditions();
		
		uEqn();

		vEqn();		

		wEqn();

		pressure_correction_east();
		
		pressure_correction_north();
		
		pressure_correction_top();
		
		pressure_correction_source();
		
		pCorrEqn();		//finding pressure correction based on initial velocities

		boundary_pressure(pp);
				
		correctFaceFlux();
	
		update_Pressure_Velocity();       //updating pressure
		
		printf("Global mass residual%e\n",maxMassResidual);
		
		if(fabs(maxMassResidual)>accuracy)		nSimpleLoops = nSimpleLoops + 1;

		//printf("u = %e\n",u[10][10][4][l]);
		
 }

	printf("No of SIMPLE loops: %d\n",nSimpleLoops);
	printf("No of pressure loops: %d\n",nPressureLoops);
	printf("No of U Velocity loops: %d\n",uVelocityLoops);
	printf("No of V Velocity loops: %d\n",vVelocityLoops);
	printf("No of W Velocity loops: %d\n",wVelocityLoops);
	printf("max residual %e\n",std::max(std::max(std::max(uMaxRes,vMaxRes),wMaxRes),pMaxRes));

	check = check + nSimpleLoops;
	printf("Total loops= %d\n",check);	
}

