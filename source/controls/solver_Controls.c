#include"solver_Controls.h"

//calculates the max courant number
void courantNumber()
{
	uMax=0;
	vMax=0;
	wMax=0;

	for  (k=1; k<=nz; k++)
	 for (j=1; j<=ny; j++)
	  for(i=1; i<=nx; i++)
	{
		if(u[i][j][k][l]>uMax)	uMax=u[i][j][k][l];
		if(v[i][j][k][l]>vMax)	vMax=v[i][j][k][l];			//finding out max velocity
		if(w[i][j][k][l]>wMax)	wMax=w[i][j][k][l];
	}

	maxCourantNumber = (uMax/dx + vMax/dy + wMax/dz)*dt;

	printf("Max courant number: %lf\n",maxCourantNumber);


}

//Increments the time variables
void incrementTime()
{

	if(maxCourantNumber>allowableCourantNumber)
	{
		dt = allowableCourantNumber/(4*(uMax/dx + vMax/dy + wMax/dz));
		courantNumberControl=0;
	}

	if(courantNumberControl==10)						dt = deltaT;


	simulationTime = simulationTime + dt;	//time increment

	t=t+1;			//output file name parameter
	courantNumberControl = courantNumberControl+1;

}



//shifts the values to l=0
void shift()
{

	for  (k=0; k<(zNumberOfCells+5); k++)
	 for (j=0; j<(yNumberOfCells+5); j++)
	  for(i=0; i<(xNumberOfCells+5); i++)
	{
	  u[i][j][k][0]=u[i][j][k][1];
	  v[i][j][k][0]=v[i][j][k][1];
	  w[i][j][k][0]=w[i][j][k][1];
	  p[i][j][k][0]=p[i][j][k][1];			//copying values to l=0
	  pp[i][j][k][0]=pp[i][j][k][1];

	  pp[i][j][k][1]=0;	                          //equating values in l=1 to zero

	}

}

//Starts CPU clock
void startCPUClock()
{

	t=1;			//output file starts from 1st time step

	simulationTime=dt;		//time counting starts from dt

	assert((start = clock())!=-1);			//starting real time

}

//Calculation of the spacial and temporal step sizes
void stepSize()
{

	dx = xDomainLength/xNumberOfCells;
	dy = yDomainLength/yNumberOfCells;
	dz = zDomainLength/zNumberOfCells;
	deltaT=(endTime-startTime)/nTimeSteps;  //initial time step

	dt=deltaT;  				//equating initial time step to time step

	nx=xNumberOfCells + 2;
	ny=yNumberOfCells + 2;
	nz=zNumberOfCells + 2;
	
	nxm = nx - 1;
	nym = ny - 1;
	nzm = nz - 1;


}


