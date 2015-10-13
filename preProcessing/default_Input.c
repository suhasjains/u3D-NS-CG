#include"preProcessing.h"

//Inputs the default value that is required for simulation
void defaultInput()
{

		xDomainLength=1;	//Length of domain along x
		yDomainLength=1;
		zDomainLength=0.1;

		xNumberOfCells=50;	//Number of cells along x
		yNumberOfCells=50;
		zNumberOfCells=5;

		startTime=0;

		endTime=50;

		nTimeSteps=1000;

		writeOutputSteps=10;

		allowableCourantNumber=10e100;

		rho=0.01;
		mu=0.00001;

		uTopWallVel=1;				//highest value of y
		vTopWallVel=0;
		wTopWallVel=0;

		uBottomWallVel=0;			//lowest value of y
		vBottomWallVel=0;
		wBottomWallVel=0;

		uLeftWallVel=0;				//Lowest value of x
		vLeftWallVel=0;
		wLeftWallVel=0;

		uRightWallVel=0;			//Highest value of x 
		vRightWallVel=0;
		wRightWallVel=0;

		uFrontWallVel=0;			//Highest value of z
		vFrontWallVel=0;
		wFrontWallVel=0;

		uBackWallVel=0;				//Lowest value of z
		vBackWallVel=0;
		wBackWallVel=0;

		pUnderRelaxCoeff=0.9;		//Pressure under relaxation factor
		wUnderRelaxCoeff=0.9;
		vUnderRelaxCoeff=0.9;
		uUnderRelaxCoeff=0.9;

		xGrav=0;					//Body force value in terms of acceleration
		yGrav=0;
		zGrav=0;

		DCvalue=1;					//Deferred Correction value

		accuracy=10e-6;

		scheme=1;	//1=UQ,2=HQ,3=PQ



}

