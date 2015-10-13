//Declaration of all variables


#ifndef VARIABLES_H
#define VARIABLES_H


	//INPUT OUTPUT VARIABLES
	double xDomainLength,yDomainLength,zDomainLength;
	int   xNumberOfCells,yNumberOfCells,zNumberOfCells,nx,ny,nz,n,nxm,nym,nzm;   //number of cells
	int   nTimeSteps;
	double startTime,endTime,deltaT,simulationTime,realTime,start,stop;
	char  adjustTimeStep;
	char  topBoundary,bottomBoundary,leftBoundary,rightBoundary,frontBoundary,backBoundary;
	double uTopWallVel,uBottomWallVel,uLeftWallVel,uRightWallVel,uFrontWallVel,uBackWallVel;
	double vTopWallVel,vBottomWallVel,vLeftWallVel,vRightWallVel,vFrontWallVel,vBackWallVel;
	double wTopWallVel,wBottomWallVel,wLeftWallVel,wRightWallVel,wFrontWallVel,wBackWallVel;


	//COEFFICIENTS OF VELOCITY EQUATION
	double dx,dy,dz,dt;
	double xGrav,yGrav,zGrav;
	double ***Fe,***Fn,***Ft;


	//COEFFICIENTS OF PRESSURE CORRECTION EQUATION
	double ***ae,***aw,***an,***as,***at,***ab,***ap;
	double ***su,***sv,***sw,maxMassResidual;
	double ***apu,***apv,***apw,apt;
	double uMaxRes,vMaxRes,wMaxRes,pMaxRes,***uRes,***vRes,***wRes,***pRes,uMaxResidual[2],vMaxResidual[2],wMaxResidual[2];
	int pMassResidual,nMassResidual,nsMassResidual,psMassResidual,nuMaxRes,nvMaxRes,nwMaxRes,puMaxRes,pvMaxRes,pwMaxRes;

	//VELOCITIES and PRESSURE
	double ****u,****v,****w;
	double pe,pw,pn,ps,pt,pb,ppo;
	double ***dpx,***dpy,***dpz;													//pressure gradient
	double ****pp,****p;
	double ***xCoord,***yCoord,***zCoord;
	double convf,convp,***flux;														//convective mass fluxes
	double diff;
	
	double vol,vole,voln,volt;
	
	
	double uMax,vMax,wMax;															//max velocities
	double maxCourantNumber,allowableCourantNumber;									//max courant number
	int nSimpleLoops,uVelocityLoops,vVelocityLoops,wVelocityLoops,nPressureLoops,nSimplerLoops;

	int i,j,k,l,t,courantNumberControl,thomasI,simpleVar,simplerVar,nonLinear;      //indices
	int ieast,iwest,jnorth,jsouth,ktop,kbottom,ieasteast,jnorthnorth,ktoptop;								
	
	double fxe,fxp,fyn,fyp,fzt,fzp,fxw,fys,fzb;

	double dxpe,dypn,dzpt,dxpw,dyps,dzpb;													//transoformation
	
	double lower[250],upper[250],c[250],x[250],m[250],cc[250],middle[250];       	//MATRIX A,X AND B

	double pUnderRelaxCoeff,vUnderRelaxCoeff,uUnderRelaxCoeff,wUnderRelaxCoeff,urecru,urecrv,urecrw,urecrp; // under relaxation coefficient for pressure

	double rho,mu; 																	//density and viscosity
	
	//area
	double s;

	FILE *fp; 																		//file pointer
	
	char str[100];

	int writeOutputSteps;


	//Deferred correction
	double DCvalue,fuuds,fvuds,fwuds,fucds,fvcds,fwcds,fuQ,fvQ,fwQ,fuH,fvH,fwH,fuQHP,fvQHP,fwQHP,fuP,fvP,fwP,f1,f2,f3,f4,f5,f6;

	//grid
	double xc[250],yc[250],zc[250],xf[250],yf[250],zf[250],fx[250],fy[250],fz[250];

	//accuracy of outer iterations and p and v loops
	double accuracy;
	
	//scheme
	double scheme;
	
	//counts no of outer iterations
	int check;



#endif

