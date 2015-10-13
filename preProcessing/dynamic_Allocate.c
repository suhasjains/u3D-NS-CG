#include"preProcessing.h"


//dynamically allocates the size of array
void dynamicAllocate()
{

   	u	 = new double***[xNumberOfCells+5];
	v	 = new double***[xNumberOfCells+5];
	w	 = new double***[xNumberOfCells+5];
	p	 = new double***[xNumberOfCells+5];
	pp	 = new double***[xNumberOfCells+5];
	xCoord	 = new double**[xNumberOfCells+5];
	yCoord	 = new double**[xNumberOfCells+5];
	zCoord	 = new double**[xNumberOfCells+5];
	Ft	 = new double**[xNumberOfCells+5];
	Fn	 = new double**[xNumberOfCells+5];
	Fe	 = new double**[xNumberOfCells+5];
	su	 = new double**[xNumberOfCells+5];
	sv	 = new double**[xNumberOfCells+5];
	sw	 = new double**[xNumberOfCells+5];
	ae	 = new double**[xNumberOfCells+5];
	aw	 = new double**[xNumberOfCells+5];
	an	 = new double**[xNumberOfCells+5];
	as	 = new double**[xNumberOfCells+5];
	at	 = new double**[xNumberOfCells+5];
	ab	 = new double**[xNumberOfCells+5];
	ap	 = new double**[xNumberOfCells+5];
	apu	 = new double**[xNumberOfCells+5];
	apv	 = new double**[xNumberOfCells+5];
	apw	 = new double**[xNumberOfCells+5];
	flux = new double**[xNumberOfCells+5];
	dpx = new double**[xNumberOfCells+5];
	dpy = new double**[xNumberOfCells+5];
	dpz = new double**[xNumberOfCells+5];
	uRes = new double**[xNumberOfCells+5];
	vRes = new double**[xNumberOfCells+5];
	wRes = new double**[xNumberOfCells+5];
	pRes = new double**[xNumberOfCells+5];
//	Sc = new double**[xNumberOfCells+5];


	for (int i=0; i<(xNumberOfCells+5); ++i)
	{
	u[i]	 = new double**[yNumberOfCells+5];
	v[i]	 = new double**[yNumberOfCells+5];
	w[i]	 = new double**[yNumberOfCells+5];
	p[i]	 = new double**[yNumberOfCells+5];
	pp[i]	 = new double**[yNumberOfCells+5];
	xCoord[i]	 = new double*[yNumberOfCells+5];
	yCoord[i]	 = new double*[yNumberOfCells+5];
	zCoord[i]	 = new double*[yNumberOfCells+5];
	Ft[i]	 = new double*[yNumberOfCells+5];
	Fn[i]	 = new double*[yNumberOfCells+5];
	Fe[i]	 = new double*[yNumberOfCells+5];
	su[i]	 = new double*[yNumberOfCells+5];
	sv[i]	 = new double*[yNumberOfCells+5];
	sw[i]	 = new double*[yNumberOfCells+5];
	ae[i]	 = new double*[yNumberOfCells+5];
	aw[i]	 = new double*[yNumberOfCells+5];
	an[i]	 = new double*[yNumberOfCells+5];
	as[i]	 = new double*[yNumberOfCells+5];
	at[i]	 = new double*[yNumberOfCells+5];
	ab[i]	 = new double*[yNumberOfCells+5];
	ap[i]	 = new double*[yNumberOfCells+5];
	apu[i]	 = new double*[yNumberOfCells+5];
	apv[i]	 = new double*[yNumberOfCells+5];
	apw[i]	 = new double*[yNumberOfCells+5];
	flux[i] = new double*[yNumberOfCells+5];
	dpx[i] = new double*[yNumberOfCells+5];
	dpy[i] = new double*[yNumberOfCells+5];
	dpz[i] = new double*[yNumberOfCells+5];
	uRes[i] = new double*[yNumberOfCells+5];
	vRes[i] = new double*[yNumberOfCells+5];
	wRes[i] = new double*[yNumberOfCells+5];
	pRes[i] = new double*[yNumberOfCells+5];
//	Sc[i] = new double*[yNumberOfCells+5];


  	  for (int j=0; j<(yNumberOfCells+5); ++j)
	  {
		u[i][j]	 = new double*[zNumberOfCells+5];
		v[i][j]	 = new double*[zNumberOfCells+5];
		w[i][j]	 = new double*[zNumberOfCells+5];
		p[i][j]	 = new double*[zNumberOfCells+5];
		pp[i][j]	 = new double*[zNumberOfCells+5];
		xCoord[i][j]	 = new double[zNumberOfCells+5];
		yCoord[i][j]	 = new double[zNumberOfCells+5];
		zCoord[i][j]	 = new double[zNumberOfCells+5];
		Ft[i][j]	 = new double[zNumberOfCells+5];
		Fn[i][j]	 = new double[zNumberOfCells+5];
		Fe[i][j]	 = new double[zNumberOfCells+5];
		su[i][j]	 = new double[zNumberOfCells+5];
		sv[i][j]	 = new double[zNumberOfCells+5];
		sw[i][j]	 = new double[zNumberOfCells+5];
		ae[i][j]	 = new double[zNumberOfCells+5];
		aw[i][j]	 = new double[zNumberOfCells+5];
		an[i][j]	 = new double[zNumberOfCells+5];
		as[i][j]	 = new double[zNumberOfCells+5];
		at[i][j]	 = new double[zNumberOfCells+5];
		ab[i][j]	 = new double[zNumberOfCells+5];
		ap[i][j]	 = new double[zNumberOfCells+5];
		apu[i][j]	 = new double[zNumberOfCells+5];
		apv[i][j]	 = new double[zNumberOfCells+5];
		apw[i][j]	 = new double[zNumberOfCells+5];
		flux[i][j] = new double[zNumberOfCells+5];
		dpx[i][j] = new double[zNumberOfCells+5];
		dpy[i][j] = new double[zNumberOfCells+5];
		dpz[i][j] = new double[zNumberOfCells+5];
		uRes[i][j] = new double[zNumberOfCells+5];
		vRes[i][j] = new double[zNumberOfCells+5];
		wRes[i][j] = new double[zNumberOfCells+5];
		pRes[i][j] = new double[zNumberOfCells+5];
//		Sc[i][j] = new double[zNumberOfCells+5];

		 for (int k=0; k<(zNumberOfCells+5); ++k)
		 {
		  u[i][j][k] 		 = new double[2];
		  v[i][j][k]		 = new double[2];
	 	  w[i][j][k]		 = new double[2];
	  	  p[i][j][k]		 = new double[2];
	  	  pp[i][j][k]	 = new double[2];
	 	 }

	  }
	}

}
