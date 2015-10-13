#include"pressure.h"

void pressure_correction_east()
{
	double d;	/* rho*dx */
  	double dpxel, uel, apue;	/* For the interpolation work */
  	double ue, dpxe;	/* Rhie chow interpolated cell face velocities */ 
  	int ie;
  
  
  	for (i=2;i<=nxm-1;i++)
	  {
		ieast = i+1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;		
		fxp = 1.0 - fxe;

	  	for (j=2;j<=nym;j++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	{
	      			/* Since, volume of the east cell equal to that of the central face, */
	      			/* interpolated cel face quantites can be written as follows */
	      
		      		s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
	      			vole = dxpe * s;
	      
	      			d = rho*s;
	      
	      			dpxel = 0.5*(dpx[i][j][k] + dpx[ieast][j][k]);
	      			uel   = fxp * u[i][j][k][l] +   fxe * u[ieast][j][k][l];
	      			apue  = fxp * apu[i][j][k] + fxe * apu[ieast][j][k];
	      
	      			/* Evaluating cell face gradient and velocity */
	      			dpxe  = (  p[ieast][j][k][l] -  p[i][j][k][l])/dxpe;
	      			ue    = uel - apue*vole*(dpxe - dpxel);
	      			Fe[i][j][k] = d*ue;	/* Correcting the fluxes with interpolated velocities */
	      
	      			/* Constructing the coefficients */
	      
	      			ae[i][j][k]  = -d*apue*s;
	      			aw[ieast][j][k]= ae[i][j][k];
              		}
	      	  }
	     }
}

void pressure_correction_north()
{
	double d;	/* rho*dx */
  	double dpynl, vnl, apvn;	/* For the interpolation work */
  	double vn, dpyn;	/* Rhie chow interpolated cell face velocities and gradients */ 
  	int jn;
  
  
  for (j=2;j<=nym-1;j++)
  {
	jnorth = j+1;
	dypn = yc[jnorth] - yc[j];
	fyn = (yf[j] - yc[j])/dypn;	
    	fyp = 1.0 - fyn;

  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
	  		
	  		s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
			voln = s * dypn;
			d = rho * s;
		
      			/* Since, volume of the east cell equal to that of the central face, */
      			/* interpolated cel face quantites can be written as follows */
      			dpynl = 0.5*(dpy[i][j][k] + dpy[i][jnorth][k]);
      			vnl   = fyp * v[i][j][k][l] +   fyn * v[i][jnorth][k][l];
      			apvn  = fyp * apv[i][j][k] + fyn * apv[i][jnorth][k];
      	
      			/* Evaluating cell face gradient and velocity */
      			dpyn  = ( p[i][jnorth][k][l] -  p[i][j][k][l])/dypn;
      			vn    = vnl - apvn*voln*(dpyn - dpynl);
      			Fn[i][j][k] = d*vn;	/* Correcting the fluxes with interpolated velocities */
      	
      			/* Constructing the coefficients */
      	
      			an[i][j][k]  = -d*apvn*s;
      			as[i][jnorth][k] = an[i][j][k];
     		}
    	}
  }
}

void pressure_correction_top()
{
	double d;	/* rho*dx */
  	double dpztl, wtl, apwt;	/* For the interpolation work */
  	double wt, dpzt;	/* Rhie chow interpolated cell face velocities and gradients */ 
  	int kt;
  
  
  for (k=2;k<=nzm-1;k++)
  {
	ktop = k+1;
	dzpt = zc[ktop]-zc[k];
	fzt = (zf[k] - zc[k])/dzpt;
	fzp = 1.0 - fzt; 
	
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
	  		
	  		s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
			volt = s * dzpt;
			d = rho * s;
		
	      		/* Since, volume of the east cell equal to that of the central face, */
	      		/* interpolated cel face quantites can be written as follows */
	      		dpztl = 0.5*(dpz[i][j][k] + dpz[i][j][ktop]);
	      		wtl   = fzp * w[i][j][k][l] +   fzt * w[i][j][ktop][l];
	      		apwt  = fzp * apw[i][j][k] + fzt * apw[i][j][ktop];
      
	      		/* Evaluating cell face gradient and velocity */
	      		dpzt  = ( p[i][j][ktop][l] -  p[i][j][k][l])/dzpt;
	      		wt    = wtl - apwt*volt*(dpzt - dpztl);
	      		Ft[i][j][k] = d*wt;	/* Correcting the fluxes with interpolated velocities */
      
	      		/* Constructing the coefficients */
	      
	      		at[i][j][k]  = -d*apwt*s;
	      		ab[i][j][ktop] = at[i][j][k];
     		}
    	}
  }
}

void pressure_correction_source()
{
	maxMassResidual=0; /* To check for continuity */
  
  	for (i=2;i<=nxm;i++)
	  {
	  	for (j=2;j<=nym;j++)
		  {
		  	for(k=2;k<=nzm;k++)
		  	 {
		  	 	ieast 	= i + 1;
      				iwest 	= i - 1;
      				jnorth	= j + 1;
      				jsouth	= j - 1;
      				ktop  	= k + 1;
      				kbottom	= k - 1;
            
      				su[i][j][k] = Fe[iwest][j][k] - Fe[i][j][k] + Fn[i][jsouth][k] - Fn[i][j][k] + Ft[i][j][kbottom] - Ft[i][j][k] ;
      				ap[i][j][k] = -(ae[i][j][k] + aw[i][j][k] + an[i][j][k] + as[i][j][k] + at[i][j][k] + ab[i][j][k]);
      				maxMassResidual	= maxMassResidual + fabs(su[i][j][k]);	/* Checking continuity */
      				pp[i][j][k][l] = 0;
    		 }
    	}
  	 }
}

