#include"boundaryConditions.h"

//Velocity boundary Conditions
void velocityBoundaryConditions()	//Check using 3D code
{
	/* Bottom boundary conditions */
  	j = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
			//diff = 2*mu*dx*dz/dy;	  		
			diff = mu*(xf[i]-xf[i-1])*(zf[k]-zf[k-1])/(yc[j]-yc[j-1]);
    		apu[i][j][k] = apu[i][j][k] + diff;
    		apw[i][j][k] = apw[i][j][k] + diff;
	    	su[i][j][k]  = su[i][j][k]  + diff*u[i][j-1][k][l];
	    	sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j-1][k][l];

		

	  	}
  	  }
  	
	    
  	  /* Top boundary conditions */
  	j = nym;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
		//diff = 2*mu*dx*dz/dy;	  		
		diff = mu*(xf[i]-xf[i-1])*(zf[k]-zf[k-1])/(yc[j+1]-yc[j]);
    		apu[i][j][k] = apu[i][j][k] + diff;
    		apw[i][j][k] = apw[i][j][k] + diff;
	    	su[i][j][k]  = su[i][j][k]  + diff*u[i][j+1][k][l];
	    	sw[i][j][k]  = sw[i][j][k]  + diff*w[i][j+1][k][l];
 		
	  	}
	  }
  
  	/* West boundary conditions */
	i = 2;
  	for (j=2;j<=nym;j++)
	  {
	  	for(k=2;k<=nzm;k++)
	  	{
		//diff = 2*mu*dy*dz/dx;	  	
		diff = mu*(yf[j]-yf[j-1])*(zf[k]-zf[k-1])/(xc[i]-xc[i-1]);
    		apv[i][j][k] = apv[i][j][k] + diff;
    		apw[i][j][k] = apw[i][j][k] + diff;
    		sv[i][j][k]  = sv[i][j][k]  + diff*v[i-1][j][k][l];
	  		sw[i][j][k]  = sw[i][j][k]  + diff*w[i-1][j][k][l];

		//printf("%e	%e\n",yf[j]-yf[j-1],dy);
		}
    }
          
  /* East boundary conditions */
  i = nxm;
  for (j=2;j<=nym;j++)
  {
  	for(k=2;k<=nzm;k++)
  	{
	//diff = 2*mu*dy*dz/dx;  	
	diff = mu*(yf[j]-yf[j-1])*(zf[k]-zf[k-1])/(xc[i+1]-xc[i]);
    	apv[i][j][k] = apv[i][j][k] + diff;
    	apw[i][j][k] = apw[i][j][k] + diff;
    	sv[i][j][k]  = sv[i][j][k]  + diff*v[i+1][j][k][l];
	  	sw[i][j][k]  = sw[i][j][k]  + diff*w[i+1][j][k][l];
  	}
  }  
  
  /* Back bounday conditions */
	k = 2;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
		//diff = 2*mu*dy*dx/dz;	  	
		diff = mu*(yf[j]-yf[j-1])*(xf[i]-xf[i-1])/(zc[k]-zc[k-1]);
//    		apv[i][j][k] = apv[i][j][k] + diff;
//    		apu[i][j][k] = apu[i][j][k] + diff;
//    		sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k-1][l];
//	  		su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k-1][l];
			apw[i][j][k] = apw[i][j][k] + diff;
			u[i][j][k-1][l] = u[i][j][k][l];
			v[i][j][k-1][l] = v[i][j][k][l]; 
				

		}
    }
          
  /* Front boundary conditions */
  k = nzm;
  	for (i=2;i<=nxm;i++)
	  {
	  	for(j=2;j<=nym;j++)
	  	{
		//diff = 2*mu*dy*dx/dz;	  	
		diff = mu*(yf[j]-yf[j-1])*(xf[i]-xf[i-1])/(zc[k+1]-zc[k]);
//    		apv[i][j][k] = apv[i][j][k] + diff;
//    		apu[i][j][k] = apu[i][j][k] + diff;
//    		sv[i][j][k]  = sv[i][j][k]  + diff*v[i][j][k+1][l];
//	  		su[i][j][k]  = su[i][j][k]  + diff*u[i][j][k+1][l];
			apw[i][j][k] = apw[i][j][k] + diff;
			u[i][j][k+1][l] = u[i][j][k][l];
			v[i][j][k+1][l] = v[i][j][k][l]; 


		}
    }
}

void setting_ulid()
{
  j=ny;					/* Top boundary */
  for (i=2;i<=nxm;i++)
   for(k=2;k<=nzm;k++)
   {
    u[i][j][k][1] = uTopWallVel;
   }
}

