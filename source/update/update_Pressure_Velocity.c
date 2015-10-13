#include"update.h"

//Updates the values of the pressure
void update_Pressure_Velocity()
{
	double ppe,ppw,ppn,pps,ppt,ppb; /* Linear interpolation of pprime */
  
  for(i=2;i<=nxm;i++)
  {
	dx = xf[i] - xf[i-1];
	ieast 	= i + 1;
      	iwest 	= i - 1;
	dxpe = xc[ieast] - xc[i];
	fxe = (xf[i]-xc[i])/dxpe;		
	fxp = 1.0 - fxe;
	dxpw = xc[i] - xc[iwest];
	fxw = (xf[iwest]-xc[iwest])/dxpw;

    	for(j=2;j<=nym;j++)
	  {
		dy = yf[j] - yf[j-1];
		jnorth	= j + 1;
      		jsouth	= j - 1;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;	
    		fyp = 1.0 - fyn;
		dyps = yc[j] - yc[jsouth];
		fys = (yf[jsouth] - yc[jsouth])/dyps;
	
	  	for(k=2;k<=nzm;k++)
		  {
			dz = zf[k] - zf[k-1];
			ktop  	= k + 1;
		  	kbottom	= k - 1;
			dzpt = zc[ktop]-zc[k];
			fzt = (zf[k] - zc[k])/dzpt;
			fzp = 1.0 - fzt;
			dzpb = zc[k]-zc[kbottom];
			fzb = (zf[kbottom] - zc[kbottom])/dzpb;
			
		  	ppe = fxe*pp[ieast][j][k][l]  		+ fxp*pp[i][j][k][l];	/* Interpolating pprime */
      			ppw = (1.0 - fxw)*pp[iwest][j][k][l]	+ fxw*pp[i][j][k][l];
      			ppn = fyn*pp[i][jnorth][k][l] 		+ fyp*pp[i][j][k][l];
      			pps = (1.0 - fys)*pp[i][jsouth][k][l] 	+ fys*pp[i][j][k][l];
      			ppt = fzt*pp[i][j][ktop][l] 		+ fzp*pp[i][j][k][l];
      			ppb = (1.0 - fzb)*pp[i][j][kbottom][l]	+ fzb*pp[i][j][k][l];		  	
		
			
      			u[i][j][k][l] = u[i][j][k][l] - (ppe - ppw)*dy*dz*apu[i][j][k];	/* Correcting velocities and pressure */
      			v[i][j][k][l] = v[i][j][k][l] - (ppn - pps)*dx*dz*apv[i][j][k];
      			w[i][j][k][l] = w[i][j][k][l] - (ppt - ppb)*dx*dy*apw[i][j][k];
      		
			p[i][j][k][l] = p[i][j][k][l] + pUnderRelaxCoeff*(pp[i][j][k][l] - ppo);
      
		  } 
      }
  }
}


