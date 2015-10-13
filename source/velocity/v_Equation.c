#include"velocity.h"

//calculates the value of initial v velocity
void vEqn()
{
	 for (i=2;i<=nxm;i++)
	 {
	 	for (j=2;j<=nym;j++)
		 {
		 	for(k=2;k<=nzm;k++)
		 	{
		 		ap[i][j][k]  = (- ae[i][j][k] - aw[i][j][k] - an[i][j][k] - as[i][j][k] - at[i][j][k] - ab[i][j][k] + apv[i][j][k])*urecrv;
      			sv[i][j][k]  = sv[i][j][k] + (1 - vUnderRelaxCoeff) * ap[i][j][k] * v[i][j][k][l] ;
      			apv[i][j][k] = 1 / ap[i][j][k];	/* For implementation in the pressure correction approach */
		 		
		 	}
 	     }
  	 }
  	 
  for(nonLinear=0;nonLinear<=vVelocityLoops;nonLinear++)
  {
  	
//  Creating matrices for ADI-TDMA algorithm
	for (k=2;k<=nzm;k++)
	{
		for (i=2;i<=nxm;i++)
		{
			for(j=2;j<=nym;j++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[j] 	= (sv[i][j][k] - (ae[i][j][k]*v[ieast][j][k][l] + aw[i][j][k]*v[iwest][j][k][l] + at[i][j][k]*v[i][j][ktop][l] + ab[i][j][k]*v[i][j][kbottom][l]))*apv[i][j][k];
				lower[j]= as[i][j][k]*apv[i][j][k];
				upper[j]= an[i][j][k]*apv[i][j][k];
				
				Thomas(nym);
				
				//Equating solved value back to velocity
				for(thomasI=2;thomasI<=nym;thomasI++)
	  				v[i][thomasI][k][l]=x[thomasI];				
			}
    	}
   }
  
  
	vMaxRes = 0;  
  	for (k=2;k<=nzm;k++)
	{
		for (i=2;i<=nxm;i++)
		{
			for(j=2;j<=nym;j++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
	  			vRes[i][j][k] 	= sv[i][j][k] - (an[i][j][k]*v[i][jnorth][k][l] + as[i][j][k]*v[i][jsouth][k][l] + at[i][j][k]*v[i][j][ktop][l] + ab[i][j][k]*v[i][j][kbottom][l] + ae[i][j][k]*v[ieast][j][k][l] + aw[i][j][k]*v[iwest][j][k][l] + ap[i][j][k]*v[i][j][k][l]);
				vMaxRes	= vMaxRes + fabs(vRes[i][j][k]);
	  			
			}
    	}
   }
   
//   printf(" v Res = %e\n",vMaxRes);
   if(fabs(vMaxRes)>accuracy)		vVelocityLoops = vVelocityLoops + 1;
   
 }
   
}
