#include"velocity.h"

//calculates the value of initial u velocity
void uEqn()
{
	
	for (i=2;i<=nxm;i++)
	{
		for (j=2;j<=nym;j++)
		{
			for(k=2;k<=nzm;k++)
			{
				ap[i][j][k]  = (- ae[i][j][k] - aw[i][j][k] - an[i][j][k] - as[i][j][k] - at[i][j][k] - ab[i][j][k] + apu[i][j][k])*urecru;
      			su[i][j][k]  = su[i][j][k] + (1 - uUnderRelaxCoeff) * ap[i][j][k] * u[i][j][k][l];
      			apu[i][j][k] = 1 / ap[i][j][k];	/* For implementation in the pressure correction approach */
				
			}
    	}
  }


for(nonLinear=0;nonLinear<=uVelocityLoops;nonLinear++)
{

// Creating matrices for ADI-TDMA algorithm  
  for (j=2;j<=nym;j++)
	{
		for (k=2;k<=nzm;k++)
		{
			for(i=2;i<=nxm;i++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				c[i] 	= (su[i][j][k] - (an[i][j][k]*u[i][jnorth][k][l] + as[i][j][k]*u[i][jsouth][k][l] + at[i][j][k]*u[i][j][ktop][l] + ab[i][j][k]*u[i][j][kbottom][l]))*apu[i][j][k];
				lower[i]= aw[i][j][k]*apu[i][j][k];
				upper[i]= ae[i][j][k]*apu[i][j][k];
				
				Thomas(nxm);
				
				//Equating solved value back to velocity
				for(thomasI=2;thomasI<=nxm;thomasI++)
	  				u[thomasI][j][k][l]=x[thomasI];				
			}
    	}
  }
  
  uMaxRes=0;
  for (j=2;j<=nym;j++)
	{
		for (k=2;k<=nzm;k++)
		{
			for(i=2;i<=nxm;i++)
			{
				ieast  = i + 1;
	  			iwest  = i - 1;
	  			jnorth = j + 1;
	  			jsouth = j - 1;
	  			ktop   = k + 1;
	  			kbottom= k - 1;
	  			
				uRes[i][j][k] 	= su[i][j][k] - (an[i][j][k]*u[i][jnorth][k][l] + as[i][j][k]*u[i][jsouth][k][l] + at[i][j][k]*u[i][j][ktop][l] + ab[i][j][k]*u[i][j][kbottom][l] + ae[i][j][k]*u[ieast][j][k][l] + aw[i][j][k]*u[iwest][j][k][l] + ap[i][j][k]*u[i][j][k][l]);
				uMaxRes	= uMaxRes + fabs(uRes[i][j][k]);
			}
    	}
   }
//   printf(" u Res = %e\n",uMaxRes);
   if(fabs(uMaxRes)>accuracy)		uVelocityLoops = uVelocityLoops + 1;
}
  
}

