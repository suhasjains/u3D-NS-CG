#include"pressure.h"
#include"../solver/solver.h"


//calculates the value of pressure correction
void pCorrEqn()
{
	
	for(nonLinear=0;nonLinear<=nPressureLoops;nonLinear++)
  {
	
	
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
	  			
				c[i] 	= (su[i][j][k] - (an[i][j][k]*pp[i][jnorth][k][l] + as[i][j][k]*pp[i][jsouth][k][l] + at[i][j][k]*pp[i][j][ktop][l] + ab[i][j][k]*pp[i][j][kbottom][l]))/ap[i][j][k];
				lower[i]= aw[i][j][k]/ap[i][j][k];
				upper[i]= ae[i][j][k]/ap[i][j][k];
				
				Thomas(nxm);
				
				//Equating solved value back to velocity
				for(thomasI=2;thomasI<=nxm;thomasI++)
	  				pp[thomasI][j][k][l]=x[thomasI];				
			}
    	}
  }
  
  pMaxRes=0;
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
	  			
	  			pRes[i][j][k] 	= su[i][j][k] - (an[i][j][k]*pp[i][jnorth][k][l] + as[i][j][k]*pp[i][jsouth][k][l] + at[i][j][k]*pp[i][j][ktop][l] + ab[i][j][k]*pp[i][j][kbottom][l] + ae[i][j][k]*pp[ieast][j][k][l] + aw[i][j][k]*pp[iwest][j][k][l] + ap[i][j][k]*pp[i][j][k][l]);
				pMaxRes	= pMaxRes + fabs(pRes[i][j][k]);
	  			
			}
    	}
    }
    
//    printf(" p Res = %e\n",pMaxRes);
    if(fabs(pMaxRes)>accuracy)		nPressureLoops = nPressureLoops + 1;
    
  }
}

