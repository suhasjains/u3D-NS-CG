#include"update.h"

void correctFaceFlux()
{
	/* Correcting the fluxes */
  ppo = pp[2][2][2][1];
  
    for (i=2;i<=nxm-1;i++)
	{
		for (j=2;j<=nym;j++)
		{
			for(k=2;k<=nzm;k++)
	  		{
	  			Fe[i][j][k] = Fe[i][j][k] + ae[i][j][k]*(pp[i+1][j][k][l] - pp[i][j][k][l]);	/* East face mass flux correction */
	  		}
        }
    }
    
    for (j=2;j<=nym-1;j++)
	{
		for (i=2;i<=nxm;i++)
		{
			for(k=2;k<=nzm;k++)
			{
				Fn[i][j][k] = Fn[i][j][k] + an[i][j][k]*(pp[i][j+1][k][l] - pp[i][j][k][l]);	/* North face mass flux correction */
			}
      	}
    }
    
    for (k=2;k<=nzm-1;k++)
	{
		for (i=2;i<=nxm;i++)
		{
			for(j=2;j<=nym;j++)
			{
				Ft[i][j][k] = Ft[i][j][k] + at[i][j][k]*(pp[i][j][k+1][l] - pp[i][j][k][l]);	/* Top face mass flux correction */
			}
      	}
    }
}

