#include"boundaryConditions.h"

void boundary_pressure(double ****p)
{
  int js = 1, jn = ny;
  int iw = 1, ie = nx;
  int kb = 1, kt = nz;
  
  /* Top and Bottom pressure boundary conditions */
  for (i=2;i<=nxm;i++)
  for (k=2;k<=nzm;k++)
  {
   p[i][js][k][l] = p[i][js+1][k][l];
   p[i][jn][k][l] = p[i][jn-1][k][l];
  }
  
  /* East adn West boundary conditions */
  for (j=2;j<=nym;j++)
  for (k=2;k<=nzm;k++)
  {
     p[iw][j][k][l] = p[iw+1][j][k][l];
     p[ie][j][k][l] = p[ie-1][j][k][l];
  }
  
  //Front and Back boundary Conditions
  for (j=2;j<=nym;j++)
  for (i=2;i<=nxm;i++)
  {
     p[i][j][kb][l] = p[i][j][kb+1][l];
     p[i][j][kt][l] = p[i][j][kt-1][l];
  }
  
}

