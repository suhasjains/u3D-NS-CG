#include"velocity.h"
#include"../pv_Coupling/pv_Coupling.h"

//Updates u Equation Coefficients using hybrid scheme
void uEqnCoeff()
{
	
	//double ue,uw,ve,vw,we,ww;
	
	//Constructing convective fluxes along east face 	
	for  (i=2; i<=nxm-1; i++)
 	{
		ieast = i+1;
		dxpe = xc[ieast] - xc[i];
		fxe = (xf[i]-xc[i])/dxpe;
		fxp = 1.0 - fxe;

	 	for (j=2; j<=nym; j++)
   		{
	  		for(k=2; k<=nzm; k++)
			{
		
		
				/* Evaluation of cell face area */
				s = (yf[j]-yf[j-1])*(zf[k]-zf[k-1]);
		
				/* Evaluation of diffusive term */	
				diff = mu*s/dxpe;
		
				ieasteast 	= i+2;
				iwest  		= i-1;
		
				upwind_implicit(Fe);

      				if(scheme==1)
				{
					/* Constructing co-efficients for upwind*/
			      		ae[i][j][k]     =  convf - diff;
			      		aw[ieast][j][k] = -convp - diff;
				}
				else if(scheme==2)
				{
					/* Constructing co-efficients for hybrid*/
			      		ae[i][j][k]     =  convf - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0);
			      		aw[ieast][j][k] = -convp - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0);
      				}
				else if(scheme==3)
				{
					/* Constructing co-efficients for power law*/
			      		ae[i][j][k]     =  convf - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0);
      					aw[ieast][j][k] = -convp - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0);
				}

//				Deferred Correction Upwind Flux
				fuuds = convp*u[i][j][k][l] + convf*u[ieast][j][k][l];
      				fvuds = convp*v[i][j][k][l] + convf*v[ieast][j][k][l];
      				fwuds = convp*w[i][j][k][l] + convf*w[ieast][j][k][l];
      	
//				Deferred Correction Central Differencing Flux      	
      				fucds = Fe[i][j][k]*(u[ieast][j][k][l]*fxe + u[i][j][k][l]*fxp);
      				fvcds = Fe[i][j][k]*(v[ieast][j][k][l]*fxe + v[i][j][k][l]*fxp);
      				fwcds = Fe[i][j][k]*(w[ieast][j][k][l]*fxe + w[i][j][k][l]*fxp);
      	
      		//		Deferred Correction Hybrid Flux
      				fuH = convp*u[i][j][k][l] + convf*u[ieast][j][k][l] - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(u[ieast][j][k][l] - u[i][j][k][l]);
      				fvH = convp*v[i][j][k][l] + convf*v[ieast][j][k][l]	- diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(v[ieast][j][k][l] - v[i][j][k][l]);
      				fwH = convp*w[i][j][k][l] + convf*w[ieast][j][k][l] - diff*std::max(1 - 0.5*fabs(Fe[i][j][k]/diff),0.0)*(w[ieast][j][k][l] - w[i][j][k][l]);
      	
		//		Deferred Correction Power Law Flux  
				fuP = convp*u[i][j][k][l] + convf*u[ieast][j][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0)*(u[ieast][j][k][l] - u[i][j][k][l]);
      				fvP = convp*v[i][j][k][l] + convf*v[ieast][j][k][l]	- diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0)*(v[ieast][j][k][l] - v[i][j][k][l]);
      				fwP = convp*w[i][j][k][l] + convf*w[ieast][j][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fe[i][j][k]/diff),5),0.0)*(w[ieast][j][k][l] - w[i][j][k][l]);  
      	
      	//			Deferred Correction QUICK Flux
				fuQ = DC_QUICK_flux_east(fuQ,u);
				fvQ = DC_QUICK_flux_east(fvQ,v);
				fwQ = DC_QUICK_flux_east(fwQ,w);

		

		//		Deferred Correction QUICK Flux when hybrid and Power Law is used to predict
				fuQHP = DC_QUICK_flux_east(fuQHP,u);
				fuQHP = fuQHP - diff*(u[ieast][j][k][l] - u[i][j][k][l]);			
		
				fvQHP = DC_QUICK_flux_east(fvQHP,v);
				fvQHP = fvQHP - diff*(v[ieast][j][k][l] - v[i][j][k][l]);
				
				fwQHP = DC_QUICK_flux_east(fwQHP,w);
				fwQHP = fwQHP - diff*(w[ieast][j][k][l] - w[i][j][k][l]);     
      	
      
      				//Adding contribution to source term
      				if(scheme==1)	v_east_source_dc(fuuds,fvuds,fwuds,fuQ,fvQ,fwQ);
      				if(scheme==2)	v_east_source_dc(fuH,fvH,fwH,fuQHP,fvQHP,fwQHP);
				if(scheme==3)	v_east_source_dc(fuP,fvP,fwP,fuQHP,fvQHP,fwQHP);
			}
  		}
	}
	
	//Constructing convective fluxes along north face
	for(j=2;j<=nym-1;j++)
	{
		jnorth = j+1;
		dypn = yc[jnorth] - yc[j];
		fyn = (yf[j] - yc[j])/dypn;
	    	fyp = 1.0 - fyn;

		for(i=2;i<=nxm;i++)
		{
			for(k=2;k<=nzm;k++)
			{
				
				s = (xf[i]-xf[i-1])*(zf[k]-zf[k-1]);
      
	      			jnorthnorth	= j+2;
	      			jsouth		= j-1;
      
		      		diff = mu*s/dypn;
      	
      				upwind_implicit(Fn);
      			
				if(scheme==1)
				{
					/* Constructing co-efficients for upwind */
      					an[i][j][k]	   	=  convf - diff;
      					as[i][jnorth][k]= -convp - diff;
				}
				else if(scheme==2)
				{
					/* Constructing co-efficients for hybrid */
      					an[i][j][k]	   	=  convf - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0);
      					as[i][jnorth][k]= -convp - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0);
      				}
				else if(scheme==3)
				{
					/* Constructing co-efficients for power law */
      					an[i][j][k]	   	=  convf - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0);
      					as[i][jnorth][k]= -convp - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0);
				}
				
				//		Deferred Correction Upwind Flux
				fuuds = convp*u[i][j][k][l] + convf*u[i][jnorth][k][l];
      				fvuds = convp*v[i][j][k][l] + convf*v[i][jnorth][k][l];
      				fwuds = convp*w[i][j][k][l] + convf*w[i][jnorth][k][l];
      			
      			//		Deferred Correction Central Differencing Flux
				fucds = Fn[i][j][k]*(u[i][jnorth][k][l]*fyn + u[i][j][k][l]*fyp);
      				fvcds = Fn[i][j][k]*(v[i][jnorth][k][l]*fyn + v[i][j][k][l]*fyp);
      				fwcds = Fn[i][j][k]*(w[i][jnorth][k][l]*fyn + w[i][j][k][l]*fyp);
      			
      			//		Deferred Correction hybrid Flux
      				fuH = convp*u[i][j][k][l] + convf*u[i][jnorth][k][l] - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0)*(u[i][jnorth][k][l]-u[i][j][k][l]);
      				fvH = convp*v[i][j][k][l] + convf*v[i][jnorth][k][l] - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0)*(v[i][jnorth][k][l]-v[i][j][k][l]);
      				fwH = convp*w[i][j][k][l] + convf*w[i][jnorth][k][l] - diff*std::max(1 - 0.5*fabs(Fn[i][j][k]/diff),0.0)*(w[i][jnorth][k][l]-w[i][j][k][l]);
      			
      			//		Deferred Correction Power Law Flux
      				fuP = convp*u[i][j][k][l] + convf*u[i][jnorth][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0)*(u[i][jnorth][k][l]-u[i][j][k][l]);
      				fvP = convp*v[i][j][k][l] + convf*v[i][jnorth][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0)*(v[i][jnorth][k][l]-v[i][j][k][l]);
      				fwP = convp*w[i][j][k][l] + convf*w[i][jnorth][k][l] - diff*std::max(pow(1 - 0.1*fabs(Fn[i][j][k]/diff),5),0.0)*(w[i][jnorth][k][l]-w[i][j][k][l]);

				//		Deferred Correction QUICK Flux
				fuQ = DC_QUICK_flux_north(fuQ,u);
				fvQ = DC_QUICK_flux_north(fvQ,v);
				fwQ = DC_QUICK_flux_north(fwQ,w);

				//		Deferred Correction QUICK Flux when hybrid and power Law is used to predict
				fuQHP = DC_QUICK_flux_north(fuQHP,u);				
				fuQHP = fuQHP - diff*(u[i][jnorth][k][l]-u[i][j][k][l]);
				
				fvQHP = DC_QUICK_flux_north(fvQHP,v);		
				fvQHP = fvQHP - diff*(v[i][jnorth][k][l]-v[i][j][k][l]);
				  
				fwQHP = DC_QUICK_flux_north(fwQHP,w);
				fwQHP = fwQHP - diff*(w[i][jnorth][k][l]-w[i][j][k][l]);    
      			
       			
      				if(scheme==1)	v_north_source_dc(fuuds,fvuds,fwuds,fuQ,fvQ,fwQ);
      				if(scheme==2)	v_north_source_dc(fuH,fvH,fwH,fuQHP,fvQHP,fwQHP);
				if(scheme==3)	v_north_source_dc(fuP,fvP,fwP,fuQHP,fvQHP,fwQHP); 
      
			}
		}
	}

	//Constructing convective fluxes along top face
	for(k=2;k<=nzm-1;k++)
	{
		ktop = k+1;
		dzpt = zc[ktop]-zc[k];
		fzt = (zf[k] - zc[k])/dzpt;
		fzp = 1.0 - fzt; 

		for(i=2;i<=nxm;i++)
		{
			for(j=2;j<=nym;j++)
			{
				s = (xf[i]-xf[i-1])*(yc[j]-yf[j-1]);
      
	      			ktoptop	= k+2;
	      			kbottom	= k-1;
      
			      	diff = mu*s/dzpt;
      
      				upwind_implicit(Ft);
      			
				if(scheme==1)
				{
					/* Constructing co-efficients for upwind */
      					at[i][j][k]	   	=  convf - diff;
      					ab[i][j][ktop]	= -convp - diff;
				}
				else if(scheme==2)
				{
					/* Constructing co-efficients for hybrid */
      				at[i][j][k]	   	=  convf - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0);
      				ab[i][j][ktop]	= -convp - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0);
      				}
				else if(scheme==3)
				{
					/* Constructing co-efficients for power law */
      					at[i][j][k]	   	=  convf - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0);
      					ab[i][j][ktop]	= -convp - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0);
				}

				
				
				//Deferred Correction Upwind Flux
				fuuds = convp*u[i][j][k][l] + convf*u[i][j][ktop][l];
      				fvuds = convp*v[i][j][k][l] + convf*v[i][j][ktop][l];
      				fwuds = convp*w[i][j][k][l] + convf*w[i][j][ktop][l];
      			
      			//		Deferred Correction Central Differencing Flux
				fucds = Ft[i][j][k]*(u[i][j][ktop][l]*fzt + u[i][j][k][l]*fzp);
      				fvcds = Ft[i][j][k]*(v[i][j][ktop][l]*fzt + v[i][j][k][l]*fzp);
      				fwcds = Ft[i][j][k]*(w[i][j][ktop][l]*fzt + w[i][j][k][l]*fzp);
      			
      			//		Deferred Correction Hybrid Flux
      				fuH = convp*u[i][j][k][l] + convf*u[i][j][ktop][l] - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0)*(u[i][j][ktop][l] - u[i][j][k][l]);
      				fvH = convp*v[i][j][k][l] + convf*v[i][j][ktop][l] - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0)*(v[i][j][ktop][l] - v[i][j][k][l]);
      				fwH = convp*w[i][j][k][l] + convf*w[i][j][ktop][l] - diff*std::max(1 - 0.5*fabs(Ft[i][j][k]/diff),0.0)*(w[i][j][ktop][l] - w[i][j][k][l]);
      			
      			//		Deferred Correction Power Law Flux
      				fuP = convp*u[i][j][k][l] + convf*u[i][j][ktop][l] - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0)*(u[i][j][ktop][l] - u[i][j][k][l]);
      				fvP = convp*v[i][j][k][l] + convf*v[i][j][ktop][l] - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0)*(v[i][j][ktop][l] - v[i][j][k][l]);
      				fwP = convp*w[i][j][k][l] + convf*w[i][j][ktop][l] - diff*std::max(pow(1 - 0.1*fabs(Ft[i][j][k]/diff),5),0.0)*(w[i][j][ktop][l] - w[i][j][k][l]);
				
				//		Deferred Correction QUICK Flux
				fuQ = DC_QUICK_flux_top(fuQ,u);
				fvQ = DC_QUICK_flux_top(fvQ,v);
				fwQ = DC_QUICK_flux_top(fwQ,w);
				
				//		Deferred Correction QUICK Flux when hybrid and Power Law is used to predict
				fuQHP = DC_QUICK_flux_top(fuQHP,u);	
				fuQHP = fuQHP - diff*(u[i][j][ktop][l] - u[i][j][k][l]);
				
				fvQHP = DC_QUICK_flux_top(fvQHP,v);					
				fvQHP = fvQHP - diff*(v[i][j][ktop][l] - v[i][j][k][l]);
				  
				fwQHP = DC_QUICK_flux_top(fwQHP,w);					
				fwQHP = fwQHP - diff*(w[i][j][ktop][l] - w[i][j][k][l]);    
      			
      			
       
      				if(scheme==1)	v_top_source_dc(fuuds,fvuds,fwuds,fuQ,fvQ,fwQ);
      				if(scheme==2)	v_top_source_dc(fuH,fvH,fwH,fuQHP,fvQHP,fwQHP);
				if(scheme==3)	v_top_source_dc(fuP,fvP,fwP,fuQHP,fvQHP,fwQHP);
      							
				
			}
		}
	}
	
	
	//Constructing Source Terms
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
				
				vol	= dx*dy*dz;
				
				
				/* Constructing the pressure gradient terms */
      				pe = fxe	*p[ieast][j][k][l]  	+ fxp*p[i][j][k][l];
      				pw = (1.0-fxw)	*p[iwest][j][k][l]	+ fxw*p[i][j][k][l];
      				pn = fyn	*p[i][jnorth][k][l] 	+ fyp*p[i][j][k][l];
      				ps = (1.0-fys)	*p[i][jsouth][k][l] 	+ fys*p[i][j][k][l];
      				pt = fzt	*p[i][j][ktop][l]	+ fzp*p[i][j][k][l];
      				pb = (1.0-fzb)	*p[i][j][kbottom][l]	+ fzb*p[i][j][k][l];
      
      				dpx[i][j][k] = (pe - pw)/dx;
      				dpy[i][j][k] = (pn - ps)/dy;
      				dpz[i][j][k] = (pt - pb)/dz;
      			
      			
      				/* Updating source terms with pressure gradient */
      				su[i][j][k] = su[i][j][k] - dpx[i][j][k]*vol;
      				sv[i][j][k] = sv[i][j][k] - dpy[i][j][k]*vol;
				sw[i][j][k] = sw[i][j][k] - dpz[i][j][k]*vol;
								
				
			}
		}
	}
	
	
//	Adding unsteady contribution to source terms
	for(i=2;i<=nxm;i++)
	{
		for(j=2;j<=nym;j++)
		{
			for(k=2;k<=nzm;k++)
			{
				
				apt	= rho*vol/dt;
				su[i][j][k]	= su[i][j][k] + apt*u[i][j][k][l-1];
				sv[i][j][k]	= sv[i][j][k] + apt*v[i][j][k][l-1]; 
				sw[i][j][k]	= sw[i][j][k] + apt*w[i][j][k][l-1];
				apu[i][j][k]= apu[i][j][k] + apt;  
				apv[i][j][k]= apv[i][j][k] + apt;
				apw[i][j][k]= apw[i][j][k] + apt;
			}
		}
	}
}

