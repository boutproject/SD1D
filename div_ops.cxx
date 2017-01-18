/*
   Copyright B.Dudson, University of York, 2016
              email: benjamin.dudson@york.ac.uk
              
    This file is part of SD1D.

    SD1D is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SD1D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SD1D.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <mpi.h>

#include "div_ops.hxx"

#include <bout/mesh.hxx>
#include <globals.hxx>
#include <derivs.hxx>
#include <output.hxx>
#include <utils.hxx>
#include <bout/assert.hxx>

#include <cmath>


const Field3D Div_par_diffusion(const Field3D &K, const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal c = 0.5*(K(i,j,k) + K(i,j+1,k)); // K at the upper boundary
        BoutReal J = 0.5*(mesh->J(i,j) + mesh->J(i,j+1)); // Jacobian at boundary
        
        //if((i == mesh->xstart) && (j == mesh->ystart) && (k==0))
        //  output << "c = " << c << endl;
        
        BoutReal g_22 = 0.5*(mesh->g_22(i,j) + mesh->g_22(i,j+1));
        
        BoutReal gradient = 2.*(f(i,j+1,k) - f(i,j,k)) / (mesh->dy(i,j) + mesh->dy(i,j+1));
        
        BoutReal flux = c * J * gradient / g_22;
        
        result(i,j,k) += flux / (mesh->dy(i,j) * mesh->J(i,j));
        result(i,j+1,k) -= flux / (mesh->dy(i,j+1) * mesh->J(i,j+1));
      }
  return result;
}

const Field3D Div_par_spitzer(BoutReal K0, const Field3D &Te, bool bndry_flux) {
  Field3D result;
  result = 0.0;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface

        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;

          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }

	BoutReal Te0 = 0.5*(Te(i,j,k) + Te(i,j+1,k)); // Te at the upper boundary
        BoutReal K = K0*pow(Te0,2.5);
        BoutReal J = 0.5*(mesh->J(i,j) + mesh->J(i,j+1)); // Jacobian at boundary

        BoutReal g_22 = 0.5*(mesh->g_22(i,j) + mesh->g_22(i,j+1));

        BoutReal gradient = 2.*(Te(i,j+1,k) - Te(i,j,k)) / (mesh->dy(i,j) + mesh->dy(i,j+1));

        BoutReal flux = K * J * gradient / g_22;

        result(i,j,k) += flux / (mesh->dy(i,j) * mesh->J(i,j));
        result(i,j+1,k) -= flux / (mesh->dy(i,j+1) * mesh->J(i,j+1));
      }
  return result;
}

const Field3D Div_par_diffusion_upwind(const Field3D &K, const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal J = 0.5*(mesh->J(i,j) + mesh->J(i,j+1)); // Jacobian at boundary
        
        BoutReal g_22 = 0.5*(mesh->g_22(i,j) + mesh->g_22(i,j+1));
        
        BoutReal gradient = 2.*(f(i,j+1,k) - f(i,j,k)) / (mesh->dy(i,j) + mesh->dy(i,j+1));
        
        BoutReal c; // K at the upper boundary
        if(gradient > 0.0) {
          c = K(i,j+1,k);
        }else {
          c = K(i,j,k);
        }
        
        BoutReal flux = c * J * gradient / g_22;
        
        result(i,j,k) += flux / (mesh->dy(i,j) * mesh->J(i,j));
        result(i,j+1,k) -= flux / (mesh->dy(i,j+1) * mesh->J(i,j+1));
      }
  return result;
}

const Field3D Div_par_diffusion_index(const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal J = 0.5*(mesh->J(i,j) + mesh->J(i,j+1)); // Jacobian at boundary
        
        BoutReal gradient = f(i,j+1,k) - f(i,j,k);
        
        BoutReal flux = J * gradient;
        
        result(i,j,k) += flux / mesh->J(i,j);
        result(i,j+1,k) -= flux / mesh->J(i,j+1);
      }
  return result;
}

const Field3D AddedDissipation(const Field3D &N, const Field3D &P, const Field3D f, bool bndry_flux) {
  Field3D result = 0.0;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j >= mesh->yend-1) && mesh->lastY(i))
            continue;
          
          if((j <= mesh->ystart) && mesh->firstY(i))
            continue;
        }
        
        // At upper boundary
        
        BoutReal d = 0.5*(1./N(i,j,k) + 1./N(i,j+1,k));
        
        // Velocity 
        BoutReal v = - 0.25*d*( (P(i,j-1,k) + P(i,j+1,k) - 2.*P(i,j,k)) - (P(i,j,k) + P(i,j+2,k)-2.*P(i,j+1,k)) );
        
        // Variable being advected. Could use different interpolation?
        BoutReal var = 0.5*(f(i,j,k) + f(i,j+1,k));

        BoutReal flux = var * v * (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j))+ sqrt(mesh->g_22(i,j+1)));
        
        result(i,j,k) -= flux / (mesh->dy(i,j) * mesh->J(i,j));
        result(i,j+1,k) += flux / (mesh->dy(i,j+1) * mesh->J(i,j+1));
      }
  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// XPPM methods

BoutReal BOUTMIN(const BoutReal &a, const BoutReal &b, const BoutReal &c, const BoutReal &d) {
  BoutReal r1 = (a < b) ? a : b;
  BoutReal r2 = (c < d) ? c : d;
  return (r1 < r2) ? r1 : r2;
}

struct Stencil1D {
  // Cell centre values
  BoutReal c, m, p, mm, pp;
  
  // Left and right cell face values
  BoutReal L, R;
};

// First order upwind for testing
void Upwind(Stencil1D &n, const BoutReal h) {
  n.L = n.R = n.c;
}

// Fromm method
void Fromm(Stencil1D &n, const BoutReal h) {
  n.L = n.c - 0.25*(n.p - n.m);
  n.R = n.c + 0.25*(n.p - n.m);
}

BoutReal minmod(BoutReal a, BoutReal b) {
  if( a*b <= 0.0 )
    return 0.0;
  
  if(fabs(a) < fabs(b))
    return a;
  return b;
}

void MinMod(Stencil1D &n, const BoutReal h) {
  // Choose the gradient within the cell
  // as the minimum (smoothest) solution
  BoutReal slope = minmod(n.p - n.c, n.c - n.m);
  n.L = n.c - 0.5*slope; //0.25*(n.p - n.m);
  n.R = n.c + 0.5*slope; //0.25*(n.p - n.m);
}

const Field3D Div_par_FV(const Field3D &f, const Field3D &v) {
  Field3D result = 0.0;
  
  // Calculate in guard cells to avoid having to communicate fluxes
  int ys = mesh->ystart-1;
  int ye = mesh->yend+1;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=ys;j<=ye;j++) {
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // Reconstruct f at the cell faces
        Stencil1D s;
        s.c  = f(i,j,  k);
        s.m  = f(i,j-1,k);
        s.p  = f(i,j+1,k);
        
        //Upwind(s, mesh->dy(i,j));  // 1st order accurate
        //Fromm(s, mesh->dy(i,j));     // 2nd order, some upwinding
        MinMod(s, mesh->dy(i,j));  // Slope limiter
        
        // Calculate velocity at right boundary (y+1/2)
        BoutReal vpar = 0.5*(v(i,j,k) + v(i,j+1,k));
        if(vpar > 0.0) {
          // Out of this cell; use s.R
	  
	  if(mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
	    // Last point in domain: Use mid-point to be consistent with boundary conditions
	    //s.R = 0.5*(s.c + s.p);
	  }
	  BoutReal flux = s.R * vpar * (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j))+ sqrt(mesh->g_22(i,j+1)));
          
          result(i,j,k)   += flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j+1,k) -= flux / (mesh->dy(i,j+1)*mesh->J(i,j+1));
          
        }
        
        // Calculate at left boundary
        vpar = 0.5*(v(i,j,k) + v(i,j-1,k));
        
        if(vpar < 0.0) {
          // Out of this cell; use s.L
	  
          
	  if(mesh->firstY(i) && (j == mesh->ystart) && !mesh->periodicY(i)) {
	    // First point in domain
	    s.L = 0.5*(s.c + s.m);
	  }
          BoutReal flux = s.L * vpar * (mesh->J(i,j) + mesh->J(i,j-1)) / (sqrt(mesh->g_22(i,j)) + sqrt(mesh->g_22(i,j-1)));
          
          result(i,j,k)   -= flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j-1,k) += flux / (mesh->dy(i,j-1)*mesh->J(i,j-1));
        } 
      }      
    }
  return result;
}

const Field3D Div_par_FV_FS(const Field3D &f, const Field3D &v, const Field3D &a) {
  // Finite volume parallel divergence
  Field3D result = 0.0;

  // Calculate in guard cells
  int ys = mesh->ystart-1;
  int ye = mesh->yend+1;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=ys;j<=ye;j++) {
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // Right boundary

        // Reconstruct f at the cell faces
        Stencil1D s;
        s.c  = f(i,j,  k);
        s.m  = f(i,j-1,k);
        s.p  = f(i,j+1,k);
        
        //Upwind(s, mesh->dy(i,j));  // 1st order accurate
        //Fromm(s, mesh->dy(i,j));     // 2nd order, some upwinding
        MinMod(s, mesh->dy(i,j));  // Slope limiter
        
        // Calculate velocity at right boundary (y+1/2)
        BoutReal vpar = 0.5*(v(i,j,k) + v(i,j+1,k));
        
        // Maximum wave speed in the two cells
        BoutReal amax = BOUTMAX(a(i,j,k), a(i,j+1,k));
        
        BoutReal flux;
        
        if(vpar > amax) {
          // Supersonic flow out of this cell
          flux = s.R * vpar;
        }else if(vpar < -amax) {
          // Supersonic flow into this cell
          flux = 0.0;
        }else {
          // Subsonic flow, so a mix of right and left fluxes
          flux = s.R * 0.5*(vpar + amax);
        }
        
        flux *= (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j))+ sqrt(mesh->g_22(i,j+1)));

        result(i,j,k)   += flux / (mesh->dy(i,j)*mesh->J(i,j));
        result(i,j+1,k) -= flux / (mesh->dy(i,j+1)*mesh->J(i,j+1));
        
        // Calculate at left boundary
        vpar = 0.5*(v(i,j,k) + v(i,j-1,k));
        
        // Maximum wave speed in the two cells
        amax = BOUTMAX(a(i,j,k), a(i,j-1,k));
        
        if(vpar < -amax) {
          // Supersonic out of this cell
          flux = s.L * vpar;
        }else if(vpar > amax) {
          // Supersonic into this cell
          flux = 0.0;
        }else {
          flux = s.L * 0.5*(vpar - amax);
        }
        
        flux *= (mesh->J(i,j) + mesh->J(i,j-1)) / (sqrt(mesh->g_22(i,j)) + sqrt(mesh->g_22(i,j-1)));
          
        result(i,j,k)   -= flux / (mesh->dy(i,j)*mesh->J(i,j));
        result(i,j-1,k) += flux / (mesh->dy(i,j-1)*mesh->J(i,j-1));
        
      }
      
    }
  return result;
}

const Field3D Div_par_FV3(const Field3D &f, const Field3D &g, const Field3D &v) { 
  Field3D result = 0.0;

  int ys = mesh->ystart-1;
  int ye = mesh->yend+1;
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=ys;j<=ye;j++) {
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // Reconstruct f at the cell faces
        Stencil1D sf;
        sf.c  = f(i,j,  k);
        sf.m  = f(i,j-1,k);
        sf.p  = f(i,j+1,k);
        
        //Upwind(sf, mesh->dy(i,j));  // 1st order accurate
        //Fromm(sf, mesh->dy(i,j));     // 2nd order, some upwinding
        MinMod(sf, mesh->dy(i,j));  // Slope limiter

        // Reconstruct g at the cell faces
        Stencil1D sg;
        sg.c  = g(i,j,  k);
        sg.m  = g(i,j-1,k);
        sg.p  = g(i,j+1,k);
        
        //Upwind(sf, mesh->dy(i,j));  // 1st order accurate
        //Fromm(sf, mesh->dy(i,j));     // 2nd order, some upwinding
        MinMod(sg, mesh->dy(i,j));  // Slope limiter

        // Calculate velocity at right boundary (y+1/2)
        BoutReal vpar = 0.5*(v(i,j,k) + v(i,j+1,k));
        
        if(vpar > 0.0) {
          // Out of this cell; use s.R
	  
	  if(mesh->lastY(i) && (j == mesh->yend) && !mesh->periodicY(i)) {
	    // Last point in domain: Use mid-point to be consistent with boundary conditions
	    //s.R = 0.5*(s.c + s.p);
	  }
          
          BoutReal flux = 0.5*(sf.R*sg.c + sg.R*sf.c) * vpar * (mesh->J(i,j) + mesh->J(i,j+1)) / (sqrt(mesh->g_22(i,j))+ sqrt(mesh->g_22(i,j+1)));
          
          result(i,j,k)   += flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j+1,k) -= flux / (mesh->dy(i,j+1)*mesh->J(i,j+1));
          
        }
        
        // Calculate at left boundary
        vpar = 0.5*(v(i,j,k) + v(i,j-1,k));
        
        if(vpar < 0.0) {
          // Out of this cell; use s.L
	  
          BoutReal flux = 0.5*(sf.L*sg.c + sf.c*sg.L) * vpar * (mesh->J(i,j) + mesh->J(i,j-1)) / (sqrt(mesh->g_22(i,j)) + sqrt(mesh->g_22(i,j-1)));
          
          result(i,j,k)   -= flux / (mesh->dy(i,j)*mesh->J(i,j));
          result(i,j-1,k) += flux / (mesh->dy(i,j-1)*mesh->J(i,j-1));
        }
        
      }
      
    }
  return result;
}

const Field3D Grad_par_FV(const Field3D &f) {
  Field3D result = 0.0;

  int ys = mesh->ystart-1;
  int ye = mesh->yend+1;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=ys;j<=ye;j++) {
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // Reconstruct f at the cell faces
        Stencil1D s;
        s.c  = f(i,j,  k);
        s.m  = f(i,j-1,k);
        s.p  = f(i,j+1,k);
        
        //Upwind(s, mesh->dy(i,j));  // 1st order accurate
        //Fromm(s, mesh->dy(i,j));     // 2nd order, some upwinding
        MinMod(s, mesh->dy(i,j));  // Slope limiter
        
        result(i,j,k) = (s.R - s.L)/(mesh->dy(i,j) * sqrt(mesh->g_22(i,j)));
      }
    }
  return result;
}

const Field3D D4DY4_FV(const Field3D &d, const Field3D &f) {
  Field3D result = 0.0;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      BoutReal dy3 = SQ(mesh->dy(i,j))*mesh->dy(i,j);
      for(int k=0;k<mesh->ngz-1;k++) {
        
        // 3rd derivative at right boundary
        
        BoutReal d3fdx3 = (
                                f(i,j+2,k)
                           - 3.*f(i,j+1,k)
                           + 3.*f(i,j,  k)
                           -    f(i,j-1,k)
                                ) / dy3;
                           
        BoutReal flux = 0.5*(d(i,j,k) + d(i,j+1,k))*(mesh->J(i,j) + mesh->J(i,j+1)) * d3fdx3;
        
        /*
        if(mesh->lastY() && (j == mesh->yend)) {
          // Boundary? Only if not periodic
          flux = 0.0;
        }
        */

        result(i,j,  k) += flux / (mesh->J(i,j) * mesh->dy(i,j));
        result(i,j+1,k) -= flux / (mesh->J(i,j+1) * mesh->dy(i,j+1));
        
        if(j == mesh->ystart && (!mesh->firstY())) {
          // Left cell boundary, no flux through boundaries
          d3fdx3 = (
                         f(i,j+1,k)
                    - 3.*f(i,j,  k)
                    + 3.*f(i,j-1,k)
                    -    f(i,j-2,k)
                    ) / dy3;
          
          flux = 0.5*(d(i,j,k) + d(i,j-1,k))*(mesh->J(i,j) + mesh->J(i,j-1)) * d3fdx3;
          
          result(i,j,  k) -= flux / (mesh->J(i,j) * mesh->dy(i,j));
          result(i,j-1,k) += flux / (mesh->J(i,j-1) * mesh->dy(i,j-1));
        }
      }
    }
  return result;
}
