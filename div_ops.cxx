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

using bout::globals::mesh;

const Field3D Div_par_diffusion(const Field3D &K, const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  Coordinates *coord = mesh->getCoordinates();
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal c = 0.5*(K(i,j,k) + K(i,j+1,k)); // K at the upper boundary
        BoutReal J = 0.5*(coord->J(i,j) + coord->J(i,j+1)); // Jacobian at boundary
        
        //if((i == mesh->xstart) && (j == mesh->ystart) && (k==0))
        //  output << "c = " << c << endl;
        
        BoutReal g_22 = 0.5*(coord->g_22(i,j) + coord->g_22(i,j+1));
        
        BoutReal gradient = 2.*(f(i,j+1,k) - f(i,j,k)) / (coord->dy(i,j) + coord->dy(i,j+1));
        
        BoutReal flux = c * J * gradient / g_22;
        
        result(i,j,k) += flux / (coord->dy(i,j) * coord->J(i,j));
        result(i,j+1,k) -= flux / (coord->dy(i,j+1) * coord->J(i,j+1));
      }
  return result;
}

const Field3D Div_par_spitzer(BoutReal K0, const Field3D &Te, bool bndry_flux) {
  Field3D result;
  result = 0.0;

  Coordinates *coord = mesh->getCoordinates();

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        // Calculate flux at upper surface

        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;

          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }

	BoutReal Te0 = 0.5*(Te(i,j,k) + Te(i,j+1,k)); // Te at the upper boundary
        BoutReal K = K0*pow(Te0,2.5);
        BoutReal J = 0.5*(coord->J(i,j) + coord->J(i,j+1)); // Jacobian at boundary

        BoutReal g_22 = 0.5*(coord->g_22(i,j) + coord->g_22(i,j+1));

        BoutReal gradient = 2.*(Te(i,j+1,k) - Te(i,j,k)) / (coord->dy(i,j) + coord->dy(i,j+1));

        BoutReal flux = K * J * gradient / g_22;

        result(i,j,k) += flux / (coord->dy(i,j) * coord->J(i,j));
        result(i,j+1,k) -= flux / (coord->dy(i,j+1) * coord->J(i,j+1));
      }
  return result;
}

const Field3D Div_par_diffusion_upwind(const Field3D &K, const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  Coordinates *coord = mesh->getCoordinates();

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal J = 0.5*(coord->J(i,j) + coord->J(i,j+1)); // Jacobian at boundary
        
        BoutReal g_22 = 0.5*(coord->g_22(i,j) + coord->g_22(i,j+1));
        
        BoutReal gradient = 2.*(f(i,j+1,k) - f(i,j,k)) / (coord->dy(i,j) + coord->dy(i,j+1));
        
        BoutReal c; // K at the upper boundary
        if(gradient > 0.0) {
          c = K(i,j+1,k);
        }else {
          c = K(i,j,k);
        }
        
        BoutReal flux = c * J * gradient / g_22;
        
        result(i,j,k) += flux / (coord->dy(i,j) * coord->J(i,j));
        result(i,j+1,k) -= flux / (coord->dy(i,j+1) * coord->J(i,j+1));
      }
  return result;
}

const Field3D Div_par_diffusion_index(const Field3D &f, bool bndry_flux) {
  Field3D result;
  result = 0.0;
  
  Coordinates *coord = mesh->getCoordinates();

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        // Calculate flux at upper surface
        
        if(!bndry_flux && !mesh->periodicY(i)) {
          if((j == mesh->yend) && mesh->lastY(i))
            continue;
        
          if((j == mesh->ystart-1) && mesh->firstY(i))
            continue;
        }
        BoutReal J = 0.5*(coord->J(i,j) + coord->J(i,j+1)); // Jacobian at boundary
        
        BoutReal gradient = f(i,j+1,k) - f(i,j,k);
        
        BoutReal flux = J * gradient;
        
        result(i,j,k) += flux / coord->J(i,j);
        result(i,j+1,k) -= flux / coord->J(i,j+1);
      }
  return result;
}

const Field3D AddedDissipation(const Field3D &N, const Field3D &P, const Field3D f, bool bndry_flux) {
  Field3D result = 0.0;
  
  Coordinates *coord = mesh->getCoordinates();

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart-1;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
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

        BoutReal flux = var * v * (coord->J(i,j) + coord->J(i,j+1)) / (sqrt(coord->g_22(i,j))+ sqrt(coord->g_22(i,j+1)));
        
        result(i,j,k) -= flux / (coord->dy(i,j) * coord->J(i,j));
        result(i,j+1,k) += flux / (coord->dy(i,j+1) * coord->J(i,j+1));
      }
  return result;
}
