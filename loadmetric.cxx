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

#include <globals.hxx>
#include <output.hxx>
#include <utils.hxx>
#include <field2d.hxx>
#include <bout/mesh.hxx>

#include "loadmetric.hxx"

using bout::globals::mesh;

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load metric coefficients from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;
  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics

  Coordinates *coord = mesh->getCoordinates();
  
  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if(!mesh->get(dx,   "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    coord->dx = dx; // Only use dpsi if found
  }else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }

  Rxy      /= Lnorm;
  hthe     /= Lnorm;
  sinty    *= SQ(Lnorm)*Bnorm;
  coord->dx /= SQ(Lnorm)*Bnorm;
  
  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  coord->Bxy  /= Bnorm;
  
  // Calculate metric components
  sinty = 0.0;  // I disappears from metric for shifted metric
  
  BoutReal sbp = 1.0; // Sign of Bp
  if(min(Bpxy, true) < 0.0)
    sbp = -1.0;
  
  coord->g11 = SQ(Rxy*Bpxy);
  coord->g22 = 1.0 / SQ(hthe);
  coord->g33 = SQ(sinty)*coord->g11 + SQ(coord->Bxy)/coord->g11;
  coord->g12 = 0.0;
  coord->g13 = -sinty*coord->g11;
  coord->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
  
  coord->J = hthe / Bpxy;
  
  coord->g_11 = 1.0/coord->g11 + SQ((sinty*Rxy));
  coord->g_22 = SQ(coord->Bxy*hthe/Bpxy);
  coord->g_33 = Rxy*Rxy;
  coord->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
  coord->g_13 = sinty*Rxy*Rxy;
  coord->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
  
  coord->geometry();
}
