/*
  Finite volume discretisations of advection and diffusion operators
 
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

#ifndef __DIV_OPS_H__
#define __DIV_OPS_H__

#include <field3d.hxx>

/*!
 * Parallel diffusion (in y)
 *
 * Calculated in terms of fluxes through cell faces. Takes the average
 * coefficient K from the cells either side, and the gradient across the 
 * boundary. Calculated flux is added to one cell, subtracted from the other.
 *
 * Div_par( K Grad_par(f) )
 *
 * @param[in] K The diffusion coefficient
 * @param[in] f The variable to be differentiated
 * @param[in] bndry_flux  Are fluxes calculated through Y boundaries?
 *
 */
const Field3D Div_par_diffusion(const Field3D &k, const Field3D &f, bool bndry_flux=true);

/*!
 * Parallel heat conduction, assuming a heat conduction coefficient
 * K which depends on the temperature Te^2.5
 * 
 * Div_par( K0 Te^2.5 Grad_par(Te) )
 *
 * To calculate K0*Te^2.5 the temperature is averaged from cell centre
 * to cell boundary.
 *
 * @param[in] K0  Constant coefficient in the conductivity
 * @param[in] Te  Temperature
 * @param[in] bndry_flux  Are fluxes through the boundary calculated?
 */
const Field3D Div_par_spitzer(BoutReal K0, const Field3D &Te, bool bndry_flux=true);

/*!
 * Diffusion using upwinding of the conduction coefficient
 *
 * Depending on the sign of the gradient, the value of K from the 
 * "upwind" side is used in calculating the flux, rather than taking 
 * the average of upstream and downstream sides. 
 *
 * Div_par( K Grad_par(f) )
 * 
 * @param[in] K  The diffusion coefficient
 * @param[in] f  The variable which is differentiated
 * @param[in] bndry_flux   Are boundary fluxes calculated?
 */
const Field3D Div_par_diffusion_upwind(const Field3D &K, const Field3D &f, bool bndry_flux=true);

/*!
 * Diffusion in index space
 * 
 * Similar to using Div_par_diffusion(SQ(mesh->dy)*mesh->g_22, f)
 *
 * @param[in] The field to be differentiated
 * @param[in] bndry_flux  Are fluxes through the boundary calculated?
 */
const Field3D Div_par_diffusion_index(const Field3D &f, bool bndry_flux=true);

/*!
 * Added Dissipation scheme (related to Momentum Interpolation)
 *
 * This uses a 3rd-order derivative of the pressure as
 * a correction to the velocity. 
 *
 * This should appear in the form
 * 
 * df/dt = ... + AddedDissipation(N, P, f);
 */
const Field3D AddedDissipation(const Field3D &N, const Field3D &P, const Field3D f, bool bndry_flux=true);

/*!
 * Finite volume parallel divergence
 *
 * Assumes there are (at least) two guard cells (MYG >= 2)
 * 
 * @param[in] f  The field being advected
 * @param[in] v  The advection velocity
 */
const Field3D Div_par_FV(const Field3D &f, const Field3D &v);

/*!
 * Parallel divergence, flux splitting version
 *
 * @param[in] f   The field being advected
 * @param[in] v   The advection velocity
 * @param[in] a   Maximum wave speed. Used to determine the amount of upwinding
 * @param[in] bndry_flux_fixed   If true, calculate the flux by interpolating f and v to the boundary
 * 
 * Split into fluxes with speed v+a and v-a
 */
const Field3D Div_par_FV_FS(const Field3D &f, const Field3D &v, const Field3D &a, bool bndry_flux_fixed=false);

/*!
 * Finite volume parallel divergence
 * 
 * Div_par( f g v )
 *
 * Both f and g are reconstructed to cell boundaries,
 * then the combination is calculated as
 * 
 * (fg)_R = (1/2) ( f_C g_R + f_R g_C )
 *
 */ 
const Field3D Div_par_FV3(const Field3D &f, const Field3D &g, const Field3D &v);

/*!
 * 4th-order derivative
 *
 * Implemented as a flux through cell boundaries, calculated
 * using one-sided 3rd derivative at the boundary.
 *
 * @param[in]  d  Coefficient, averaged from neighbouring cells
 * @param[in]  f  The field being differentiated
 */
const Field3D D4DY4_FV(const Field3D &d, const Field3D &f);

#endif //  __DIV_OPS_H__
