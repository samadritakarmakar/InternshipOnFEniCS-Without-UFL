// Copyright (C) 2006-2017 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.

#include <vector>
#include <ufc.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/Cell.h>
#include "PlasticityModel.h"
#include "utils.h"
#include "ConstitutiveUpdate.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
ConstitutiveUpdate::ConstitutiveUpdate(
  std::shared_ptr<const dolfin::Function> u,
  std::shared_ptr<const dolfin::FiniteElement> sigma_element,
  std::shared_ptr<const dolfin::GenericDofMap> stress_dofmap,
  std::shared_ptr<const PlasticityModel> plastic_model)
  : _gdim(u->function_space()->mesh()->geometry().dim()),
    _u(u),
    _plastic_model(plastic_model),
    _expansion_coeffs(u->function_space()->dofmap()->max_cell_dimension()),
    _De(_plastic_model->elastic_tangent)
{
    _eps_p = std::make_shared<HistoryData>(u->function_space()->mesh(), sigma_element, 6);
    _eps_p_equiv = std::make_shared<HistoryData>(u->function_space()->mesh(), sigma_element, 1);
  // Get stress UFC element
  auto ufc_element_sigma = sigma_element->ufc_element();
  dolfin_assert(ufc_element_sigma);

  // Get stress dof dimension data
  const std::size_t dim = ufc_element_sigma->space_dimension();
  const std::size_t gdim = ufc_element_sigma->geometric_dimension();
  const std::size_t tdim = ufc_element_sigma->topological_dimension();

  // Get quadrature point coordinates on reference element
  //boost::multi_array<double, 2> ip_coordinates(boost::extents[dim][tdim]);
  _X.resize(boost::extents[dim][tdim]);
  ufc_element_sigma->tabulate_reference_dof_coordinates(_X.data());

  // Get displacement UFC element (single component)
  const dolfin::FiniteElement& u_element_new = *(*_u)[0].function_space()->element();
  auto ufc_element_u = u_element_new.ufc_element();
  dolfin_assert(ufc_element_u);

  // Compute basis function derivatives on reference element and store
  const std::size_t dim_u = ufc_element_u->space_dimension();
  const std::size_t gdim_u = ufc_element_u->geometric_dimension();
  const std::size_t num_points = _X.shape()[0];
  const std::size_t value_size = ufc_element_u->reference_value_size();

  //boost::multi_array<double, 3> derivatives(boost::extents[num_points][dim_u][tdim]);
  this->_derivs_reference.resize(boost::extents[num_points][dim_u][tdim]);
  ufc_element_u->evaluate_reference_basis_derivatives(this->_derivs_reference.data(),
                                                      1, _X.shape()[0], _X.data());

  // Compute number of quadrature points per cell
  const std::size_t num_ip_dofs = sigma_element->value_dimension(0);
  _num_ip_per_cell = sigma_element->space_dimension()/num_ip_dofs;

  // Resize _w_stress/tangent
  const::std::size_t sigma_element_dim = sigma_element->space_dimension();
  _w_stress.resize(num_ip_dofs*_num_ip_per_cell);
  _w_tangent.resize(num_ip_dofs*num_ip_dofs*_num_ip_per_cell);

  const std::size_t num_cells = u->function_space()->mesh()->num_cells();
  _plastic_last.resize(boost::extents[num_cells][_num_ip_per_cell]);
  std::fill(_plastic_last.data(),
            _plastic_last.data() + _plastic_last.num_elements(),
            false);
}
//-----------------------------------------------------------------------------
ConstitutiveUpdate::~ConstitutiveUpdate()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ConstitutiveUpdate::update(const dolfin::Cell& cell,
                                const double* vertex_coordinates)
{
  const std::size_t cell_index = cell.index();

  // Get solution dofs on cell
  const auto dofs = _u->function_space()->dofmap()->cell_dofs(cell_index);

  // Get expansion coefficients on cell
  dolfin_assert(_expansion_coeffs.size()
                ==  _u->function_space()->dofmap()->max_cell_dimension());
  _u->vector()->get_local(_expansion_coeffs.data(), dofs.size(), dofs.data());

  Eigen::Matrix<double, 6, 6> cons_tangent;
  Eigen::Matrix<double, 6, 1> strain, strain_p, trial_stress;
  Eigen::Matrix<double, 1, 1> strain_p_eq;

  // Call functions on UFC coordinate mapping object
  /*
  void compute_reference_geometry( double * X, double * J, double * detJ, double * K,
                                   std::size_t num_points,
                                   const double * x, const double * coordinate_dofs,
                                   int cell_orientation) const final override
  */

  // Compute geometry mapping
  double J[9];
  double detJ;
  double K[9];
  if (_gdim == 2)
  {
    compute_jacobian_triangle_2d(J, vertex_coordinates);
    compute_jacobian_inverse_triangle_2d(K, detJ, J);
  }
  else if (_gdim == 3)
  {
    compute_jacobian_tetrahedron_3d(J, vertex_coordinates);
    compute_jacobian_inverse_tetrahedron_3d(K, detJ, J);
  }

  // Duplicate data at each point
  auto num_points = _derivs_reference.shape()[0];
  boost::multi_array<double, 2> _J(boost::extents[num_points][9]);
  std::vector<double> _detJ(num_points, detJ);
  boost::multi_array<double, 2> _K(boost::extents[num_points][9]);
  for (std::size_t i = 0; i < num_points; ++i)
  {
    for (std::size_t j = 0; j < 9; ++j)
    {
      _J[i][j] = J[j];
      _K[i][j] = K[j];
    }
  }

  // Get displacement UFC element (single component)
  const dolfin::FiniteElement& u_element_new = *(*_u)[0].function_space()->element();
  auto ufc_element_u = u_element_new.ufc_element();
  dolfin_assert(ufc_element_u);

  // Push derivatives forward to current physical cell
  auto dim_u = _derivs_reference.shape()[1];
  auto gdim = _gdim;
  boost::multi_array<double, 3> derivs_physical(boost::extents[num_points][dim_u][gdim]);
  ufc_element_u->transform_reference_basis_derivatives(derivs_physical.data(), 1, num_points,
                                                       _derivs_reference.data(), _X.data(),
                                                       _J.data(), _detJ.data(), _K.data(), 0);

  // Loop over quadrature points
  for (std::size_t ip = 0; ip < _num_ip_per_cell; ip++)
  {
    // Assume elastic tangent
    cons_tangent = _De;

    // Get plastic strain from previous converged time step
    _eps_p->get_old_values(cell_index, ip, strain_p);

    // Compute strain on physical cell (Voigt notation)
    compute_voigt_strain(strain, derivs_physical[ip], _expansion_coeffs);

    // Compute trial stress
    trial_stress = _De*(strain - strain_p);

    // Get equivalent plastic strain from previous converged time
    // step
    _eps_p_equiv->get_old_values(cell_index, ip, strain_p_eq);

    // Testing trial stresses, if yielding occurs the stresses are
    // mapped back onto the yield surface, and the updated parameters
    // are returned.
    const bool active = _plastic_last[cell_index][ip];
    return_mapping.closest_point_projection(_plastic_model, cons_tangent,
                                            trial_stress, strain_p,
                                            strain_p_eq(0),
                                            active);
    _plastic_last[cell_index][ip] = false;

    // Update plastic strain history for current load step
    _eps_p->set_new_values(cell_index, ip, strain_p);

    // Update equivalent plastic strain for current load step
    _eps_p_equiv->set_new_values(cell_index, ip, strain_p_eq);

    // Copy data into structures
    if (_gdim == 3)
    {
      for (std::size_t d = 0; d < 6; ++d)
        _w_stress[_num_ip_per_cell*d  + ip] = trial_stress(d);

      for (std::size_t d0 = 0; d0 < 6; ++d0)
      {
        for (std::size_t d1 = 0; d1 < 6; ++d1)
        {
          const std::size_t pos = d0*6 + d1;
          _w_tangent[_num_ip_per_cell*pos  + ip] = cons_tangent(d0, d1);
        }
      }
    }
    else
    {
      _w_stress[_num_ip_per_cell*0  + ip] = trial_stress(0);
      _w_stress[_num_ip_per_cell*1  + ip] = trial_stress(1);
      _w_stress[_num_ip_per_cell*2  + ip] = trial_stress(3);

      _w_tangent[_num_ip_per_cell*0  + ip] = cons_tangent(0, 0);
      _w_tangent[_num_ip_per_cell*1  + ip] = cons_tangent(0, 1);
      _w_tangent[_num_ip_per_cell*2  + ip] = cons_tangent(0, 3);

      _w_tangent[_num_ip_per_cell*3  + ip] = cons_tangent(1, 0);
      _w_tangent[_num_ip_per_cell*4  + ip] = cons_tangent(1, 1);
      _w_tangent[_num_ip_per_cell*5  + ip] = cons_tangent(1, 3);

      _w_tangent[_num_ip_per_cell*6  + ip] = cons_tangent(3, 0);
      _w_tangent[_num_ip_per_cell*7  + ip] = cons_tangent(3, 1);
      _w_tangent[_num_ip_per_cell*8  + ip] = cons_tangent(3, 3);
    }
  }
}
//-----------------------------------------------------------------------------
void ConstitutiveUpdate::update_history()
{
  // Update plastic elements
  const boost::multi_array<double, 3>& old_eps = _eps_p_equiv->old_data();
  const boost::multi_array<double, 3>& new_eps = _eps_p_equiv->current_data();

  const std::size_t num_cells = _plastic_last.shape()[0];
  const std::size_t ip_per_cell = _plastic_last.shape()[1];
  dolfin_assert(old_eps.shape()[0] == num_cells);

  for (std::size_t c = 0; c < num_cells; ++c)
  {
    for (std::size_t p = 0; p < ip_per_cell; ++p)
    {
      if ((new_eps[c][p][0] - old_eps[c][p][0] > 0.0))
        _plastic_last[c][p] = true;
      else
        _plastic_last[c][p] = false;
    }
  }

  _eps_p->update_history();
  _eps_p_equiv->update_history();
}
//-----------------------------------------------------------------------------
