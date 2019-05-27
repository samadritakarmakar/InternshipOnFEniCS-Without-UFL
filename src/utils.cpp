// Copyright (C) 2009-2012 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2009-10-02
// Last changed: 2012-07-17

#include <iostream>

#include <boost/multi_array.hpp>
#include <ufc.h>
#include <ufc_geometry.h>
#include <dolfin/fem/FiniteElement.h>
#include "utils.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
std::string fenicssolid::git_commit_hash()
{
  return FENICSSOLID_GIT_COMMIT_HASH;
}
//-----------------------------------------------------------------------------
void fenicssolid::compute_voigt_strain(Eigen::Matrix<double, 6, 1>& strain,
                                       const boost::multi_array<double, 2>& derivatives,
                                       const std::vector<double>& expansion_coeff)
{
  // Zero strain vector
  strain.setZero();

  const std::size_t gdim = derivatives.shape()[1];
  const std::size_t space_dim = derivatives.shape()[0];

  for (unsigned int dim = 0; dim < space_dim; dim++)
  {
    // Ux,x (eps_xx)
    strain(0) += derivatives[dim][0]*expansion_coeff[dim];

    // Uy,y (eps_yy)
    strain(1) += derivatives[dim][1]*expansion_coeff[space_dim + dim];

    // Ux,y + Uy,x (gamma_xy)
    strain(3) += derivatives[dim][1]*expansion_coeff[dim]
      + derivatives[dim][0]*expansion_coeff[space_dim + dim];

    // Add 3D strains
    if (gdim == 3)
    {
      // Uz,z (eps_zz)
      strain(2) += derivatives[dim][2]*expansion_coeff[2*space_dim + dim];

      // Ux,z + Uz,x (gamma_xz)
      strain(4) += derivatives[dim][2]*expansion_coeff[dim]
        + derivatives[dim][0]*expansion_coeff[2*space_dim + dim];

      // Uy,z + Uz,y (gamma_yz)
      strain(5) += derivatives[dim][2]*expansion_coeff[space_dim + dim]
        + derivatives[dim][1]*expansion_coeff[2*space_dim + dim];
    }
  }
}
//-----------------------------------------------------------------------------
