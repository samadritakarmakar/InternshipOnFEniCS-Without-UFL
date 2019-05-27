// Copyright (C) 2009-2017 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.

#ifndef __FENICS_SOLID_UTILS_H
#define __FENICS_SOLID_UTILS_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>

namespace fenicssolid
{

  /// Return git commit hash for library
  std::string git_commit_hash();

  // Compute strain (Voigt notation)
  void compute_voigt_strain(Eigen::Matrix<double, 6, 1>& strain,
                            const boost::multi_array<double, 2>& derivatives,
                            const std::vector<double>& expansion_coeffs);

};

#endif
