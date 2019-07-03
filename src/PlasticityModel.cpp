// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include "PlasticityModel.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
PlasticityModel::PlasticityModel(double E, double nu)
  : _hardening_parameter(0.0)
{
  // Lame coefficients
  const double lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
  const double mu = E/(2.0*(1.0 + nu));

  // Create elastic tangent
  elastic_tangent = Eigen::Matrix<double, 6, 6>::Constant(0.0);
  elastic_tangent(0, 0) = lambda + 2.0*mu;
  elastic_tangent(1, 1) = lambda + 2.0*mu;
  elastic_tangent(2, 2) = lambda + 2.0*mu;
  elastic_tangent(3, 3) = mu;
  elastic_tangent(4, 4) = mu;
  elastic_tangent(5, 5) = mu;
  elastic_tangent(0, 1) = lambda;
  elastic_tangent(0, 2) = lambda;
  elastic_tangent(1, 0) = lambda;
  elastic_tangent(1, 2) = lambda;
  elastic_tangent(2, 0) = lambda;
  elastic_tangent(2, 1) = lambda;
}
//-----------------------------------------------------------------------------
PlasticityModel::~PlasticityModel()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
double PlasticityModel::hardening_parameter(double const eps_p_eq) const
{
  return _hardening_parameter;
}
// Edited by SAM---------------------------------------------------------------
double PlasticityModel::kappa(double eps_p_eq,
                              const Eigen::Matrix<double, 6, 1>& stress,
                              double lambda_dot) const
{
  return eps_p_eq += lambda_dot;
}
//-----------------------------------------------------------------------------
void PlasticityModel::dg(Eigen::Matrix<double, 6, 1>& dg_dsigma,
                         const Eigen::Matrix<double, 6, 1>& stress) const
{
  // Assume associative flow (dg/dsigma = df/dsigma)
  df(dg_dsigma, stress);
}
//ADDED BY SAM-------------------------------------------------------
void PlasticityModel::df_dq(double &df_dQ,const double &q) const
{
    df_dQ=0;
}

void PlasticityModel::ddg_dsigma_dq(Eigen::Matrix<double, 6, 1> &ddg_dsgma_dq,
                                    const Eigen::Matrix<double, 6, 1> &stress,
                                    const double &q) const
{
    ddg_dsgma_dq.setZero();
}

void PlasticityModel::M(double &m,
                        const Eigen::Matrix<double, 6, 1> &stress,
                        double q) const
{
    m=0;
}

void PlasticityModel::dM_dsigma(Eigen::Matrix<double, 6, 1> &dm_dsgma,
                                const Eigen::Matrix<double, 6, 1> &stress, double q) const
{
    dm_dsgma.setZero();
}

void PlasticityModel::dM_dq(double &dm_dQ,
                            const double &q) const
{
    dm_dQ=0;
}

double PlasticityModel:: q_0() const
{
    return q_0_default;
}

 double PlasticityModel::miscellaneous() const
 {
     return 0.0;
 }

void PlasticityModel:: set_q_0(double q0)
{
    q_0_default=q0;
}
