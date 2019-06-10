// Copyright (C) 2006-2010 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2010-01-03

#include <iostream>

#include <dolfin/common/constants.h>
#include "PlasticityModel.h"
#include "ReturnMapping.h"

using namespace fenicssolid;

//-----------------------------------------------------------------------------
ReturnMapping::ReturnMapping(const unsigned int maxit) : _maxit(maxit)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ReturnMapping::~ReturnMapping()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::pair<bool, unsigned int>
ReturnMapping::closest_point_projection(std::shared_ptr<const PlasticityModel> plastic_model,
                                        Eigen::Matrix<double, 6, 6>& D,
                                        Eigen::Matrix<double, 6, 1>& trial_stresses,
                                        Eigen::Matrix<double, 6, 1>& plastic_strain,
                                        double& equivalent_plastic_strain, double &qn,
                                        bool use_plastic_tangent)
{
  // Initialise variables
  double delta_lambda = 0.0;
  bool plastic_flag = false;
  unsigned int num_iterations = 0;

  // Work arrays (statically allocated)
  Eigen::Matrix<double, 6, 6> Q1;
  Eigen::Matrix<double, 6, 6>  /*R,*/ ddg_ddsigma, inverse_Q;
  // ADDED BY SAM
  Eigen::Matrix<double, 7, 7> Q, XI_2, C_dud, D_dud;
  Eigen::Matrix<double, 7, 7> Qinv;
  Eigen::Matrix<double, 1, 1> df_dQ, q_current, q_dot, q_n, M_current, dM_dQ, q_residual, Q4;
  Eigen::Matrix<double, 6, 1> ddg_dsgma_dq, Q2,dM_dsgma, Q3, zeroRow;
  Eigen::Matrix<double, 7, 1> residualCombined, df_dCombined, f_residualCombined, C_combined, delta_sigma_q, XI, zero;
  Eigen::Matrix<double, 1, 7> Rn;
  //----------------------------------------

  Eigen::Matrix<double, 6, 1> Rm, RinvQ;

  // Variables for return mapping
  Eigen::Matrix<double, 6, 1> df_dsigma, dg_dsigma,
    sigma_current,
    sigma_dot, sigma_residual;


  // Trial stress
  sigma_current = trial_stresses;


  // Auxiliary variables to speed up return mapping

  // Elastic tangent
  const Eigen::Matrix<double, 6, 6>& De = plastic_model->elastic_tangent;

  // ADDED BY SAM---------------------------------------------------------
  zero.setZero(); zeroRow.setZero();
  C_dud<< De, zeroRow,
          zeroRow.transpose(), zero;

  // Evaluate hardening parameter
  double hardening_parameter
    = plastic_model->hardening_parameter(equivalent_plastic_strain);

  // Evaluate yield function (trial stresses)
  double residual_f = plastic_model->f(sigma_current, equivalent_plastic_strain);

  // Check for yielding
  if (residual_f/sigma_current.norm() > 1.0e-12)
  {
    // Compute normal vectors to yield surface and plastic potential
    plastic_model->df(df_dsigma, sigma_current);
    plastic_model->dg(dg_dsigma, sigma_current);
    plastic_flag = true;
    // ADDED BY SAM------------------------------------------------------------
    q_n(0,0)=qn;
    q_current=q_n;
    //------------------------------------------------------------------------

    // Perform Newton iterations to project stress onto yield surface
    while (std::abs(residual_f)/sigma_current.norm() > 1e-12)
    {
      num_iterations++;
      if (num_iterations > _maxit)
        dolfin::error("Return mapping iterations > %d.", _maxit);

      // Reset sigma_residual in first step
      if (num_iterations == 1)
      {
          sigma_residual.setZero();
          q_residual.setZero();
      }

      // Compute second derivative (with respect to stresses) of
      // plastic potential
      plastic_model->ddg(ddg_ddsigma, sigma_current);
      // ADDED BY SAM
      plastic_model->M(M_current, sigma_current, q_current);
      plastic_model->df_dq(df_dQ, q_current);
      plastic_model->ddg_dsigma_dq(ddg_dsgma_dq, sigma_current, q_current);
      plastic_model->dM_dsigma(dM_dsgma, sigma_current);
      plastic_model->dM_dq(dM_dQ, q_current);

      //q_n=q_current(0,0);
      //-------------------------------------------------------------

      // Compute auxiliary matrix Q
      Q1 = Eigen::Matrix<double, 6, 6>::Identity() + delta_lambda*De*ddg_ddsigma;

      // ADDED BY SAM
      Q2 = delta_lambda*De*ddg_dsgma_dq;
      Q3 = delta_lambda*hardening_parameter*dM_dsgma;
      Q4 = Eigen::Matrix<double, 1, 1>::Identity()+delta_lambda*hardening_parameter*dM_dQ;
      Qinv<<Q1,Q2
      ,Q3.transpose(),Q4;
      //------------------------------------------------------------------------
      //Edited by SAM------------------------------------------------------------
      // Invert Qinv
      Q = Qinv.inverse();
      //ADDED BY SAM
      C_combined<<De*dg_dsigma,
                  hardening_parameter*M_current;
      // Compute auxiliary matrix R
      XI = Q*C_combined;
      // ADDED BY SAM------------------------------------------------------------
      df_dCombined<<df_dsigma,
              df_dQ;
      residualCombined<<sigma_residual,
                        q_residual;
      // lambda_dot, rate of plastic multiplier
      //const double residual_tmp = residual_f - sigma_residual.dot(Qinv*df_dsigma);
      //double lambda_dot = residual_tmp/(df_dsigma.dot(XI*dg_dsigma)
      //                                  + hardening_parameter);
      const double residual_tmp = residual_f- df_dCombined.dot(Q*residualCombined);
      double lambda_dot =residual_tmp/(df_dCombined.dot(Q*C_combined));
      //---------------------------------------------------------------------------

      // Compute stress increment
      // FIXME: (GNW) is the below correction? In uBLAS it was prod(x, A)
      //sigma_dot = Q.transpose()*(-lambda_dot*De*dg_dsigma -sigma_residual);
      // ADDED BY SAM------------------------------------------------------------
      delta_sigma_q=-Q*(residualCombined+C_combined*lambda_dot);
      sigma_dot=delta_sigma_q.block(0,0,5,0);
      q_dot=delta_sigma_q.block(6,0,6,0);
      //-------------------------------------------------------------------------


      // Increment plastic multiplier
      delta_lambda += lambda_dot;

      // Update current stress state
      sigma_current += sigma_dot;
      // ADDED BY SAM------------------------------------------------------------
      q_current += q_dot;
      //-------------------------------------------------------------------------

      //Edited by SAM------------------------------------------------------------
      // Update equivalent plastic strain
      equivalent_plastic_strain
        = plastic_model->kappa(equivalent_plastic_strain, sigma_current,
                              lambda_dot);
      //-------------------------------------------------------------------------

      // Compute hardening parameter
      hardening_parameter
        = plastic_model->hardening_parameter(equivalent_plastic_strain);

      // Evaluate yield function at new stress state
      residual_f = plastic_model->f(sigma_current, equivalent_plastic_strain);

      // Compute normal to yield surface at new stress state
      plastic_model->df(df_dsigma, sigma_current);

      // Compute normal normal to plastic potential at new stress state
      plastic_model->dg(dg_dsigma, sigma_current);

      // Compute residual vector
      sigma_residual = sigma_current
        - (trial_stresses - delta_lambda*De*dg_dsigma);
      // ADDED BY SAM------------------------------------------------------------
      q_residual = q_current -q_n+ delta_lambda*hardening_parameter*M_current;
      //-------------------------------------------------------------------------
    }

    // Update matrices
    plastic_model->ddg(ddg_ddsigma, sigma_current);
    // ADDED BY SAM------------------------------------------------------------
    plastic_model->M(M_current, sigma_current, q_current);
    plastic_model->df_dq(df_dQ, q_current);
    plastic_model->ddg_dsigma_dq(ddg_dsgma_dq, sigma_current, q_current);
    plastic_model->dM_dsigma(dM_dsgma, sigma_current);
    plastic_model->dM_dq(dM_dQ, q_current);
    //-------------------------------------------------------------------------
    // Compute matrix Q
    Q1 = Eigen::Matrix<double, 6, 6>::Identity() + delta_lambda*De*ddg_ddsigma;
    // ADDED BY SAM
    Q2 = delta_lambda*De*ddg_dsgma_dq;
    Q3 = delta_lambda*hardening_parameter*dM_dsgma;
    Q4 = Eigen::Matrix<double, 1, 1>::Identity()+delta_lambda*hardening_parameter*dM_dQ;

    Qinv<<Q1,Q2
    ,Q3.transpose(),Q4;

    // Invert Q
    Q = Qinv.inverse();

    //ADDED BY SAM
    C_combined<<De*dg_dsigma,
                hardening_parameter*M_current;

    // Compute auxiliary matrix XI
    XI = Q*C_combined;
    df_dCombined<<df_dsigma,
                  df_dQ;
    // Compute matrix R and vector Rn
    XI_2 = Q*C_dud;
    Rn = df_dCombined.transpose()*XI_2;


    // Compute consistent tangent operator
    //D = XI - XI*(dg_dsigma*Rn.transpose())/(df_dsigma.dot(XI*dg_dsigma)
    //                                      + hardening_parameter);
    D_dud = XI_2 - (XI*Rn)/(df_dCombined.dot(XI));
    D=D_dud.block(0,0,6,6);

    // Stresses for next Newton iteration, trial stresses are
    // overwritten by current stresses
    trial_stresses = sigma_current;
    //ADDED BY SAM-------------------------------------------
    q_n=q_current;
    qn=q_n(0,0);
    //-------------------------------------------------------

    // Update plastic strains
    plastic_strain += delta_lambda*dg_dsigma;
  }
  else if (use_plastic_tangent)
  {
    // Compute continuum tangent operator
    plastic_model->df(df_dsigma, sigma_current);
    plastic_model->dg(dg_dsigma, sigma_current);
    //ADDED BY SAM
    plastic_model->M(M_current, sigma_current, q_current);
    plastic_model->df_dq(df_dQ, q_current);
    //Rn = De*df_dsigma;
    //EDITED BY SAM
    const double denom = df_dsigma.dot(De*dg_dsigma) + hardening_parameter*df_dQ*M_current;
    //---------------------------------------------------------------------------
    D = De - De*(dg_dsigma*(De*df_dsigma).transpose())/denom;
  }
  else
  {
    // Use elastic tangent
    D = De;
  }

  return std::make_pair(plastic_flag, num_iterations);
}
//-----------------------------------------------------------------------------
