// Copyright (C) 2006-2011 Kristian Oelgaard and Garth N. Wells.
// Licensed under the GNU LGPL Version 3.
//
// First added:  2006-11-13
// Last changed: 2011-02-06

#include <dolfin.h>
#include <FenicsSolidMechanics.h>
#include <WriteToOctaveASCII.hpp>

// Switch between linear and quadratic elements
//#include "../forms/p2_forms/Plas3DWithForce.h"
//#include "../forms/p1_forms/Plas3DWithForce.h"
#include "Plas3DWithForce.h"
using namespace dolfin;
// Displacement right end
/*class DBval : public Expression
{
  public:

    DBval(const double& t) : t(t), Expression(3) {}

    void eval(Array<double>& values, const Array<double>& x) const
    {
      // Stretch
      values[0] = 0;
      values[1] = 0;
      values[2] = -t;

    }
    const double& t;
};*/


class Traction : public Expression
{
  public:

    Traction(const double& traction) : traction(traction), Expression(3) {}

    void eval(Array<double>& values, const Array<double>& x) const
    {
      // Stretch
      values[0] = 0;
      values[1] = 0;
      values[2] = -traction*(x[2] >= 0.05 - DOLFIN_EPS);

    }
    const double& traction;
};

class NeumannBC : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      //return std::abs(x[2] - 0.05) < DOLFIN_EPS;
      return x[2] >= 0.05 - DOLFIN_EPS;

  }
};

// Sub domain for Dirichlet boundary condition
/*class DirichletBoundaryX1 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      //return std::abs(x[2] - 0.05) < DOLFIN_EPS;
      return x[2] >= 0.05 - DOLFIN_EPS;

  }
};*/

// Sub domain for Dirichlet boundary condition
class DirichletBoundaryX0 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      return (x[2] < DOLFIN_EPS );
  }
};

class DirichletBoundaryX2 : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
      return on_boundary && (std::sqrt(x[0]*x[0]+x[1]*x[1]) > 0.0495-DOLFIN_EPS);
      //return on_boundary &&(x[2]>DOLFIN_EPS && x[2]<0.05-DOLFIN_EPS);
  }
};

int main()
{
  Timer timer("Total plasicity solver time");

  //dolfin::parameters["reorder_dofs_serial"] = false;
  //dolfin::parameters["dof_ordering_library"] = "SCOTCH";
  //dolfin::parameters["linear_algebra_backend"] = "Epetra";

  // Create mesh
  auto mesh=std::make_shared<Mesh>("cylinder.xml");

  // Young's modulus and Poisson's ratio
  double E=3.46e6, nu=0.2;

  // Source term, RHS
  auto f = std::make_shared<Constant>(0.0, 0.0, 0.0);

  // Function spaces
  auto V = std::make_shared<Plas3DWithForce::CoefficientSpace_f>(mesh);
  dolfin::cout << "Number of dofs: " << V->dim() << dolfin::endl;

  // Extract elements for stress and tangent
  std::shared_ptr<const FiniteElement> element_t;
  {
    Plas3DWithForce::CoefficientSpace_t Vt(mesh);
    element_t = Vt.element();
  }

  auto Vs = std::make_shared<Plas3DWithForce::CoefficientSpace_s>(mesh);
  auto element_s = Vs->element();

  // Load-disp info
  unsigned int step = 0;
  //unsigned int steps = dt0_steps + dt1_steps + dt2_steps + 1;
  //double u_target=0.0035;
  unsigned int steps = 60;
  int startSuffusionAtStep=20;
  int stopSuffusionAtStep=40;
  double tractionTarget=6.3e4;
  // Time parameter
  double t = 0.0;
  // Elastic time step, always one step.
  //double Edt  = u_target/steps;
  // Load region 0, time step and number of steps
 // double dt0 = Edt;
  // Traction Step
  double Tdt = tractionTarget/(steps-double(stopSuffusionAtStep-startSuffusionAtStep));
  double Tdt0  = Tdt;

  // Create boundary conditions (use SubSpace to apply simply
  // supported BCs)
  auto zero = std::make_shared<Constant>(0.0, 0.0, 0.0);
  //auto val =std::make_shared<DBval>(t);
  auto TractionObject = std::make_shared<Traction>(t);
  auto zero0 = std::make_shared<Constant>(0.0);
  auto Vx= V->sub(0);
  auto Vy= V->sub(1);
  auto Vz= V->sub(2);

  auto dbX0 = std::make_shared<DirichletBoundaryX0>();
  //auto dbX1 = std::make_shared<DirichletBoundaryX1>();
  auto dbX2 = std::make_shared<DirichletBoundaryX2>();
  auto bcX0 = std::make_shared<DirichletBC>(V, zero, dbX0);
  //auto bcX1 = std::make_shared<DirichletBC>(V, val, dbX1);
  auto bcX2 = std::make_shared<DirichletBC>(Vx, zero0, dbX2);
  auto bcX3 = std::make_shared<DirichletBC>(Vy, zero0, dbX2);

  //std::vector<std::shared_ptr<const DirichletBC>> bcs = {bcX0, bcX1, bcX2, bcX3};
  std::vector<std::shared_ptr<const DirichletBC>> bcs = {bcX0, bcX2, bcX3};

  // Slope of hardening (linear) and hardening parameter
  const double E_t = 0.1*E;

  // Yield stress
  //const double yield_stress = 235.0e6;

  // Solution function
  auto u = std::make_shared<Function>(V);

  // Object of class von Mises
  double beta=1.2, phiDegree=38.5, betaP=0.01, varKappa=0, Pc=50e3, varP=0;
  double Pdash0 = 0.2*Pc;

  double Phi_er_dot=0.001;
  auto SncF = std::make_shared<fenicssolid::SinfoniettaClassica>(E, nu, beta, phiDegree, betaP, varKappa, Pc, varP, Pdash0);

  // Constituive update
  auto constitutive_update
    = std::make_shared<fenicssolid::ConstitutiveUpdate>(u, element_s,
                                                       Vs->dofmap(), SncF);

  // Create forms and attach functions
  auto tangent
    = std::make_shared<fenicssolid::QuadratureFunction>(mesh, element_t,
                                                        constitutive_update,
                                                        constitutive_update->w_tangent());

  auto a = std::make_shared<Plas3DWithForce::Form_a>(V, V);
  a->t = tangent;

  auto L = std::make_shared<Plas3DWithForce::Form_L>(V);
  L->f = f;
  auto stress = std::make_shared<fenicssolid::QuadratureFunction>(mesh, element_s,
                                         constitutive_update->w_stress());
  L->s = stress;
  L->traction= TractionObject;
  L->ds;


  //ADDED by SAM
  //auto V2 = std::make_shared<Plas3DWithForce::CoefficientSpace_s>(mesh);
  //auto V2 =std::make_shared<Plas3DWithForce::CoefficientSpace_fstrss>(mesh); //ADDED BY Q
  //auto V3 =std::make_shared<Plas3DWithForce::CoefficientSpace_s2>(mesh);
  //auto strss =std::make_shared<Function>(V2);
  //auto aStrss =std::make_shared<Plas3DWithForce::Form_aStrss>(V2,V2);
  //auto LStrss = std::make_shared<Plas3DWithForce::Form_L_Strss>(V2);
  //LStrss->s2 =stress;
  //auto f2 = std::make_shared<Constant>(0.0, 0.0, 0.0, );
  std::vector<std::size_t> value_shape;
  value_shape.push_back(3);
  value_shape.push_back(3);
  std::vector<double> values;
  for (int i=0; i<9; i++)
      values.push_back(0.0);

  auto f2=std::make_shared<Constant>(value_shape, values); //ADDED BY Q
  //LStrss->fstrss =f2 ; //ADDED BY Q
//
//ADDED BY SAM For Eps
  auto V3= std::make_shared<Plas3DWithForce::CoefficientSpace_fEps>(mesh);
  auto Eps= std::make_shared<Function>(V3);
  auto aEps = std::make_shared<Plas3DWithForce::Form_aEps>(V3,V3);
  auto LEps = std:: make_shared<Plas3DWithForce::Form_LEps>(V3);
  LEps->u2=u;
  auto fEps= std::make_shared<Constant>(value_shape, values);
  LEps->fEps =fEps;
//

//ADDED BY SAM For Eps_p
  auto V4= std::make_shared<Plas3DWithForce::CoefficientSpace_fEps_p>(mesh);
  auto Eps_p= std::make_shared<Function>(V4);
  auto aEps_p = std::make_shared<Plas3DWithForce::Form_aEps_p>(V4,V4);
  auto LEps_p =std::make_shared<Plas3DWithForce::Form_LEps_p>(V4);
  LEps_p->eps_p=constitutive_update->eps_p();
  auto fEps_p= std::make_shared<Constant>(value_shape, values);
  LEps_p->fEps_p =fEps_p;
  //
//ADDED BY SAM For stress
  auto V5= std::make_shared<Plas3DWithForce::CoefficientSpace_fstrss>(mesh);
  auto strss2 =std::make_shared<Function>(V5);
  auto aStrss =std::make_shared<Plas3DWithForce::Form_aStrss2>(V5, V5);
  auto LStress = std::make_shared<Plas3DWithForce::Form_L_Strss2>(V5);
  LStress->u2 = u;
  LStress->eps_p=constitutive_update->eps_p();
  LStress->fstrss=f2;

  auto ConstZeroScalar=std::make_shared<Constant>(0);
  //ADDED BY Q for p dash

  auto V6=std::make_shared<Plas3DWithForce::CoefficientSpace_f_zero>(mesh);
  auto Pdash= std::make_shared<Function>(V6);
  auto aPdash = std::make_shared<Plas3DWithForce::Form_aPdash>(V6,V6);
  auto LPdash =std::make_shared<Plas3DWithForce::Form_LPdash>(V6);
  LPdash->u2 = u;
  LPdash->eps_p=constitutive_update->eps_p();
  LPdash->f_zero=std::make_shared<Constant>(0);

  auto Qdash= std::make_shared<Function>(V6);
  auto aQdash = std::make_shared<Plas3DWithForce::Form_aQdash>(V6,V6);
  auto LQdash =std::make_shared<Plas3DWithForce::Form_LQdash>(V6);
  LQdash->u2 = u;
  LQdash->eps_p=constitutive_update->eps_p();
  LQdash->f_zero=std::make_shared<Constant>(0);

  // ADDED BY Q for assemble mean P_Dash and Q_Dash
  std::vector<double> Pdash_global;
  Pdash_global.push_back(Pdash0);
  auto PdashForm = std::make_shared<Plas3DWithForce::Form_PdashForm>(mesh);
  PdashForm->u2 = u;
  PdashForm->eps_p=constitutive_update->eps_p();
  PdashForm->f_zero=std::make_shared<Constant>(0);

  std::vector<double> Qdash_global;
  Qdash_global.push_back(0.0);
  auto QdashForm = std::make_shared<Plas3DWithForce::Form_QdashForm>(mesh);
  QdashForm->u2 = u;
  QdashForm->eps_p=constitutive_update->eps_p();
  QdashForm->f_zero=std::make_shared<Constant>(0);

  //---- New Posprocessing Added by Sam---------------------------------
  std::vector<double> eps_v_global;
  eps_v_global.push_back(0.0);
  auto eps_v = std::make_shared<Function>(V6);
  auto aEps_vForm=std::make_shared<Plas3DWithForce::Form_aEps_v>(V6,V6);
  auto LEps_vForm=std::make_shared<Plas3DWithForce::Form_LEps_v>(V6);
  LEps_vForm->u2=u;
  LEps_vForm->f_zero=ConstZeroScalar;
  auto Eps_vForm = std::make_shared<Plas3DWithForce::Form_Eps_vForm>(mesh);
  Eps_vForm->u2= u;
  Eps_vForm->f_zero=ConstZeroScalar;

  std::vector<double> epsP_v_global;
  epsP_v_global.push_back(0.0);
  auto epsP_v = std::make_shared<Function>(V6);
  auto aEpsP_vForm=std::make_shared<Plas3DWithForce::Form_aEpsP_v>(V6, V6);
  auto LEpsP_vForm=std::make_shared<Plas3DWithForce::Form_LEpsP_v>(V6);
  LEpsP_vForm->eps_p=constitutive_update->eps_p();
  LEpsP_vForm->f_zero=ConstZeroScalar;
  auto EpsP_vForm = std::make_shared<Plas3DWithForce::Form_EpsP_vForm>(mesh);
  EpsP_vForm->eps_p=constitutive_update->eps_p();
  EpsP_vForm->f_zero=ConstZeroScalar;

  std::vector<double> eps_zz_global;
  eps_zz_global.push_back(0.0);
  auto eps_zz = std::make_shared<Function>(V6);
  auto aEps_zz_Form = std::make_shared<Plas3DWithForce::Form_aEps_zz>(V6, V6);
  auto LEps_zz_Form = std::make_shared<Plas3DWithForce::Form_LEps_zz>(V6);
  LEps_zz_Form->u2=u;
  LEps_zz_Form->f_zero=ConstZeroScalar;
  auto Eps_zz = std::make_shared<Plas3DWithForce::Form_Eps_zzForm>(mesh);
  Eps_zz->u2=u;
  Eps_zz->f_zero=ConstZeroScalar;


  std::vector<double> eps_dev_global;
  eps_dev_global.push_back(0.0);
  auto eps_dev = std::make_shared<Function>(V6);
  auto aEps_devForm=std::make_shared<Plas3DWithForce::Form_aEps_dev>(V6, V6);
  auto LEps_devForm=std::make_shared<Plas3DWithForce::Form_LEps_dev>(V6);
  LEps_devForm->u2=u;
  LEps_devForm->f_zero=ConstZeroScalar;
  auto Eps_devForm = std::make_shared<Plas3DWithForce::Form_Eps_devForm>(mesh);
  Eps_devForm->u2=u;
  Eps_devForm->f_zero=ConstZeroScalar;

  std::vector<double> epsP_dev_global;
  epsP_dev_global.push_back(0.0);
  auto epsP_dev = std::make_shared<Function>(V6);
  auto aEpsP_devForm=std::make_shared<Plas3DWithForce::Form_aEpsP_dev>(V6, V6);
  auto LEpsP_devForm=std::make_shared<Plas3DWithForce::Form_LEpsP_dev>(V6);
  LEpsP_devForm->eps_p=constitutive_update->eps_p();
  LEpsP_devForm->f_zero=std::make_shared<Constant>(0);
  auto EpsP_devForm = std::make_shared<Plas3DWithForce::Form_EpsP_devForm>(mesh);
  EpsP_devForm->eps_p=constitutive_update->eps_p();
  EpsP_devForm->f_zero=ConstZeroScalar;

  //--------------------------------------------------------------------

  auto DomainVolume = std::make_shared<Plas3DWithForce::Form_DomainVolume>(mesh);
  DomainVolume->f_dom=std::make_shared<Constant>(1.0);
  double VOLUME = dolfin::assemble(*DomainVolume);

  auto nonlinear_problem = std::make_shared<fenicssolid::PlasticityProblem>(a, L, u, tangent, stress, bcs, SncF);

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver nonlinear_solver;
  nonlinear_solver.parameters["convergence_criterion"] = "incremental";
  nonlinear_solver.parameters["maximum_iterations"]    = 100;
  nonlinear_solver.parameters["relative_tolerance"]    = 1.0e-6;
  nonlinear_solver.parameters["absolute_tolerance"]    = 1.0e-12;

  // File names for output
  File file1("outputWithForce/disp.pvd");
  File file2("outputWithForce/eq_plas_strain.pvd");
  File file3("outputWithForce/eps.pvd");
  File file4("outputWithForce/eps_p.pvd");
  File file5("outputWithForce/stress.pvd");

  File file8("outputWithForce/eps_v.pvd");
  File file12("outputWithForce/eps_zz.pvd");
  File file9("outputWithForce/epsP_v.pvd");
  File file10("outputWithForce/eps_dev.pvd");
  File file11("outputWithForce/epsP_dev.pvd");

  File file6("outputWithForce/pdash.pvd");
  File file7("outputWithForce/qdash.pvd");

  // Equivalent plastic strain for visualisation
  auto eps_eq = std::make_shared<MeshFunction<double>>(mesh, mesh->topology().dim());

  while (step < steps)
  {
    // Use elastic tangent for first time step
    if (step == 0)
      t += Tdt;
    else if(step<startSuffusionAtStep || step>stopSuffusionAtStep)
    {
        t += Tdt0;
        SncF->phi_er_dot=0.0;
    }
    else
    {
        SncF->phi_er_dot=Phi_er_dot;
    }


    step++;
    std::cout << "step begin: " << step << std::endl;
    std::cout << "time: " << t << std::endl;

    // Solve non-linear problem
    nonlinear_solver.solve(*nonlinear_problem, *u->vector());

    // Update variables
    constitutive_update->update_history();
    //ADDED BY SAM
    //dolfin::solve(*aStrss==*LStrss, *strss);
    dolfin::solve(*aEps==*LEps, *Eps);
    dolfin::solve(*aEps_p==*LEps_p, *Eps_p);
    dolfin::solve(*aStrss==*LStress, *strss2);

    dolfin::solve(*aEps_vForm==*LEps_vForm, *eps_v);
    auto Eps_vFormTemp =dolfin::assemble(*Eps_vForm)/VOLUME;
    eps_v_global.push_back(Eps_vFormTemp);

    dolfin::solve(*aEpsP_vForm==*LEpsP_vForm, *epsP_v);
    auto EpsP_vFormTemp =dolfin::assemble(*EpsP_vForm)/VOLUME;
    epsP_v_global.push_back(EpsP_vFormTemp);

    dolfin::solve(*aEps_zz_Form==*LEps_zz_Form, *eps_zz);
    auto Eps_zzFormTemp = dolfin::assemble(*Eps_zz)/VOLUME;
    eps_zz_global.push_back(Eps_zzFormTemp);

    dolfin::solve(*aEps_devForm==*LEps_devForm, *eps_dev);
    auto Eps_devFormTemp =dolfin::assemble(*Eps_devForm)/VOLUME;
    eps_dev_global.push_back(Eps_devFormTemp);

    dolfin::solve(*aEpsP_devForm==*LEpsP_devForm, *epsP_dev);
    auto EpsP_devFormTemp =dolfin::assemble(*EpsP_devForm)/VOLUME;
    epsP_dev_global.push_back(EpsP_devFormTemp);

    dolfin::solve(*aPdash==*LPdash, *Pdash);
    dolfin::solve(*aQdash==*LQdash, *Qdash);
    auto PdashTemp = dolfin::assemble(*PdashForm)/VOLUME+Pdash0;
    auto QdashTemp = dolfin::assemble(*QdashForm)/VOLUME;
    Pdash_global.push_back(PdashTemp);
    Qdash_global.push_back(QdashTemp);
    // Write output to files
    file1 << *u;
    constitutive_update->eps_p_eq()->compute_mean(eps_eq);
    file2 << *eps_eq;
    //file3 << *strss;
    file3 << *Eps;
    file4 << *Eps_p;
    file5 << *strss2;

    file8<< *eps_v;
    file12<< *eps_zz;
    file9<< *epsP_v;
    file10<< *eps_dev;
    file11<< *epsP_dev;

    file6 << *Pdash;
    file7 << *Qdash;

    //file4.write(*Eps,t);
    //file5.write(*Eps_p,t);
    //file6.write(*strss2,t);
  }
  WriteToOctaveASCII("outputWithForce/eps_v.txt", eps_v_global);
  WriteToOctaveASCII("outputWithForce/eps_zz.txt", eps_zz_global);
  WriteToOctaveASCII("outputWithForce/epsP_v.txt", epsP_v_global);
  WriteToOctaveASCII("outputWithForce/eps_dev.txt", eps_dev_global);
  WriteToOctaveASCII("outputWithForce/epsP_dev.txt", epsP_dev_global);
  WriteToOctaveASCII("outputWithForce/pdash.txt", Pdash_global);
  WriteToOctaveASCII("outputWithForce/qdash.txt", Qdash_global);
  cout << "Solution norm: " << u->vector()->norm("l2") << endl;

  timer.stop();
  dolfin::list_timings(dolfin::TimingClear::clear, {dolfin::TimingType::wall});

  return 0;
}
