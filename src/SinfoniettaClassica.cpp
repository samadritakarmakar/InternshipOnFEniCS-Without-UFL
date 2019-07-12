#include <iostream>
#include <math.h>
#include "SinfoniettaClassica.h"
using namespace fenicssolid;

SinfoniettaClassica::SinfoniettaClassica(double E, double nu, double beta, double phiDegree, double betaP, double varKappa,
                                         double Pc=0, double varP=0, double Pdash0=0):PlasticityModel(E, nu),
        E_Internal(E), nu_Internal(nu), beta_Internal(beta), betaP_Internal(betaP), varKappa_Internal(varKappa), Pc_0(Pc),
        varP_Internal(varP), Pdash0_Internal(Pdash0)
    {
        hardeningParameter=1.0/betaP;
        q_0_default=-std::log(Pc_0);
        double phiRad=phiDegree/180.0*M_PI;
        Z=6.0*sin(phiRad)/(3.0-sin(phiRad));
        mu=(9.0-Z*Z)/(3.0-Z*Z+2.0*pow(Z,3.0)/9.0);
        //std::cout<<"Z="<<Z<<"\n";
    }

double SinfoniettaClassica::hardening_parameter(double q) const
{
    return hardeningParameter;
}

double SinfoniettaClassica::f(const Eigen::Matrix<double, 6, 1> &stress, double q) const
{
    Eigen::Matrix<double, 3, 3> stressTensor, onesT;
    stressTensor=VoigtTo3By3Tensor(stress);
    //std::cout<<"stress =\n"<<stress<<"\n";
    //std::cout<<"stressTensor =\n"<<stressTensor<<"\n";
    Eigen::Matrix<double, 3, 3> stressDeviatoric;
    stressDeviatoric=GetDeviatoricQuantity(stressTensor);
    double eps=std::numeric_limits<double>::epsilon();
    eps=std::sqrt(eps);
    double p_dash=Get_MeanStress(stressTensor);
    Eigen::Matrix<double, 3, 3> xi=stressDeviatoric/p_dash;
    double J_2xi=pow(xi.norm(),2.0);
    double J_3xi=3.0*xi.determinant();
    double f=3*beta_Internal*(mu-3.0)*(log(p_dash)+q)+9.0/4.0*(mu-1.0)*J_2xi+mu*J_3xi;
    //std::cout<<"f= "<<f <<"\n";
    if(isnan(f)&& p_dash<=0)
    {
        std::cout<<"p_dash="<<p_dash<<"\n";
    }
    return f;
}

void SinfoniettaClassica::df(Eigen::Matrix<double, 6, 1>& df_dsigma,
                  const Eigen::Matrix<double, 6, 1>& stress) const
{
    Eigen::Matrix<double, 3, 3> stressTensor;
    stressTensor=VoigtTo3By3Tensor(stress);
    Eigen::Matrix<double, 3, 3> stressDeviatoric;
    stressDeviatoric=GetDeviatoricQuantity(stressTensor);
    double p_dash=Get_MeanStress(stressTensor);
    Eigen::Matrix<double, 3, 3> df_dsigma_temp1=-(3.0*beta_Internal*(mu-3.0)/p_dash
                                -9.0/2.0*(mu-1.0)*pow(stressDeviatoric.norm(),2.0)/pow(p_dash,3.0)
                               -9.0*mu*stressDeviatoric.determinant()/pow(p_dash,4.0))*1.0/3.0*Eigen::Matrix<double, 3, 3>::Identity();

    Eigen::Matrix<double, 6, 6> P_D4=Get_P4_D();
    Eigen::Matrix<double, 3, 3> df_dsigma_temp2=(9.0/2.0*(mu-1)*stressDeviatoric/pow(p_dash,2.0)
                                               +3*mu*(stressDeviatoric*stressDeviatoric)/pow(p_dash,3.0));
    //Eigen::Matrix<double, 3, 3> df_dsigma_temp2=(9.0/2.0*(mu-1)*stressDeviatoric/pow(p_dash,2.0)
                                                  // +3*mu*(stressDeviatoric.determinant()*(stressDeviatoric.inverse().transpose()))/pow(p_dash,3.0));

    df_dsigma=Tensor3by3ToVoigt(df_dsigma_temp1)+(Tensor3by3ToVoigt(df_dsigma_temp2).transpose()*P_D4).transpose();
}

void SinfoniettaClassica::dg(Eigen::Matrix<double, 6, 1>& dg_dsigma,
                  const Eigen::Matrix<double, 6, 1>& stress) const
{
    Eigen::Matrix<double, 3, 3> stressTensor;
    stressTensor=VoigtTo3By3Tensor(stress);
    Eigen::Matrix<double, 3, 3> stressDeviatoric;
    stressDeviatoric=GetDeviatoricQuantity(stressTensor);
    double p_dash=Get_MeanStress(stressTensor);
    Eigen::Matrix<double, 3, 3> dg_dsigma_temp1=-(9.0*(mu-3.0)/p_dash
                                -9.0/2.0*(mu-1.0)*pow(stressDeviatoric.norm(),2)/pow(p_dash,3)
                               -9.0*mu*stressDeviatoric.determinant()/pow(p_dash,4))*1.0/3.0*Eigen::Matrix<double, 3, 3>::Identity();

    Eigen::Matrix<double, 6, 6> P_D4=Get_P4_D();
    Eigen::Matrix<double, 3, 3> dg_dsigma_temp2=(9.0/2.0*(mu-1.0)*stressDeviatoric/pow(p_dash,2.0)
                                               +3.0*mu*(stressDeviatoric*stressDeviatoric)/pow(p_dash,3.0));
    //Eigen::Matrix<double, 3, 3> dg_dsigma_temp2=(9.0/2.0*(mu-1.0)*stressDeviatoric/pow(p_dash,2.0)
                                                  // +3.0*mu*(stressDeviatoric.determinant()*(stressDeviatoric.inverse().transpose()))/pow(p_dash,3.0));
    //dg_dsigma=Tensor3by3ToVoigt(dg_dsigma_temp1)+P_D4*Tensor3by3ToVoigt(dg_dsigma_temp2);
    dg_dsigma=Tensor3by3ToVoigt(dg_dsigma_temp1)+
            (Tensor3by3ToVoigt(dg_dsigma_temp2).transpose()*P_D4).transpose();
}

void SinfoniettaClassica::ddg(Eigen::Matrix<double, 6, 6>& ddg_ddsigma,
                     const Eigen::Matrix<double, 6, 1>& stress) const
{
    Eigen::Matrix<double, 6, 1> dg_dsigma;
    Eigen::Matrix<double, 6, 1> dg_dsigma_h;
    double eps=std::numeric_limits<double>::epsilon();
    eps=std::sqrt(eps);
    //cout<<"eps = "<<eps<<"\n";
    Eigen::Matrix<double, 6, 1> onesVec;
    onesVec.setOnes();
    Eigen::Matrix<double, 6, 6> stress_eps=onesVec*stress.transpose()+eps*Eigen::Matrix<double, 6, 6>::Identity();
     dg(dg_dsigma, stress);
     //std::cout<<"dg_dsigma=\n"<<dg_dsigma<<"\n";
     //cout<<"stress_eps=\n"<<stress_eps<<"\n";
    for(int col=0; col<6; col++)
    {
        dg(dg_dsigma_h, stress_eps.row(col).transpose());
        //cout<<"dg_dsigma_h\n"<<dg_dsigma_h<<" for stress_eps"<<stress_eps.row(col)<<"at col "<<col<<"\n";
        ddg_ddsigma.col(col)=(dg_dsigma_h-dg_dsigma)/eps;
    }
}

void SinfoniettaClassica::df_dq(double &df_dQ,
                       const double &q) const
{
    df_dQ=3*beta_Internal*(mu-3.0);
}

void SinfoniettaClassica::M(double &m,
                        const Eigen::Matrix<double, 6, 1> &stress,
                        double q) const
{
    Eigen::Matrix<double, 6, 1> dg_dsigma;
    dg(dg_dsigma, stress);
    Eigen::Matrix<double, 6, 1> Voigt_Delta_ij_over_Delta_mn=GetVoigt_Delta_ij_over_Delta_mn();
    Eigen::Matrix<double, 6, 6> P4_D=Get_P4_D();
    //m=-dg_dsigma.transpose()*Voigt_Delta_ij_over_Delta_mn+varKappa_Internal*(dg_dsigma.transpose()*P4_D).norm();
            //+pow(varP_Internal,3.0)*sqrt(VoigtTo3By3Tensor((dg_dsigma.transpose()*P4_D).transpose()).determinant());
    Eigen::Matrix<double, 3, 3> DGDS = VoigtTo3By3Tensor(dg_dsigma);
    m=-DGDS.trace()+varKappa_Internal*sqrt((GetDeviatoricQuantity(DGDS)*GetDeviatoricQuantity(DGDS).transpose()).trace());
}

void SinfoniettaClassica::dM_dsigma(Eigen::Matrix<double, 6, 1>& dM_dsgma,
    const Eigen::Matrix<double, 6, 1>& stress, double q) const
{
    double m;
    M(m, stress, q);
    double eps=std::numeric_limits<double>::epsilon();
    eps=std::sqrt(eps);
    Eigen::Matrix<double, 6, 1> onesVec;
    onesVec.setOnes();
    Eigen::Matrix<double, 6, 6> stress_eps=onesVec*stress.transpose()+eps*Eigen::Matrix<double, 6, 6>::Identity();
    double m_h;
    //cout<<"stress_eps =\n"<<stress_eps<<"\n";
    for(int col=0; col<6; col++)
    {
        M(m_h, stress_eps.row(col), q);
        //cout<<"m_h = "<<m_h<<" at col "<<col<<"\n";
        dM_dsgma(col,0)=(m_h-m)/eps;
    }
}

double SinfoniettaClassica::q_0() const
{
    return q_0_default;
}

double SinfoniettaClassica::miscellaneous() const
{
    return phi_er_dot/betaP_Internal;
}

void SinfoniettaClassica::set_q_0(const double q0)
{
    q_0_default=q0;
}


Eigen::Matrix<double, 6, 6> SinfoniettaClassica::Get_stress_eps(double& eps, Eigen::Matrix<double, 6, 1> stress) const
{
    eps=std::numeric_limits<double>::epsilon();
    eps=std::sqrt(eps);
    Eigen::Matrix<double, 6, 1> onesVec;
    onesVec.setOnes();
    Eigen::Matrix<double, 6, 6> stress_eps=onesVec*stress.transpose()+eps*Eigen::Matrix<double, 6, 6>::Identity();
    return stress_eps;
}

double SinfoniettaClassica::Get_MeanStress(const Eigen::Matrix<double, 3, 3>& stressTensor) const
{
    return -1.0/3.0*stressTensor.trace()+Pdash0_Internal;
}

Eigen::Matrix<double, 6, 1> SinfoniettaClassica::GetVoigt_Delta_ij_over_Delta_mn() const
{
    Eigen::Matrix<double, 6, 1> Voigt_Delta_ij_over_Delta_mn;
        Voigt_Delta_ij_over_Delta_mn.setZero();
        //Voigt_Delta_ij_over_Delta_mn.block(0,0,2,0)=Eigen::Matrix<double, 3, 1>::Ones();
        Voigt_Delta_ij_over_Delta_mn(0,0)=1.0;
        Voigt_Delta_ij_over_Delta_mn(1,0)=1.0;
        Voigt_Delta_ij_over_Delta_mn(2,0)=1.0;
        return Voigt_Delta_ij_over_Delta_mn;
}

Eigen::Matrix<double, 6, 6> SinfoniettaClassica::Get_P4_D() const
{

    Eigen::Matrix<double, 6, 1> Voigt_Delta_ij_over_Delta_mn=GetVoigt_Delta_ij_over_Delta_mn();
    Eigen::Matrix<double, 6, 6> DevIdentity=Voigt_Delta_ij_over_Delta_mn*Voigt_Delta_ij_over_Delta_mn.transpose();
    Eigen::Matrix<double, 6, 6> P_D4=Eigen::Matrix<double, 6, 6>::Identity()-1.0/3.0*DevIdentity;
    return P_D4;
}

Eigen::Matrix<double, 3, 3> SinfoniettaClassica::VoigtTo3By3Tensor(const Eigen::Matrix<double, 6, 1>& VoigtQuantity) const
{
    Eigen::Matrix<double, 3, 3> TensorQuantity;
    TensorQuantity(0,0)=VoigtQuantity(0,0);
    TensorQuantity(1,1)=VoigtQuantity(1,0);
    TensorQuantity(2,2)=VoigtQuantity(2,0);
    TensorQuantity(0,1)=VoigtQuantity(3,0);
    TensorQuantity(1,0)=VoigtQuantity(3,0);
    TensorQuantity(0,2)=VoigtQuantity(4,0);
    TensorQuantity(2,0)=VoigtQuantity(4,0);
    TensorQuantity(1,2)=VoigtQuantity(5,0);
    TensorQuantity(2,1)=VoigtQuantity(5,0);
    return TensorQuantity;
}

Eigen::Matrix<double, 6, 1> SinfoniettaClassica::Tensor3by3ToVoigt(const Eigen::Matrix<double, 3, 3>& TensorQuantity) const
{
    Eigen::Matrix<double, 6, 1> VoigtQuantity;
    VoigtQuantity(0,0)=TensorQuantity(0,0);
    VoigtQuantity(1,0)=TensorQuantity(1,1);
    VoigtQuantity(2,0)=TensorQuantity(2,2);
    VoigtQuantity(3,0)=TensorQuantity(0,1);
    VoigtQuantity(4,0)=TensorQuantity(0,2);
    VoigtQuantity(5,0)=TensorQuantity(1,2);
    return VoigtQuantity;
}

Eigen::Matrix<double, 3, 3>  SinfoniettaClassica::GetDeviatoricQuantity(const Eigen::Matrix<double, 3, 3>& Quantity) const
{
    const Eigen::Matrix<double, 3, 3> DeviatoricQuantity=Quantity-1.0/3.0*Quantity.trace()*Eigen::Matrix<double, 3, 3>::Identity();
    return DeviatoricQuantity;
}

Eigen::Matrix<double, 3, 3> SinfoniettaClassica::VoigtTo3by3DeviatoricTensor(const Eigen::Matrix<double, 6, 1>& VoigtQuantity) const
{
    Eigen::Matrix<double, 3, 3> TensorQuantity;
    TensorQuantity=VoigtTo3By3Tensor(VoigtQuantity);
    Eigen::Matrix<double, 3, 3> stressDeviatoric;
    stressDeviatoric=GetDeviatoricQuantity(TensorQuantity);
    return stressDeviatoric;
}

