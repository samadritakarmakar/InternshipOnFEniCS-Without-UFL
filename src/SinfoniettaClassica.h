#ifndef SINFONIETTACLASSICA_H
#define SINFONIETTACLASSICA_H
#include <Eigen/Dense>
#include "PlasticityModel.h"


namespace fenicssolid
{
class SinfoniettaClassica: public PlasticityModel
{
public:

    SinfoniettaClassica();
    SinfoniettaClassica(double E, double nu, double beta, double phiDegree, double betaP,
                            double varKappa, double Pc, double varP);
    /// Delete copy constructor and assignement
    SinfoniettaClassica& operator=(const SinfoniettaClassica&) = delete;  // Disallow copying
    SinfoniettaClassica(const SinfoniettaClassica&) = delete;

    double hardening_parameter(double q) const;

    double f(const Eigen::Matrix<double, 6, 1>& stress,
             double q) const;

    void df(Eigen::Matrix<double, 6, 1>& df_dsigma,
                    const Eigen::Matrix<double, 6, 1>& stress) const;

    void dg(Eigen::Matrix<double, 6, 1>& dg_dsigma,
                  const Eigen::Matrix<double, 6, 1>& stress) const;

    void ddg(Eigen::Matrix<double, 6, 6>& ddg_ddsigma,
                     const Eigen::Matrix<double, 6, 1>& stress) const;

    void df_dq(double &df_dQ,
                       const double &q) const;

    void M(double &m, const Eigen::Matrix<double, 6, 1> &stress,
                        double q) const;

    void dM_dsigma(Eigen::Matrix<double, 6, 1>& dM_dsgma,
    const Eigen::Matrix<double, 6, 1>& stress,double q) const;

    double q_0() const;

    void set_q_0(double q0);

    Eigen::Matrix<double, 6, 6> random;

    Eigen::Matrix<double, 6, 1> GetVoigt_Delta_ij_over_Delta_mn() const;
    Eigen::Matrix<double, 6, 6> Get_P4_D() const;
    Eigen::Matrix<double, 3, 3> VoigtTo3By3Tensor(const Eigen::Matrix<double, 6, 1>& VoigtQuantity) const;
    Eigen::Matrix<double, 6, 1> Tensor3by3ToVoigt(const Eigen::Matrix<double, 3, 3>& TensorQuantity) const;
    Eigen::Matrix<double, 3, 3> GetDeviatoricQuantity(const Eigen::Matrix<double, 3, 3>& Quantity) const;
    Eigen::Matrix<double, 3, 3> VoigtTo3by3DeviatoricTensor(const Eigen::Matrix<double, 6, 1>& VoigtQuantity) const;



protected:
    double q_0_default=0.0;
private:
    double Z;
    double mu;
    double E_Internal;
    double nu_Internal;
    double hardeningParameter;
    double beta_Internal;
    double betaP_Internal;
    double varKappa_Internal;
    double Pc_0;
    double varP_Internal;
};

}
#endif // SINFONIETTACLASSICA_H
