#pragma once

#include <iostream>

#include <gsl/gsl_spline.h>

#include "Particles.hpp"

namespace kariba {

//! Class for non-thermal particles, inherited from the generic Particles class
//! in Particles.hpp note: ndens is number density per unit momentum
class Powerlaw : public Particles {
  protected:
    double pspec, plnorm;
    double pmin, pmax;
    bool isEfficient;    //!< Proton acceleration

  public:
    Powerlaw(size_t size);

    void set_p(double min, double ucom, double bfield, double betaeff, double r, double fsc);
    void set_p(double min, double gmax);
    void set_ndens();
    void set_pspec(double s1);
    void set_norm(double n);

    bool get_Efficiency() const { return isEfficient; }

    void cooling_steadystate(double ucom, double n0, double bfield, double r, double tshift);
    double max_p(double ucom, double bfield, double betaeff, double r, double fsc);

    // some functions for protons
    void set_energy(double gpmin, double fsc, double f_beta, double bfield, double r_g, double z,
                    double r, int infosw, double protdens, double nwind, double Uradjet,
                    const std::string& outputConfiguration, const std::string& source);
    void check_secondary_charged_syn(double bfield, double gpmax);
    void ProtonTimescales(double& logdgp, double fsc, double f_beta, double bfield, double gpmin,
                          double& gpmax, double r_g, double z, double r, int infosw, double nwind,
                          double Uradjet, const std::string& outputConfiguration,
                          const std::string& source);
    double sigma_pp(double Ep);
    double set_normprot(double nprot);
    void set_gdens(double& plfrac_p, double Up, double protdens);
    void set_gdens(double r, double protdens, double nwind, double bfield, double plfrac,
                   double Urad);
    void set_gdens_pdens(double r, double beta, double Ljet, double ep, double pspec,
                         double& protdens);

    // secondary electrons from pp
    void set_pp_elecs(gsl_interp_accel* acc_Jp, gsl_spline* spline_Jp, double ntot_prot,
                      double nwind, double plfrac, double gammap_min, double gammap_max,
                      double bfield, double r);
    // convert secondary electrons from pg from Neutrinos units proper units for
    // synchrotron radiation
    void set_pg_electrons(const std::vector<double>& energy, const std::vector<double>& density,
                          double f_beta, double r, double vol, double B);

    // function for electrons from γγ annihilation
    void Qggeefunction(double r, double vol, double bfield, size_t phot_number,
                       const std::vector<double>& en_perseg, const std::vector<double>& lum_perseg,
                       double gmax);

    void test();
};

}    // namespace kariba
