#include <cmath>
#include <iostream>

#include <gsl/gsl_integration.h>

#include "kariba/ShSDisk.hpp"
#include "kariba/constants.hpp"

namespace kariba {

//! Constructor for the disk. No input parameters because the sized of the arrays
//! is set automatically, and every other property needs to be handled by the
//! setters below
ShSDisk::ShSDisk(size_t size) : Radiation(size) {}    // size = 50 in the declaration

//! return SD spectrum over a given radius, frequency to be integrated over
//! radius
double disk_int(double lr, void* pars) {
    DiskObsParams* params = static_cast<DiskObsParams*>(pars);
    double tin = params->tin;
    double rin = params->rin;
    double nu = params->nu;
    double r, temp, fac, bb;

    r = std::exp(lr);
    temp = tin * std::pow(rin / r, 0.75);
    fac = constants::herg * nu / (constants::kboltz * temp);

    if (fac < 1.e-3) {
        bb = 2. * constants::herg * std::pow(nu, 3.) / (std::pow(constants::cee, 2.) * fac);
    } else {
        bb = 2. * constants::herg * std::pow(nu, 3.) /
             (std::pow(constants::cee, 2.) * (std::exp(fac) - 1.));
    }

    return 2. * constants::pi * std::pow(r, 2.) * bb;
}

void ShSDisk::disk_spectrum() {
    double result, error;

    for (size_t k = 0; k < en_phot_obs.size(); k++) {
        gsl_integration_workspace* w1;
        w1 = gsl_integration_workspace_alloc(100);
        gsl_function F1;
        auto F1params = DiskObsParams{Tin, r, en_phot_obs[k] / constants::herg};
        F1.function = &disk_int;
        F1.params = &F1params;
        gsl_integration_qag(&F1, std::log(r), std::log(z), 0, 1e-2, 100, 2, w1, &result, &error);
        gsl_integration_workspace_free(w1);

        num_phot[k] = result;
        num_phot_obs[k] = cos(angle) * result;
    }
}

//! this method removes a fraction f from the observed disk luminosity, assuming
//! that it was absorbed and reprocessed into some other unspecified radiative
//! mechanism
void ShSDisk::cover_disk(double f) {
    for (size_t k = 0; k < num_phot_obs.size(); k++) {
        num_phot_obs[k] = num_phot_obs[k] * (1. - f);
    }
}

//! Simple integration method to integrate num_phot_obs and get the total
//! luminosity of the disk
double ShSDisk::total_luminosity() {
    double temp = 0.;
    for (size_t i = 0; i < en_phot_obs.size() - 1; i++) {
        temp =
            temp + (1. / 2.) *
                       (en_phot_obs[i + 1] / constants::herg - en_phot_obs[i] / constants::herg) *
                       (num_phot_obs[i + 1] / cos(angle) + num_phot_obs[i] / cos(angle));
    }
    return temp;
}

void ShSDisk::set_mbh(double M) {
    Mbh = M;
    Rg = constants::gconst * Mbh * constants::msun / (constants::cee * constants::cee);
}

//! Note: for disk r is inner radius, z is the outer radius, both are input in cm
void ShSDisk::set_rin(double R) { r = R; }

void ShSDisk::set_rout(double R) { z = R; }

//! NOTE: the constant factor to go between inner temperature Tin and disk
//! lunminosity Ldisk is 2 because the model uses the zero torque inner boundary
//! condition, Kubota et al. 1998, hence the factor 2 rather than 4pi when
//! converting between luminosity and temperature
void ShSDisk::set_luminosity(double L) {
    double emin, emax, einc;

    Ldisk = L;
    Tin = std::pow(Ldisk * 1.25e38 * Mbh / (2. * constants::sbconst * std::pow(r, 2.)), 0.25);
    Hratio = std::max(0.1, Ldisk);
    emin = 0.0001 * constants::kboltz * Tin;
    emax = 30. * constants::kboltz * Tin;
    einc = (std::log10(emax) - std::log10(emin)) / static_cast<double>(en_phot.size() - 1);

    for (size_t i = 0; i < en_phot.size(); i++) {
        en_phot_obs[i] = std::pow(10., std::log10(emin) + static_cast<double>(i) * einc);
        en_phot[i] = en_phot_obs[i];
        num_phot[i] = 0.;
        num_phot_obs[i] = 0.;
    }
}

void ShSDisk::set_tin_kev(double T) {
    double emin, emax, einc;

    // note: 1 keV = constants::kboltz_kev2erg/constants::kboltz keV
    Tin = T * constants::kboltz_kev2erg / constants::kboltz;
    Ldisk = 2. * constants::sbconst * std::pow(Tin, 4.) * std::pow(r, 2.) / (1.25e38 * Mbh);
    Hratio = std::max(0.1, Ldisk);
    emin = 0.0001 * constants::kboltz * Tin;
    emax = 30. * constants::kboltz * Tin;
    einc = (std::log10(emax) - std::log10(emin)) / static_cast<double>(en_phot.size() - 1);

    for (size_t i = 0; i < en_phot.size(); i++) {
        en_phot_obs[i] = std::pow(10., std::log10(emin) + static_cast<double>(i) * einc);
        en_phot[i] = en_phot_obs[i];
        num_phot[i] = 0.;
        num_phot_obs[i] = 0.;
    }
}

void ShSDisk::set_tin_k(double T) {
    double emin, emax, einc;

    Tin = T;
    Ldisk = 2. * constants::sbconst * std::pow(Tin, 4.) * std::pow(r, 2.) / (1.25e38 * Mbh);
    Hratio = std::max(0.1, Ldisk);
    emin = 0.0001 * constants::kboltz * Tin;
    emax = 30. * constants::kboltz * Tin;
    einc = (std::log10(emax) - std::log10(emin)) / static_cast<double>(en_phot.size() - 1);

    for (size_t i = 0; i < en_phot.size(); i++) {
        en_phot_obs[i] = std::pow(10., std::log10(emin) + static_cast<double>(i) * einc);
        en_phot[i] = en_phot_obs[i];
        num_phot[i] = 0.;
        num_phot_obs[i] = 0.;
    }
}

void ShSDisk::test() {
    std::cout << "Inner disk radius: " << r << " cm, " << r / Rg << " rg; outer disk radius: " << z
              << " cm, " << z / Rg << " rg; disk scale height: " << Hratio << std::endl;
    std::cout << "Inner disk temperature: " << Tin * constants::kboltz / constants::kboltz_kev2erg
              << " kev; emitted  disk luminosity in Eddington units: " << Ldisk
              << " and erg s^-1: " << total_luminosity() << std::endl;
}

}    // namespace kariba
