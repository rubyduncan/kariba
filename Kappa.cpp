#include <cmath>
#include <iostream>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include "kariba/Kappa.hpp"
#include "kariba/Particles.hpp"
#include "kariba/constants.hpp"

namespace kariba {

//! Class constructor to initialize object
Kappa::Kappa(size_t size) : Particles(size) {

    knorm = 1.;

    mass_gr = constants::emgm;
    mass_kev = constants::emgm * constants::gr_to_kev;
}

//! Method to set the temperature, using ergs as input
void Kappa::set_temp_kev(double T) {
    theta = T * constants::kboltz_kev2erg / (mass_gr * constants::cee * constants::cee);

    double emin = (1. / 100.) * T;    //!< minimum energy in kev, 1/100 lower than peak
    double emax = 20. * T;            //!< maximum energy in kev, 20 higher than peak
    double gmin, gmax;

    gmin = emin / mass_kev + 1.;
    gmax = emax / mass_kev + 1.;
    pmin = std::pow(std::pow(gmin, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;
    pmax = std::pow(std::pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;
}

void Kappa::set_kappa(double k) { kappa = k; }

//! Methods to set momentum/energy arrays and number density arrays
void Kappa::set_p(double ucom, double bfield, double betaeff, double r, double fsc) {
    pmax = std::max(max_p(ucom, bfield, betaeff, r, fsc), pmax);

    double pinc = (std::log10(pmax) - std::log10(pmin)) / static_cast<double>(p.size() - 1);

    for (size_t i = 0; i < p.size(); i++) {
        p[i] = std::pow(10., std::log10(pmin) + static_cast<double>(i) * pinc);
        gamma[i] = std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    }
}

//! Same as above, but assuming a fixed maximum Lorentz factor
void Kappa::set_p(double gmax) {
    pmax = std::pow(std::pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;

    double pinc = (std::log10(pmax) - std::log10(pmin)) / static_cast<double>(p.size() - 1);

    for (size_t i = 0; i < p.size(); i++) {
        p[i] = std::pow(10., std::log10(pmin) + static_cast<double>(i) * pinc);
        gamma[i] = std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    }
}

void Kappa::set_ndens() {
    for (size_t i = 0; i < gdens.size(); i++) {
        gdens[i] = knorm * gamma[i] * std::pow(std::pow(gamma[i], 2.) - 1., 1. / 2.) *
                   std::pow(1. + (gamma[i] - 1.) / (kappa * theta), -kappa - 1.);
    }
    initialize_pdens();
    differentiate();
}

//! Methods to calculate the normalization of the function
double norm_kappa_int(double x, void* pars) {
    KParams* params = static_cast<KParams*>(pars);
    double t = params->t;
    double k = params->k;

    return x * std::pow(std::pow(x, 2.) - 1., 1. / 2.) * std::pow(1. + (x - 1.) / (k * t), -k - 1.);
}

void Kappa::set_norm(double n) {
    double norm_integral, error, min, max;

    min = std::pow(std::pow(pmin / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    max = std::pow(std::pow(pmax / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);

    gsl_function F1;
    auto params = KParams{theta, kappa};
    gsl_integration_workspace* w1;
    w1 = gsl_integration_workspace_alloc(100);
    F1.function = &norm_kappa_int;
    F1.params = &params;
    gsl_integration_qag(&F1, min, max, 0, 1e-7, 100, 1, w1, &norm_integral, &error);
    gsl_integration_workspace_free(w1);

    knorm = n / norm_integral;
}

//! Method to solve steady state continuity equation. NOTE: KN cross section not
//! included in IC cooling
double injection_kappa_int(double x, void* pars) {
    InjectionKappaParams* params = static_cast<InjectionKappaParams*>(pars);
    double t = params->t;
    double k = params->k;
    double n = params->n;
    double m = params->m;

    double mom = std::pow(std::pow(x, 2.) - 1., 1. / 2.) * m * constants::cee;
    double diff = mom / (std::pow(m * constants::cee, 2.) *
                         std::pow(std::pow(mom / (m * constants::cee), 2.) + 1., 1. / 2.));

    return diff * n * x * std::pow(std::pow(x, 2.) - 1., 1. / 2.) *
           std::pow(1. + (x - 1.) / (k * t), -k - 1.);
}

void Kappa::cooling_steadystate(double ucom, double n0, double bfield, double r, double betaeff) {
    double Urad = std::pow(bfield, 2.) / (8. * constants::pi) + ucom;
    double pdot_ad = betaeff * constants::cee / r;
    double pdot_rad = (4. * constants::sigtom * constants::cee * Urad) /
                      (3. * mass_gr * std::pow(constants::cee, 2.));
    double tinj = r / (constants::cee);

    double integral, error;
    gsl_function F1;
    auto params = InjectionKappaParams{theta, kappa, knorm, mass_gr};
    F1.function = &injection_kappa_int;
    F1.params = &params;

    for (size_t i = 0; i < ndens.size(); i++) {
        if (i < ndens.size() - 1) {
            gsl_integration_workspace* w1;
            w1 = gsl_integration_workspace_alloc(100);
            gsl_integration_qag(&F1, gamma[i], gamma[i + 1], 1e1, 1e1, 100, 1, w1, &integral,
                                &error);
            gsl_integration_workspace_free(w1);

            ndens[i] =
                (integral / tinj) / (pdot_ad * p[i] / (mass_gr * constants::cee) +
                                     pdot_rad * (gamma[i] * p[i] / (mass_gr * constants::cee)));
        } else {
            ndens[ndens.size() - 1] =
                ndens[ndens.size() - 2] * std::pow(p[p.size() - 1] / p[p.size() - 2], -kappa);
        }
    }
    // the last bin is set by arbitrarily assuming cooled distribution; this is
    // necessary because the integral
    // above is undefined for the last bin

    // The last step requires a renormalization. The reason is that the result
    // of gsl_integration_qag strongly depends on the value of "size". Without
    // doing anything fancy, this can be fixed simply by ensuring that the total
    // integrated number of density equals n0 (which we know), and rescaling the
    // array ndens[i] by the appropriate constant.
    double renorm = count_particles() / n0;

    for (size_t i = 0; i < ndens.size(); i++) {
        ndens[i] = ndens[i] / renorm;
    }

    initialize_gdens();
    differentiate();
}

//! Method to calculate maximum momentum of non thermal particles based on
//! acceleration and cooling timescales
double Kappa::max_p(double ucom, double bfield, double betaeff, double r, double fsc) {
    double Urad, escom, accon, syncon, b, c, gmax;
    Urad = std::pow(bfield, 2.) / (8. * constants::pi) + ucom;
    escom = betaeff * constants::cee / r;
    syncon = (4. * constants::sigtom * Urad) / (3. * mass_gr * constants::cee);
    accon = (3. * fsc * constants::charg * bfield) / (4. * mass_gr * constants::cee);

    b = escom / syncon;
    c = accon / syncon;

    gmax = (-b + std::pow(std::pow(b, 2.) + 4. * c, 1. / 2.)) / 2.;

    return std::pow(std::pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;
}

void Kappa::test() {
    std::cout << "Kappa distribution;" << std::endl;
    std::cout << "kappa index: " << kappa << std::endl;
    std::cout << "dimensionless temperature: " << theta << std::endl;
    std::cout << "Array size: " << p.size() << std::endl;
    std::cout << "Default normalization: " << knorm << std::endl;
    std::cout << "Particle mass in grams: " << mass_gr << std::endl;
}

}    // namespace kariba
