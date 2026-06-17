#include <cmath>
#include <iostream>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include "kariba/Bknpower.hpp"
#include "kariba/Particles.hpp"
#include "kariba/constants.hpp"

namespace kariba {

Bknpower::Bknpower(size_t size) : Particles(size) {
    norm = 1.;

    mass_gr = constants::emgm;
    mass_kev = constants::emgm * constants::gr_to_kev;
}

//! Methods to set momentum/energy arrays
void Bknpower::set_p(double min, double brk, double ucom, double bfield, double betaeff, double r,
                     double fsc) {
    pmin = min;
    pbrk = brk;
    pmax = max_p(ucom, bfield, betaeff, r, fsc);
    const double p_grid_max = 10. * pmax;

    double pinc = (std::log10(p_grid_max) - std::log10(pmin)) / static_cast<double>(p.size() - 1);

    for (size_t i = 0; i < p.size(); i++) {
        p[i] = std::pow(10., std::log10(pmin) + static_cast<double>(i) * pinc);
        gamma[i] = std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    }
}

void Bknpower::set_p(double min, double brk, double gmax) {
    pmin = min;
    pbrk = brk;
    pmax = std::pow(std::pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;
    const double p_grid_max = 10 * pmax; 

    double pinc = (std::log10(p_grid_max) - std::log10(pmin)) / static_cast<double>(p.size() - 1);

    for (size_t i = 0; i < p.size(); i++) {
        p[i] = std::pow(10., std::log10(pmin) + static_cast<double>(i) * pinc);
        gamma[i] = std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    }
}

//! Method to set differential electron number density from known pspec,
//! normalization, and momentum array
void Bknpower::set_ndens() {
    for (size_t i = 0; i < p.size(); i++) {
        const double x = p[i] / pmax;
        const double C = Particles::cutoff_factor(x, cutoff_type);
        ndens[i] = norm * std::pow(p[i] / pbrk, -pspec1) /
                   (1. + std::pow(p[i] / pbrk, -pspec1 + pspec2)) * C;
    }
    initialize_gdens();
    gdens_differentiate();
}

//! methods to set the slopes, break and normalization
void Bknpower::set_pspec1(double s1) { pspec1 = s1; }

void Bknpower::set_pspec2(double s2) { pspec2 = s2; }

void Bknpower::set_brk(double brk) { pbrk = brk; }

//! Methods to calculate the normalization of the function
double norm_bkn_int(double x, void* pars) {
    BknParams* params = static_cast<BknParams*>(pars);

    double s1 = params->s1;
    double s2 = params->s2;
    double brk = params->brk;
    double max = params->max;
    double m = params->m;
    const int cutoff_type=params->cutoff_type; 

    double mom_int = std::pow(std::pow(x, 2.) - 1., 1. / 2.) * m * constants::cee;
    const double C = Particles::cutoff_factor(mom_int / max, cutoff_type);
    return std::pow(mom_int / brk, -s1) / (1. + std::pow(mom_int / brk, -s1 + s2)) *C;
}

void Bknpower::set_norm(double n) {
    double norm_integral, error, min, max;
    const double p_grid_max = 10. * pmax;

    min = std::pow(std::pow(pmin / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    max = std::pow(std::pow(p_grid_max / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);

    gsl_integration_workspace* w1;
    w1 = gsl_integration_workspace_alloc(100);
    gsl_function F1;
    auto params = BknParams{pspec1, pspec2, pbrk, pmax, mass_gr, cutoff_type};
    F1.function = &norm_bkn_int;
    F1.params = &params;
    gsl_integration_qag(&F1, min, max, 0, 1e-7, 100, 1, w1, &norm_integral, &error);
    gsl_integration_workspace_free(w1);

    norm = n / (norm_integral * mass_gr * constants::cee);
}

//! Injection function to be integrated in cooling
double injection_bkn_int(double x, void* pars) {
    InjectionBknParams* params = static_cast<InjectionBknParams*>(pars);
    double s1 = params->s1;
    double s2 = params->s2;
    double brk = params->brk;
    double max = params->max;
    double m = params->m;
    double n = params->n; 
    const int cutoff_type=params->cutoff_type; 

    double mom_int = std::pow(std::pow(x, 2.) - 1., 1. / 2.) * m * constants::cee;
    const double C = Particles::cutoff_factor(mom_int / max, cutoff_type);

    return n * std::pow(mom_int / brk, -s1) / (1. + std::pow(mom_int / brk, -s1 + s2)) *C;
}

//! Method to solve steady state continuity equation. NOTE: KN cross section not
//! included in IC cooling
void Bknpower::cooling_steadystate(double ucom, double n0, double bfield, double r,
                                   double betaeff) {
    double Urad = std::pow(bfield, 2.) / (8. * constants::pi) + ucom;
    double pdot_ad = betaeff * constants::cee / r;
    double pdot_rad = (4. * constants::sigtom * constants::cee * Urad) /
                      (3. * mass_gr * std::pow(constants::cee, 2.));
    double tinj = r / (constants::cee);

    double integral, error;
    gsl_function F1;
    auto params = InjectionBknParams{pspec1, pspec2, pbrk, pmax, mass_gr, n0, cutoff_type};
    F1.function = &injection_bkn_int;
    F1.params = &params;

    for (size_t i = 0; i < p.size(); i++) {
        if (i < p.size() - 1) {
            gsl_integration_workspace* w1;
            w1 = gsl_integration_workspace_alloc(100);
            gsl_integration_qag(&F1, gamma[i], gamma[i + 1], 1e1, 1e1, 100, 1, w1, &integral,
                                &error);
            gsl_integration_workspace_free(w1);

            ndens[i] =
                (integral / tinj) / (pdot_ad * p[i] / (mass_gr * constants::cee) +
                                     pdot_rad * (gamma[i] * p[i] / (mass_gr * constants::cee)));
        } else {
            ndens[p.size() - 1] =
                ndens[p.size() - 2] * std::pow(p[p.size() - 1] / p[p.size() - 2], -pspec2 - 1);
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
    gdens_differentiate();
}

//! Method to calculate maximum momentum of non thermal particles based on
//! acceleration and cooling timescales The estimate is identical to the old
//! agnjet but in momentum space; see Lucchini et al. 2019 for the math of the
//! old version
double Bknpower::max_p(double ucom, double bfield, double betaeff, double r, double fsc) {
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

//! simple method to check quantities.
void Bknpower::test() {
    std::cout << "Broken power-law distribution;" << std::endl;
    std::cout << "pspec1: " << pspec1 << std::endl;
    std::cout << "pspec2: " << pspec2 << std::endl;
    std::cout << "pbreak: " << pbrk << std::endl;
    std::cout << "Array size: " << p.size() << std::endl;
    std::cout << "Default normalization: " << norm << std::endl;
    std::cout << "Particle mass: " << mass_gr << std::endl;
}

}    // namespace kariba
