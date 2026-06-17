#include <cmath>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "kariba/Thermal.hpp"
#include "kariba/constants.hpp"

namespace kariba {

//! Class constructor to initialize object
Thermal::Thermal(size_t size) : Particles(size) {
    thnorm = 1.;

    mass_gr = constants::emgm;
    mass_kev = constants::emgm * constants::gr_to_kev;
}

//! Method to initialize momentum array a with default interval
void Thermal::set_p() {
    double emin = (1. / 100.) * Temp;    // minimum energy in kev, 1/100 lower than peak
    double emax = 20. * Temp;            // maximum energy in kev, 20 higher than peak
    double gmin, gmax, pmin, pmax, pinc;

    gmin = emin / mass_kev + 1.;
    gmax = emax / mass_kev + 1.;

    pmin = std::pow(std::pow(gmin, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;
    pmax = std::pow(std::pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;
    pinc = (std::log10(pmax) - std::log10(pmin)) / static_cast<double>(p.size() - 1);

    for (size_t i = 0; i < p.size(); i++) {
        p[i] = std::pow(10., std::log10(pmin) + static_cast<double>(i) * pinc);
        gamma[i] = std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    }
}

//! Method to set differential electron number density from known temperature,
//! normalization, and momentum array
void Thermal::set_ndens() {
    for (size_t i = 0; i < p.size(); i++) {
        ndens[i] = thnorm * std::pow(p[i], 2.) * std::exp(-gamma[i] / theta);
    }
    initialize_gdens();
    differentiate();
}

//! methods to set the temperature and normalization. NOTE: temperature must be
//! in ergs, no factor kb
void Thermal::set_temp_kev(double T) {
    Temp = T;
    theta = (Temp * constants::kboltz_kev2erg) / (mass_gr * constants::cee * constants::cee);
}

void Thermal::set_norm(double n) {
    thnorm = n / (std::pow(mass_gr * constants::cee, 3.) * theta * K2(1. / theta));
}

//! Evaluate Bessel function as in old agnjet
double Thermal::K2(double x) {
    double res;

    if (x < 0.1) {
        res = 2. / x / x;
    } else {
        res = gsl_sf_bessel_Kn(2, x);
    }

    return res;
}

//! simple method to check quantities.
void Thermal::test() {
    std::cout << "Thermal distribution;\n";
    std::cout << "Temperature: " << Temp << " erg, " << Temp / constants::kboltz_kev2erg
              << " kev\n";
    ;
    std::cout << "Array size: " << p.size() << "\n";
    std::cout << "Normalization: " << thnorm << "\n";
    std::cout << "Particle mass in grams: " << mass_gr << "\n";
    std::cout << "Particle mass in keV: " << mass_kev << "\n";
    std::cout << "kT/mc^2: " << theta << "\n\n";
}

}    // namespace kariba
