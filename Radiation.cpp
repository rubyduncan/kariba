#include <cmath>
#include <iostream>

#include "kariba/Radiation.hpp"
#include "kariba/constants.hpp"

namespace kariba {

Radiation::Radiation(size_t size)
    : en_phot(size, 0.0), num_phot(size, 0.0), en_phot_obs(size, 0.0), num_phot_obs(size, 0.0) {}

//! Methods to set viewing angle, beaming and geometry of emission region
void Radiation::set_beaming(double theta, double speed, double doppler) {
    angle = theta * constants::pi / 180.;
    beta = speed;
    dopfac = doppler;
}

void Radiation::set_inclination(double theta) { angle = theta * constants::pi / 180.; }

void Radiation::set_geometry(const std::string& geom, double l1, double l2) {
    if (geom == "cylinder") {
        r = l1;
        z = l2;
        vol = constants::pi * std::pow(r, 2.) * z;
        dopnum = 2.;
        geometry = geom;
    } else if (geom == "sphere") {
        r = l1;
        z = r;
        vol = (4. / 3.) * constants::pi * std::pow(r, 3.);
        dopnum = 3.;
        geometry = geom;
    } else {
        std::cout << "Input wrong, assuming sphere of radius 1 cm" << std::endl;
        std::cout << "Choose either sphere or cylinder!" << std::endl;
        r = 1.;
        z = 1.;
        vol = 4. * constants::pi / 3.;
        dopnum = 3.;
        geometry = "wrong!";
    }
}

void Radiation::set_geometry(const std::string& geom, double l1) {
    if (geom == "cylinder") {
        std::cout << "Only one length input, assuming the radius and height of "
                     "the cylinder are the same"
                  << std::endl;
        r = l1;
        z = l1;
        vol = constants::pi * std::pow(r, 2.) * z;
        dopnum = 2.;
        geometry = geom;
    } else if (geom == "sphere") {
        r = l1;
        z = r;
        vol = (4. / 3.) * constants::pi * std::pow(r, 3.);
        dopnum = 3.;
        geometry = geom;
    } else {
        std::cout << "Input wrong, assuming sphere of radius 1 cm" << std::endl;
        std::cout << "Choose either sphere or cylinder!" << std::endl;
        r = 1.;
        z = 1.;
        vol = 4. * constants::pi / 3.;
        dopnum = 3.;
        geometry = "wrong!";
    }
}

//! Simple integration method to integrate num_phot_obs and get the luminosity
//! between numin and numax
double Radiation::integrated_luminosity(double numin, double numax) {
    double temp = 0.;
    for (size_t i = 0; i < en_phot_obs.size() - 1; i++) {
        if (en_phot_obs[i] / constants::herg > numin &&
            en_phot_obs[i + 1] / constants::herg < numax) {
            temp = temp +
                   (1. / 2.) *
                       (en_phot_obs[i + 1] / constants::herg - en_phot_obs[i] / constants::herg) *
                       (num_phot_obs[i + 1] + num_phot_obs[i]);
        }
    }
    return temp;
}

//! Method to include a counterjet in cyclosycnchrotron/Compton classes
void Radiation::set_counterjet(bool flag) { counterjet = flag; }

//! Simple method to check arrays; only meant for debugging
void Radiation::test_arrays() {
    for (size_t i = 0; i < en_phot.size(); i++) {
        std::cout << en_phot[i] << " " << num_phot[i] << " "
                  << num_phot[i] * en_phot[i] * constants::herg << std::endl;
    }
}

}    // namespace kariba
