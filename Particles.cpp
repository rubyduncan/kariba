#include <cmath>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "kariba/Particles.hpp"
#include "kariba/constants.hpp"

namespace kariba {

Particles::Particles(size_t size)
    : p(size, 0.0), ndens(size, 0.0), gamma(size, 0.0), gdens(size, 0.0), pdensp2_diff_logp(size, 0.0) {}

//! Simple numerical integrals /w trapeze method
double Particles::count_particles() {
    double temp = 0.0;
    for (size_t i = 0; i < ndens.size() - 1; i++) {
        temp = temp + (1. / 2.) * (p[i + 1] - p[i]) * (ndens[i + 1] + ndens[i]);
    }
    return temp;
}

double Particles::count_particles_energy() {
    double temp = 0.0;
    for (size_t i = 0; i < gamma.size() - 1; i++) {
        temp = temp + (1. / 2.) * (gamma[i + 1] - gamma[i]) * (gdens[i + 1] + gdens[i]);
    }
    return temp;
}

double Particles::av_p() {
    double temp = 0.;
    for (size_t i = 0; i < p.size() - 1; i++) {
        temp = temp + (1. / 2.) * (p[i + 1] - p[i]) * (p[i + 1] * ndens[i + 1] + p[i] * ndens[i]);
    }
    return temp / count_particles();
}

double Particles::av_gamma() {
    return std::pow(std::pow(av_p() / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
}

double Particles::av_psq() {
    double temp = 0.;
    for (size_t i = 0; i < p.size() - 1; i++) {
        temp = temp + (1. / 2.) * (p[i + 1] - p[i]) *
                          (std::pow(p[i + 1], 2.) * ndens[i + 1] + std::pow(p[i], 2.) * ndens[i]);
    }
    return temp / count_particles();
}

double Particles::av_gammasq() {
    return std::pow(av_psq() / std::pow(mass_gr * constants::cee, 2.) + 1., 1. / 2.);
}

//! Methods to set up energy space number density, as a function of momentum
//! space number density
void Particles::initialize_gdens() {
    for (size_t i = 0; i < gdens.size(); i++) {
        gdens[i] = ndens[i] * gamma[i] * mass_gr * constants::cee /
                   (std::pow(std::pow(gamma[i], 2.) - 1., 1. / 2.));
    }
}

//! Same as above but the other way around
void Particles::initialize_pdens() {
    for (size_t i = 0; i < gdens.size(); i++) {
        ndens[i] = gdens[i] * p[i] /
                   (std::pow(mass_gr * constants::cee, 2.) *
                    std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.));
    }
}

// void Particles::gdens_differentiate() {
//     std::vector<double> temp;
//     size_t size = gdens.size();

//     for (size_t i = 0; i < size; i++) {
//         temp.push_back(gdens[i] / (std::pow(gamma[i], 1.)));
//     }

//     for (size_t i = 0; i < size - 1; i++) {
//         gdens_diff[i] = (temp[i + 1] - temp[i]) /
//                         (mass_gr * std::pow(constants::cee, 2.) * (gamma[i + 1] - gamma[i]));
//     }

//     gdens_diff[size - 1] = gdens_diff[size - 2];
// }

void Particles::differentiate() {
    const size_t size = gdens.size();

    // Work in log space.  Guard against zero/negative gdens with a floor
    // (should not occur for physical distributions, but avoids log(0)).
    const double floor_val = 1e-300;

    auto safe_log = [&](double x) -> double {
        return std::log(std::max(x, floor_val));
    };

    for (size_t i = 0; i < size; i++) {
        double d_lnf_d_lnp; // f = dn/dp * p^-2, as in syn. absorption coefficient
        double power = 2.;

        if (i == 0) {
            // Forward difference
            d_lnf_d_lnp = safe_log(ndens[1]/ndens[0] * std::pow(p[0]/p[1], power)) /
                          safe_log(p[1]/p[0]);
        } else if (i == size - 1) {
            // Backward difference
            d_lnf_d_lnp = safe_log(ndens[size-1]/ndens[size-2] * std::pow(p[size-2]/p[size-1], power)) /
                          safe_log(p[size-1]/p[size-2]);
        } else {
            // Centered difference in log-log space
            d_lnf_d_lnp = safe_log(ndens[i+1]/ndens[i-1] * std::pow(p[i-1]/p[i+1], power)) /
                          safe_log(p[i+1]/p[i-1]);
        }

        pdensp2_diff_logp[i] = (ndens[i]/std::pow(p[i], power)) * d_lnf_d_lnp;
        // multiply later, this is in p = rho * m_e * c
        // pdensp2_diff_logp[i] *= std::pow(mass_gr * constants::cee, 3);
    }
}


void Particles::set_mass(double m) {
    mass_gr = m;
    mass_kev = m * constants::gr_to_kev;
}

//! simple method to check arrays; only meant for debugging
void Particles::test_arrays() {
    for (size_t i = 0; i < p.size(); i++) {
        std::cout << p[i] << " " << gamma[i] << " " << ndens[i] << " " << ndens[i] * p[i]
                  << std::endl;
    }
}

}    // namespace kariba
