#include <cmath>

#include "kariba/Electrons.hpp"
#include "kariba/Particles.hpp"
#include "kariba/constants.hpp"

//************************************************************************************************************

namespace kariba {

double multiplicity(double pspec) {
    // n tilde from Kelner et al. 2006. It is the number of produced pions for a
    // given proton distribution:
    double ntilde;
    if (pspec <= 2.25) {
        ntilde = .77;
    } else if (pspec >= 2.75) {
        ntilde = .67;
    } else {
        ntilde = .62;
    }
    return ntilde;
}

double prob() {

    int N = 20, i;                  // The steps of integration.
    double r = .573;                // r = 1-λ = m_μ^2/m_p^2 = 0.573.
    double xmin = 0., xmax = 1.;    // x = E_{particle}/E_p.
    double x, dx = (xmax - xmin) / (N - 1);
    double gn, hn1, hn2, fe;    // The energy distribution/probability f_e of electrons as a
                                // function of gn, hn1, hn2 (equations 36-39 from Kelner et al.
                                // 2006).
    double sum = 0.;
    double Bprob;

    for (i = 0; i < N; i++) {
        x = xmin + i * dx;
        gn = (3. - 2. * r) / (9. * (1. - r) * (1. - r)) *
             (9. * x * x - 6. * std::log(x) - 4. * x * x * x - 5.);
        hn1 = (3. - 2. * r) / (9. * (1. - r) * (1. - r)) *
              (9. * r * r - 6. * std::log(r) - 4. * r * r * r - 5.);
        hn2 = (1. + 2. * r) * (r - x) / 9. / r / r * (9. * (r + x) - 4. * (r * r + r * x + x * x));

        // The function f_e is given by equation 36 from Kelner et al. 2006 and
        // is: f_e(x) = g_v*H(x-r) + (h_v1(x) + h_v2(x))*H(r-x), with H(y) the
        // Heaviside function.
        if (x >= r) {    // H(x-r) = 1 and H(r-x) = 0.
            fe = gn;
        } else {    // H(x-r) = 0 and H(r-x) = 1.
            fe = hn1 + hn2;
        }
        sum = sum + fe * dx;
    }
    Bprob = 1. / sum;

    return Bprob;
}

double elec_dist_pp(double zen, double w) {

    double rmasses = .573;      // r = 1-λ = m_μ^2/m_p^2 = 0.573.The ratio of muon
                                // and proton energies.
    double yk;                  // x = E_e / E_pion  from Kelner et al. 2006.
    double gn, hn1, hn2, fe;    // The energy distribution/probability f_e of electrons as a
                                // function of gn, hn1, hn2 (eqs. 36-39 from  Kelner et al. 2006).
    yk = std::pow(10., (zen - w));    // x = E_e/E_pion from Kelner.
    gn = (3. - 2. * rmasses) / (9. * (1. - rmasses) * (1. - rmasses)) *
         (9. * yk * yk - 6. * std::log(yk) - 4. * yk * yk * yk -
          5.);    // Equation 37 from Kelner et al. 2006.
    hn1 = (3. - 2. * rmasses) / (9. * (1. - rmasses) * (1. - rmasses)) *
          (9. * rmasses * rmasses - 6. * std::log(rmasses) - 4. * rmasses * rmasses * rmasses -
           5.);    // Eq 38 from K06.
    hn2 = (1. + 2. * rmasses) * (rmasses - yk) / 9. / rmasses * rmasses *
          (9. * (rmasses + yk) -
           4. * (rmasses * rmasses + rmasses * yk + yk * yk));    // Eq 39 from K06.

    // The function f_e is given by equation 36 from Kelner et al. 2006 and is:
    // f_e(x) = g_v*H(x-r) + (h_v1(x) + h_v2(x))*H(r-x), with H(y) the Heaviside
    // function.
    if (yk >= rmasses) {    // H(x-r) = 1 and H(r-x) = 0.
        fe = gn;
    } else {    // H(x-r) = 0 and H(r-x) = 1.
        fe = hn1 + hn2;
    }

    return fe;
}

double elec_spec_pp(double Ep, double y) {

    double L = std::log(Ep);    // L = ln(Ep/1TeV) as definied in Kelner et al. 2006
                                // for the cross section
    double Betae, be, yke;      // The sub-functions that describe the function F_e(x,E_p) that
                                // implies the number of electrons in the interval (x,x+dx) per
                                // collision. In particular, eqs. 63-65 from Kelner et al. 2006:
    double Fespec;              // The spectrum of secondary electrons from pion decay.
                                // Equation 62 from Kelner et al. 2006

    // The sub-functions that describe the function F_e(x,E_p) that implies the
    // number of electrons in the interval (x,x+dx) per collision. In particular,
    // eq. 63 from Kelner et al. 2006.
    Betae = 1. / (69.5 + 2.65 * L + .3 * L * L);
    be = 1.0 / std::pow((.201 + .062 * L + .00041 * L * L), .25);    // Eq. 64 from K06.

    yke = (.279 + .141 * L + .0172 * L * L) / (.3 + (2.3 + L) * (2.3 + L));    // Eq. 65 from k06

    Fespec = Betae * (1. + yke * (y * std::log(10.)) * (y * std::log(10.))) *
             (1. + yke * (y * std::log(10.)) * (y * std::log(10.))) *
             (1. + yke * (y * std::log(10.)) * (y * std::log(10.))) /
             (std::pow(10., y) * (1. + .3 / std::pow(10., (y * be)))) * (-y * std::log(10.)) *
             (-y * std::log(10.)) * (-y * std::log(10.)) * (-y * std::log(10.)) *
             (-y * std::log(10.));
    return Fespec;
}

//***********************************************************************************************************

double production_rate(double ge, double x) {    // from Coppi & Blandford 1990
    double y = 2. * ge * x;
    double Heaviside;
    double R_gg;

    if (y - 1.0 >= 0.0) {
        Heaviside = 1.0;
    } else {
        Heaviside = 1.0e-200;
    }
    R_gg = constants::sigtom * constants::cee * 0.652 * (y * y - 1.0) / (y * y * y) * std::log(y) *
           Heaviside;
    return R_gg;
}

//***********************************************************************************************************

}    // namespace kariba
