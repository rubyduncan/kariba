/*************************************************************************************************************
Gamma-rays from neutral pion decay, products of inelastic pp and pγ collisions
*************************************************************************************************************/

#include <array>
#include <cmath>

#include <gsl/gsl_integration.h>

#include "kariba/GammaRays.hpp"
#include "kariba/constants.hpp"

namespace kariba {

static const size_t nitems = 22;
static const std::array<double, nitems> etagTable = {1.1, 1.2,  1.3,  1.4, 1.5,  1.6,  1.7, 1.8,
                                                     1.9, 2.0,  3.0,  4.0, 5.0,  6.0,  7.0, 8.0,
                                                     9.0, 10.0, 20.0, 30., 40.0, 100.0};
static const std::array<double, nitems> sgTable = {
    0.0768, 0.106, 0.182, 0.201, 0.219,  0.216,  0.233,  0.233, 0.248,  0.244,  0.188,
    0.131,  0.120, 0.107, 0.102, 0.0932, 0.0838, 0.0761, 0.107, 0.0928, 0.0722, 0.0479};
static const std::array<double, nitems> deltagTable = {
    0.544, 0.540, 0.750, 0.791, 0.788, 0.831, 0.839, 0.825, 0.805, 0.779, 1.23,
    1.82,  2.05,  2.19,  2.23,  2.29,  2.37,  2.43,  2.27,  2.33,  2.42,  2.59};
static const std::array<double, nitems> BetagTable = {
    2.86e-19, 2.24e-18, 5.61e-18, 1.02e-17, 1.60e-17, 2.23e-17, 3.10e-17, 4.07e-17,
    5.30e-17, 6.74e-17, 1.51e-16, 1.24e-16, 1.37e-16, 1.62e-16, 1.71e-16, 1.78e-16,
    1.84e-16, 1.93e-16, 4.74e-16, 7.70e-16, 1.06e-15, 2.73e-15};

Grays::Grays(size_t size, double numin, double numax) : Radiation(size) {
    en_phot_obs.resize(en_phot_obs.size() * 2, 0.0);
    num_phot_obs.resize(num_phot_obs.size() * 2, 0.0);

    size_t lsize = en_phot.size();
    double nuinc = (std::log10(numax) - std::log10(numin)) / static_cast<double>(lsize - 1);
    for (size_t i = 0; i < lsize; i++) {
        en_phot[i] =
            std::pow(10., std::log10(numin) + static_cast<double>(i) * nuinc) * constants::herg;
        en_phot_obs[i] = en_phot[i];
        en_phot_obs[i + lsize] = en_phot[i];
    }
}

//************************************************************************************************************
void Grays::set_grays_pp(double p, double gammap_min, double gammap_max, double ntot_prot,
                         double ntargets, double plfrac, gsl_interp_accel* acc_Jp,
                         gsl_spline* spline_Jp) {

    double ntilde =
        set_ntilde(p);    // The number of produced pions for a given proton distribution
    double pp_targets = target_protons(ntot_prot, ntargets, plfrac);
    double Epcode_max = gammap_max * constants::pmgm * constants::cee * constants::cee *
                        constants::erg * 1.0e-12;    // The proton energy in TeV

    int N = 60;                        // The number of steps of the secondary particle (e.g., pion)
                                       // energy.
    double Ep;                         // The energy of the proton in TeV
    double Eg;                         // The energy of the gamma-ray photon in TeV
    double Jp;                         // The number of non-thermal protons/cm3/TeV
    double xmin = 1.e-3, xmax = 1.;    // The min/max energy of product as a
                                       // function of proton energy.
    double ymin, ymax, dy, y;          // The exponent of the above
    double Epimin, Epimax;             // Min/Max pion energy in TeV for the integral
                                       // over all pion energies.
    double dw, w;                      // exponent of the above
    double sinel;                      // The inelastic part of the total cross-section of pp
                                       // collisions.
    double qpi;                        // The production rate of pions.
    double Fpi;           /* The integrated quantity from eq. 78 from Kelner et al. 2006,
                             namely, the emissivity of electrons for one particular pion
                             energy.*/
    double sum;           // For the integral over all pion energies.
    double Phig;          /* Φ_gamma the energy spectral distribution in #/cm3/TeV/sec
                             for   gamma-rays produced by neutral pion decay.	*/
    double Fg;            // Spectrum of γ rays from pion decay. Eq 58 from Kelner et
                          // al. 2006
    double transition;    // The transition between delta approximation and
                          // distributions in TeV

    double dopfac_cj;
    dopfac_cj = dopfac * (1. - beta * cos(angle)) / (1. + beta * cos(angle));

    ymin = std::log10(xmin);
    ymax = std::log10(xmax);
    dy = (ymax - ymin) / (N - 1);

    transition = 0.10;    // The transition between delta approximation and
                          // distributions.

    size_t size = en_phot.size();
    for (size_t j = 0; j < size; j++) {
        Eg = en_phot[j] * 1.0e-12 * constants::erg;
        if (Eg < transition) {    // Delta approximation for distribution
            Epimin = Eg + constants::mpionTeV * constants::mpionTeV /
                              (4. * Eg);    // The min pion energy in TeV for the integral
            Epimax = 1.e6;                  // The max pion energy in Tev for the integral
            dw = std::log10(Epimax / Epimin) / (N - 1);    // The logarithmic step with which the
                                                           // pion energy increases.
            sum = 0.;
            for (int i = 1; i < N; i++) {    // I am using 1 instead of 0 because I get a weird
                                             // spine otherwise
                                             // at x~10^-4 (production of particles from protons
                                             // with total energy less than the rest mass --
                                             // impossible)
                w = std::log10(Epimin) + i * dw;
                Ep = constants::mprotTeV + std::pow(10., w) / constants::Kpi;
                sinel = sigma_pp(Ep);
                Jp = proton_dist(gammap_min, Ep, Epcode_max, spline_Jp, acc_Jp);
                qpi = 2. * ntilde / constants::Kpi * sinel *
                      Jp;    // The production rate of neutral pions
                Fpi = qpi * std::pow(10., w) /
                      sqrt(std::pow(10., (2. * w)) -
                           constants::mpionTeV * constants::mpionTeV);    // eq 78 in Kelner+2006
                sum += dw * Fpi;
            }
            Phig = constants::cee * pp_targets * sum * constants::mbarn *
                   std::log(10.);    // dNg/dEg in #/TeV/cm3/sec
        } else if ((Eg > transition) && (Eg <= Epcode_max)) {
            sum = 0.;
            for (int i = 1; i < N; i++) {
                y = ymin + (i - 1) * dy;
                Ep = std::pow(10., (std::log10(Eg) - y));
                if ((Ep >= 0.1) && (Ep <= Epcode_max)) {
                    sinel = sigma_pp(Ep);
                    Jp = proton_dist(gammap_min, Ep, Epcode_max, spline_Jp, acc_Jp);
                    Fg = gspec_pp(Ep, y);
                    sum += dy * (sinel * Jp * Fg);
                }
            }
            Phig = constants::cee * pp_targets * sum * constants::mbarn *
                   std::log(10.);    // dNg/dEg in #/TeV/cm3/sec
        }    // end of if argument for Eg<0.1 or Eg>0.1.
        else {    // in case the energy of the photon fails the limits ==> no
                  // gamma-rays are produced
            Phig = 1.e-50;
        }
        num_phot[j] = Phig * constants::herg * vol * Eg;             // erg/s/Hz per segment
        en_phot_obs[j] = en_phot[j] * dopfac;                        //*dopfac;
        num_phot_obs[j] = num_phot[j] * std::pow(dopfac, dopnum);    //*dopfac;	//L'_v' -> L_v
        if (counterjet == true) {
            en_phot_obs[j + size] = en_phot[j] * dopfac_cj;
            num_phot_obs[j + size] = num_phot[j] * std::pow(dopfac_cj, dopnum);
        }
    }    // end of loop for photon energies
}    // End of function that produces the gamma-rays produced by neutral pion
     // decay from pp interactions

//************************************************************************************************************

// The following are common with electrons/neutrinos:

double set_ntilde(double p) {
    // n tilde from Kelner et al. 2006. It is the number of produced pions for a
    // given proton distribution:
    double ntilde;
    if (p <= 2.25) {
        ntilde = 1.10;
    } else if (p >= 2.75) {
        ntilde = 0.86;
    } else {
        ntilde = 0.91;
    }
    return ntilde;
}

double target_protons(double ntot_prot, double ntargets, double plfrac) {
    double pp_targets;
    // if ((1.-plfrac)*ntot_prot>=ntargets || (ntargets>(1.-plfrac)*ntot_prot &&
    // (1.-plfrac)*ntot_prot>0.01*ntargets)){

    pp_targets = ntargets;
    if ((1. - plfrac) * ntot_prot > 0.01 * ntargets) {
        // in case of a dense jet
        // pp_targets = (1.-plfrac)*ntot_prot + ntargets;
        pp_targets += ntot_prot;
    }
    return pp_targets;
}

double sigma_pp(double Ep) {    // cross section of pp in mb (that's why I
                                // multiply with 1.e-27 at Phie)

    double Ethres = 1.22e-3;    // threshold energy (in TeV) for pp interactions (= 1.22GeV)
    double L = std::log(Ep);    // L = ln(Ep/1TeV) as definied in Kelner et al. 2006
                                // for the cross section
    double sinel;               // σ_inel in mb

    sinel = 1.e-50;
    if (Ep >= Ethres) {
        sinel = (34.3 + 1.88 * L + 0.25 * L * L) * (1. - std::pow((Ethres / Ep), 4)) *
                (1. - std::pow((Ethres / Ep), 4));
    }
    return sinel;
}

double proton_dist(double gpmin, double Ep, double Epcode_max, gsl_spline* spline_Jp,
                   gsl_interp_accel* acc_Jp) {
    double fp, gp;

    fp = 1.e-80;
    if ((Ep / constants::mprotTeV >= gpmin) && (Ep <= Epcode_max)) {
        gp = Ep / constants::mprotTeV;                  // from TeV to γ
        fp = gsl_spline_eval(spline_Jp, gp, acc_Jp);    // dn/dγ
    }
    return fp / constants::mprotTeV;    // the distribution of protons in #/cm3/TeV
}

double gspec_pp(double Ep, double y) {

    double L = std::log(Ep);    // L = ln(Ep/1TeV) as definied in Kelner et al. 2006
                                // for the cross section
    double Betag, bg, ykg;      /* The sub-functions that describe the function F_e(x,E_p) that
                                   implies the number of electrons in the interval (x,x+dx) per
                                   collision. In particular, eqs. 59-61 from Kelner et al. 2006.*/
    double Fg;                  // The spectrum of gamma-rays produced from pion decay. Eq 58
                                // from Kelner +06

    Betag = 1.30 + 0.14 * L + 0.011 * L * L;
    bg = 1.0 / (1.79 + 0.11 * L + 0.008 * L * L);
    ykg = 1.0 / (0.801 + 0.049 * L + 0.014 * L * L);

    Fg = Betag * y * std::log(10.) / (std::pow(10., y)) *
         std::pow(((1. - std::pow(10., (y * bg))) /
                   (1. + ykg * std::pow(10., (y * bg)) * (1. - std::pow(10., (y * bg))))),
                  4) *
         (1. / (y * std::log(10.)) -
          (4. * bg * std::pow(10., (y * bg)) / (1. - std::pow(10., (y * bg)))) -
          (4. * ykg * bg * std::pow(10., (y * bg)) * (1. - 2. * std::pow(10., (y * bg)))) /
              (1. + ykg * std::pow(10., (y * bg)) * (1. - std::pow(10., (y * bg)))));

    return Fg;
}

//************************************************************************************************************
void sum_photons(size_t nphot, std::vector<double>& en_perseg, std::vector<double>& lum_perseg,
                 size_t ntarg, const std::vector<double>& targ_en,
                 const std::vector<double>& targ_lum) {

    std::vector<double> lx(ntarg, 0.0);    // std::log10 of targ_en[]/emerg
    std::vector<double> lL(ntarg, 0.0);    // std::log10 of Luminosity of targets in erg/s/Hz

    double logx;

    for (size_t i = 0; i < ntarg;
         i++) {    // lx = std::log10(hv/mec2) of target photons with energy hv
        lx[i] = std::log10(targ_en[i] / constants::emerg);
        if (targ_lum[i] == 0.) {
            lL[i] = -100.;
        } else {
            lL[i] = std::log10(targ_lum[i]);
        }
    }

    // We interpolate over the targets
    gsl_interp_accel* acc_targ = gsl_interp_accel_alloc();
    gsl_spline* spline_targ = gsl_spline_alloc(gsl_interp_akima, ntarg);
    gsl_spline_init(spline_targ, lx.data(), lL.data(), ntarg);

    for (size_t i = 0; i < nphot; i++) {
        logx = std::log10(en_perseg[i] / constants::emerg);
        if (logx >= lx[0] && logx <= lx[ntarg - 5]) {
            lum_perseg[i] +=
                std::max(1.e-200, std::pow(10., gsl_spline_eval(spline_targ, logx, acc_targ)));
        }
    }

    gsl_spline_free(spline_targ), gsl_interp_accel_free(acc_targ);
}

void sum_photons(size_t nphot, const std::vector<double>& en_perseg,
                 std::vector<double>& lum_perseg, size_t ntarg, const std::vector<double>& targ_en,
                 const std::vector<double>& targ_lum) {

    std::vector<double> lx(ntarg, 0.0);    // std::log10 of targ_en[]/emerg
    std::vector<double> lL(ntarg, 0.0);    // std::log10 of Luminosity of targets in erg/s/Hz

    double logx;

    for (size_t i = 0; i < ntarg;
         i++) {    // lx = std::log10(hv/mec2) of target photons with energy hv
        lx[i] = std::log10(targ_en[i] / constants::emerg);
        if (targ_lum[i] == 0.) {
            lL[i] = -100.;
        } else {
            lL[i] = std::log10(targ_lum[i]);
        }
    }

    // We interpolate over the targets
    gsl_interp_accel* acc_targ = gsl_interp_accel_alloc();
    gsl_spline* spline_targ = gsl_spline_alloc(gsl_interp_akima, ntarg);
    gsl_spline_init(spline_targ, lx.data(), lL.data(), ntarg);

    for (size_t i = 0; i < nphot; i++) {
        logx = std::log10(en_perseg[i] / constants::emerg);
        if (logx >= lx[0] && logx <= lx[ntarg - 5]) {
            lum_perseg[i] +=
                std::max(1.e-200, std::pow(10., gsl_spline_eval(spline_targ, logx, acc_targ)));
        }
    }

    gsl_spline_free(spline_targ), gsl_interp_accel_free(acc_targ);
}

//************************************************************************************************************
void Grays::set_grays_pg(double gp_min, double gp_max, gsl_interp_accel* acc_Jp,
                         gsl_spline* spline_Jp, std::vector<double>& en_perseg,
                         std::vector<double>& lum_perseg, size_t nphot) {

    size_t N = 10;
    double mpion =
        137.5e6 / constants::erg / (constants::cee * constants::cee);    // mass of pion in g
    double eta, deta;           // eta parameter: η = 4εE_p/(m_p^2*c^4) and its step
    double eta_zero = 0.313;    // eq 16 from Kelner & Aharonian 08
    double Hg;                  // the quantity that we intergrate for every eta
    double dNdEg;               // spectrum of gamma-rays in #/erg/cm3/sec
    double eta_max = 99.99;     // max η
    double eta_min = 1.10;      // min η
    double nu_min = en_perseg[0] / constants::herg;            // the min freq of photon targets
    double nu_max = en_perseg[nphot - 1] / constants::herg;    // the max freq of photon targets
    std::vector<double> freq(nphot, 0.0);     // frequency of photons per segment in Hz
    std::vector<double> Uphot(nphot, 0.0);    // diff energy density per segment in #/cm3/erg

    double dopfac_cj;
    dopfac_cj = dopfac * (1. - beta * cos(angle)) / (1. + beta * cos(angle));

    for (size_t k = 0; k < nphot; k++) {
        freq[k] = en_perseg[k] / constants::herg;    // Hz from erg
        Uphot[k] =
            lum_perseg[k] * (r / constants::cee /
                             (constants::herg * constants::herg * freq[k] * vol));    // #/cm3/erg
    }

    // Interpolation for jet photon distribution
    gsl_interp_accel* acc_ng = gsl_interp_accel_alloc();
    gsl_spline* spline_ng = gsl_spline_alloc(gsl_interp_akima, nphot);
    gsl_spline_init(spline_ng, freq.data(), Uphot.data(), nphot);

    deta = std::log10(eta_max / eta_min) / static_cast<double>(N - 1);
    size_t size = en_phot.size();
#pragma omp parallel for private(eta, Hg, dNdEg)    // possibly lost: 9,424 bytes in 31 blocks
    for (size_t i = 0; i < size; i++) {             // for every produced γ ray energy
        double Eg = en_phot[i];                     // in erg
        if (Eg > mpion * constants::cee * constants::cee) {
            double sum = 0.0;
            gsl_integration_workspace* w1 = gsl_integration_workspace_alloc(100);
            double result1, error1;
            gsl_function F1;
            for (size_t j = 0; j < N; j++) {
                eta =
                    eta_zero * (std::pow(10., std::log10(eta_min) + static_cast<double>(j) * deta));
                auto F1params = HetagParams{eta,    eta_zero,  Eg,     gp_min,
                                            gp_max, spline_Jp, acc_Jp,    // product,
                                            acc_ng, spline_ng, nu_min, nu_max};
                F1.function = &Hetag;
                F1.params = &F1params;
                double max =
                    std::log10(Eg / (gp_min * constants::pmgm * constants::cee * constants::cee));
                double min =
                    std::log10(Eg / (gp_max * constants::pmgm * constants::cee * constants::cee));
                gsl_integration_qag(&F1, min, max, 1e0, 1e0, 100, 1, w1, &result1, &error1);
                Hg = std::pow(constants::pmgm * constants::cee * constants::cee, 2) / 4. * result1;
                sum += Hg * deta * eta *
                       std::log(10.);    // todo: replace std::log(10) with std::M_LN10
            }
            dNdEg = sum;    // in #/erg/cm3/sec
            gsl_integration_workspace_free(w1);
        } else {
            dNdEg = 1.e-100;    // in #/erg/cm3/sec
        }

        num_phot[i] = dNdEg * constants::herg * en_phot[i] * vol;    // erg/sec/Hz
        en_phot_obs[i] = en_phot[i] * dopfac;
        num_phot_obs[i] = num_phot[i] * std::pow(dopfac, dopnum);    // L'_v' -> L_v
        if (counterjet == true) {
            en_phot_obs[i + size] = en_phot[i] * dopfac_cj;
            num_phot_obs[i + size] = num_phot[i] * std::pow(dopfac_cj, dopnum);
        }
    }

    gsl_spline_free(spline_ng);
    gsl_interp_accel_free(acc_ng);
}

double Hetag(double x, void* pars) {
    // eq 70 from KA08 for photons and writen as 0< x=Eg/Ep <1
    HetagParams* params = static_cast<HetagParams*>(pars);
    double eta = params->eta;
    double eta_zero = params->eta_zero;
    double Eg = params->Eg;
    double gp_min = params->gp_min;
    double gp_max = params->gp_max;
    gsl_spline* spline_Jp = params->spline_Jp;
    gsl_interp_accel* acc_Jp = params->acc_Jp;
    gsl_interp_accel* acc_ng = params->acc_ng;
    gsl_spline* spline_ng = params->spline_ng;
    double nu_min = params->nu_min;
    double nu_max = params->nu_max;

    double fp;         // number density of accelerated protons in #/cm3/erg
    double fph;        // number density of target photons in #/cm3/Hz
    double Phig;       // spectrum of gamma rays
    double Ep;         // proton energy Ep=x/Eg
    double fph_jet;    // number density of target photons in #/cm3/Hz

    Ep = Eg / std::pow(10., x);
    fp = colliding_protons(spline_Jp, acc_Jp, gp_min, gp_max, Ep);
    fph_jet = photons_jet(eta, Ep, spline_ng, acc_ng, nu_min, nu_max);
    fph = fph_jet;    // all the target photons are now included into the
                      // fph_jet function
    Phig = PhiFunc_gamma(eta, eta_zero, std::pow(10., x));
    return fp * fph * Phig * std::log(10.) / Eg;
}

//************************************************************************************************************

// The following are common with electrons & neutrinos:

double colliding_protons(gsl_spline* spline_Jp, gsl_interp_accel* acc_Jp, double gp_min,
                         double gp_max,
                         double Ep) {    // the number density of the protons
    double Jp;    // number density of protons from interpolation for energy Ep
    double gp;    // the Lorentz factor of energy Ep
    gp = Ep / (constants::pmgm * constants::cee * constants::cee);

    Jp = 1.0e-100;
    if (gp >= gp_min && gp <= gp_max) {
        Jp = gsl_spline_eval(spline_Jp, gp, acc_Jp);    // in #/γ/cm3
    }
    return Jp / (constants::pmgm * constants::cee * constants::cee);    // in #/erg/cm3
}

double photons_jet(double eta, double Ep, gsl_spline* spline_ng, gsl_interp_accel* acc_ng,
                   double nu_min, double nu_max) {
    // dif number density of photons inside the jet (from other physical
    // processes)
    double nu_g;    // the freq of the colliding photon

    nu_g = eta * std::pow(constants::pmgm * constants::cee * constants::cee, 2) / (4. * Ep) /
           constants::herg;    // epsilon from KA08 over herg
    if (nu_g >= nu_min && nu_g <= nu_max) {
        return gsl_spline_eval(spline_ng, nu_g, acc_ng);
    } else {
        return 1.e-100;
    }
}

double PhiFunc_gamma(double eta, double eta0, double x) {
    // eqs 27,28,29 etc for gamma rays with interpolation in the tables given by
    // KA08 and eqs 31,32,33,34,35,36,37,38,39,40,41 and tables for leptons

    double Phi;                  // spectrum of products
    double xminus, xplus;        // eq. 19 from KA08 for min/max energy of pion
    double r = 0.146;            // mpion/mproton = 0.146 for above expressions
    double s, delta, Beta;       // the parameters for spectrum
    double xeta = eta / eta0;    // x_eta to distinguish from x = Ep/Ei -- it's
                                 // ρ sometimes in KA08

    tables_photomeson_gamma(s, delta, Beta, xeta);

    xplus = 1. / (2. + 2. * eta) *
            (eta + r * r + sqrt((eta - r * r - 2. * r) * (eta - r * r + 2. * r)));
    xminus = 1. / (2. + 2. * eta) *
             (eta + r * r - sqrt((eta - r * r - 2. * r) * (eta - r * r + 2. * r)));

    double y;
    y = (x - xminus) / (xplus - xminus);
    if (x > xminus && x < xplus) {
        Phi = Beta * std::exp(-s * std::pow(std::log(x / xminus), delta)) *
              std::pow(std::log(2. / (1. + y * y)), 2.5 + 0.4 * std::log(eta / eta0));
    } else if (x <= xminus) {
        Phi = Beta * std::pow(std::log(2.), 2.5 + 0.4 * std::log(eta / eta0));
    } else {
        Phi = 1.e-100;
    }
    return Phi;
}

//! The tables from KA08 for photomeson that give s,δ and B
void tables_photomeson_gamma(double& s, double& delta, double& Beta, double xeta) {
    // Gamma rays from neutral pion decay:
    size_t size = etagTable.size();
    gsl_interp_accel* acc_sigma = gsl_interp_accel_alloc();
    gsl_spline* spline_sigma = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_interp_accel* acc_delta = gsl_interp_accel_alloc();
    gsl_spline* spline_delta = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_interp_accel* acc_Beta = gsl_interp_accel_alloc();
    gsl_spline* spline_Beta = gsl_spline_alloc(gsl_interp_cspline, size);

    gsl_spline_init(spline_sigma, etagTable.data(), sgTable.data(), size);
    gsl_spline_init(spline_delta, etagTable.data(), deltagTable.data(), size);
    gsl_spline_init(spline_Beta, etagTable.data(), BetagTable.data(), size);

    s = gsl_spline_eval(spline_sigma, xeta, acc_sigma);
    delta = gsl_spline_eval(spline_delta, xeta, acc_delta);
    Beta = gsl_spline_eval(spline_Beta, xeta, acc_Beta);

    gsl_spline_free(spline_sigma), gsl_interp_accel_free(acc_sigma);
    gsl_spline_free(spline_delta), gsl_interp_accel_free(acc_delta);
    gsl_spline_free(spline_Beta), gsl_interp_accel_free(acc_Beta);
}

}    // namespace kariba
