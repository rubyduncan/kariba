#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "kariba/Neutrinos_pp.hpp"
#include "kariba/Radiation.hpp"
#include "kariba/constants.hpp"

namespace kariba {

Neutrinos_pp::Neutrinos_pp(size_t size, double Emin, double Emax) : Radiation(size) {

    en_phot_obs.resize(2 * en_phot_obs.size(), 0.0);
    num_phot_obs.resize(2 * num_phot_obs.size(), 0.0);

    double einc = std::log10(Emax / Emin) / static_cast<double>(size - 1);
    for (size_t i = 0; i < size; i++) {
        en_phot[i] = std::pow(10., std::log10(Emin) + static_cast<double>(i) * einc);
        en_phot_obs[i] = en_phot[i];
        en_phot_obs[i + size] = en_phot[i];
    }
}

void Neutrinos_pp::set_neutrinos_pp(double pspec, double gammap_min, double gammap_max,
                                    double ntot_prot, double nwind, double plfrac,
                                    gsl_interp_accel* acc_Jp, gsl_spline* spline_Jp,
                                    const std::string& outputConfiguration,
                                    const std::string& flavor, int infosw,
                                    std::string_view source) {

    std::ofstream NeutrinosppFile;    // for plotting
    if (infosw >= 2) {
        std::string filepath;
        if (source.compare("JET") == 0) {
            filepath = outputConfiguration + "/Output/Neutrinos/" + flavor + "_pp.dat";
        } else {
            std::cerr << "Wrong source; cannot be " << source << " but rather JET!" << std::endl;
            exit(1);
        }
        NeutrinosppFile.open(filepath, std::ios::app);
    }

    double ntilde = multiplicity(pspec);    // The number of produced pions for
                                            // a given proton distribution
    double pp_targets = target_protons(ntot_prot, nwind, plfrac);
    double Epcode_max = gammap_max * constants::pmgm * constants::cee * constants::cee *
                        constants::erg * 1.e-12;    // The proton energy in TeV

    int N = 60;    // Steps of the secondary particle (e.g., pion) energy
    double xmin = 1.e-3,
           xmax = 1.;                // min/max energy of decaying particle in Ep
    double dy, y;                    // std::log of above
    double Bprob = 0.;               // The probapility of production for pions
    double Ev;                       // The energy of the neutrino in TeV
    double Epimin, Epimax = 1.e6;    // The min/max energy of the pions and the
                                     // respective exponents
    double lEpi;                     // std::log10 of pion energy in TeV
    double dw;                       // The logarithmic step with which the pion energy increases
    double Ep;                       // The energy of the proton	in TeV
    double sum;                      // for the integrals
    double sinel;                    // The inelastic cross-section of pp collisions
    double Jp;                       // The number density of protons from interpolation
    double qpi;                      // The production rate of pions
    double fv;                       // for the neutrino production (eq. 36 from Kelner+2006)
    double Fv;                       // for the neutrino production (eq. 66 from Kelner+2006)
    double Phiv;                     // The neutrino rate in #/TeV/cm3/sec
    double Fnuspec;                  // Spectrum of muon neutrinos eq66 from Kelner+2006
    double transition = 0.0;         // the transition from delta fuctions to distribution in TeV
    int i_init = 0;                  // the first nerutrino energy

    if (flavor.compare("muon") == 0) {
        transition = 0.01;
        i_init = 3;    // from 3 otherwise I get to <Ep=1GeV
        Bprob = prob();
    } else if (flavor.compare("electron") == 0) {
        transition = 0.05;
        Bprob = prob_fve();
    }
    dy = std::log10(xmax / xmin) / (N - 1);

    for (size_t j = 0; j < en_phot.size(); j++) {     // for every single Ev
        Ev = en_phot[j] * constants::erg * 1.e-12;    // in TeV
        if (Ev <= transition) {
            Epimin = Ev + constants::mpionTeV * constants::mpionTeV / (4. * Ev);
            dw = (std::log10(Epimax / Epimin)) / (N - 1);
            sum = 1.e-100;
            for (int i = i_init; i < N; i++) {
                lEpi = std::log10(Epimin) + i * dw;
                Ep = constants::mprotTeV + std::pow(10., lEpi) / constants::Kpi;
                sinel = sigma_pp(Ep);
                Jp = proton_dist(gammap_min, Ep, Epcode_max, spline_Jp, acc_Jp);
                qpi = 2. * ntilde / constants::Kpi * sinel *
                      Jp;    // The production rate of neutral pions
                fv = distr_pp(std::log10(Ev), lEpi, flavor);
                // Fv =
                // qpi*std::pow(10.,lEpi)/sqrt(std::pow(10.,(2.*lEpi))-mpionTeV*mpionTeV)*fv*Bprob;
                Fv = qpi * std::pow(10., lEpi) / sqrt(std::pow(10., (2. * lEpi))) * fv * Bprob;
                sum += dw * Fv;
            }    // end of if statement for energies greater than Ep_min
            Phiv = constants::cee * pp_targets * sum * 1.e-27 * std::log(10.);
        }    // end of for loop for all the pions
        else if ((Ev > transition) && (Ev <= Epcode_max)) {
            sum = 1.e-100;
            for (int i = 0; i < N; i++) {
                y = std::log10(xmin) + i * dy;
                Ep = std::pow(10., (std::log10(Ev) - y));
                if (Ep >= .1 && Ep <= Epcode_max) {
                    sinel = sigma_pp(Ep);
                    Jp = proton_dist(gammap_min, Ep, Epcode_max, spline_Jp, acc_Jp);
                    Fnuspec = secondary_spectrum(Ep, y, flavor);
                    sum += dy * (sinel * Jp * Fnuspec);
                }
            }    // end of for loop for all the pions
            Phiv = constants::cee * pp_targets * sum * 1.e-27 * std::log(10.);
        } else {
            Phiv = 1.e-100;
        }
        num_phot[j] = Phiv * constants::herg * vol * Ev;    // erg/s/Hz per segment
        en_phot_obs[j] = en_phot[j] * dopfac;               // *dopfac;
        num_phot_obs[j] =
            num_phot[j] *
            std::pow(dopfac,
                     dopnum);    // dopfac*dopfac;			//L'_v' -> L_v

        if (infosw >= 2) {
            NeutrinosppFile << std::left << std::setw(15) << Ev * 1.e12 / constants::erg
                            << std::setw(25) << Phiv / (1.e12 / constants::erg) << std::setw(25)
                            << num_phot[j] / (constants::herg * en_phot[j]) << std::endl;
        }
    }    // End of for loop for all the neutrino energies.
    if (infosw >= 2) {
        NeutrinosppFile.close();
    }
}    // End of function that produces the neutrinos from pp

//***********************************************************************************************************

double prob_fve() {    // it is the same as of electrons

    size_t N = 20;                  // The steps of integration.
    double r = .573;                // r = 1-λ = m_μ^2/m_p^2 = 0.573.
    double xmin = 0., xmax = 1.;    // x = E_{particle}/E_p.
    double x, dx = (xmax - xmin) / static_cast<double>(N - 1);
    double gn, hn1, hn2, fve;    // (equations 40-43 from Kelner et al. 2006)
    double sum = 0.;
    double Bprob;

    for (size_t i = 0; i < N; i++) {
        x = xmin + static_cast<double>(i) * dx;
        gn = 2. / (3. * (1. - r) * (1. - r)) *
             ((1. - x) * (6. * (1. - x) * (1. - x) + r * (5. + 5. * x - 4. * x * x)) +
              6. * r * std::log(x));
        hn1 = 2. / (3. * (1. - r) * (1. - r)) *
              ((1. - r) * (6. - 7. * r + 11. * r * r - 4. * r * r * r) + 6. * r * std::log(r));
        hn2 = 2. * (r - x) / (3. * r * r) *
              (7. * r * r - 4. * r * r * r + 7. * x * r - 4. * x * r * r - 2. * x * x -
               4. * x * x * r);

        // The function f_ve is given by equation 40 from Kelner et al. 2006
        if (x >= r) {    // H(x-r) = 1 and H(r-x) = 0
            fve = gn;
        } else {    // H(x-r) = 0 and H(r-x) = 1
            fve = hn1 + hn2;
        }
        sum += fve * dx;
    }
    Bprob = 1. / sum;
    return Bprob;
}

double distr_pp(double lEv, double lEpi, std::string_view flavor) {
    double rmasses = .573;                     // r = 1-λ = m_μ^2/m_p^2 = 0.573.The ratio of muon
                                               // and proton energies
    double k = std::pow(10., (lEv - lEpi));    // x=Ev/Epion
    double Fvespec = 0.;                       // The spectrum of secondary electrons from pion
                                               // decay. Eq 62 from Kelner+06
    if (flavor.compare("muon") == 0) {
        double lamda = 1. - rmasses;    // 1-rmasses
        double g0, gn, hn1, h0, hn2;    // for the neutrino production (eq.
                                        // 37-39 from Kelner et al. 2006)
        double fn2, fn1, fn;

        g0 = (3. - 2. * rmasses) / (9. * (1. - rmasses) * (1. - rmasses));
        gn = g0 * (9. * k * k - 6. * std::log(k) - 4. * k * k * k - 5.);
        hn1 = g0 * (9. * rmasses * rmasses - 6. * std::log(rmasses) -
                    4. * rmasses * rmasses * rmasses - 5.);
        h0 = (1. + 2. * rmasses) * (rmasses - k) / (9. * rmasses * rmasses);
        hn2 = h0 * (9. * (rmasses + k) - 4. * (rmasses * rmasses + rmasses * k + k * k));

        if (k >= rmasses) {
            fn2 = gn;
        } else {
            fn2 = hn1 + hn2;
        }
        if (k <= lamda) {
            fn1 = 1. / lamda;
            fn = fn1 + fn2;
        } else {
            fn = fn2;
        }
        Fvespec = fn;
    } else if (flavor.compare("electron") == 0) {
        double gn, hn1, hn2, fve;    //(eqs. 40-43 from  Kelner et al. 2006)

        gn = 2. / (3. * (1. - rmasses) * (1. - rmasses)) *
             ((1. - k) * (6. * (1. - k) * (1. - k) + rmasses * (5. + 5. * k - 4. * k * k)) +
              6. * rmasses * std::log(k));
        hn1 = 2. / (3. * (1. - rmasses) * (1. - rmasses)) *
              ((1. - rmasses) * (6. - 7. * rmasses + 11. * rmasses * rmasses -
                                 4. * rmasses * rmasses * rmasses) +
               6. * rmasses * std::log(rmasses));
        hn2 = 2. * (rmasses - k) / (3. * rmasses * rmasses) *
              (7. * rmasses * rmasses - 4. * rmasses * rmasses * rmasses + 7. * k * rmasses -
               4. * k * rmasses * rmasses - 2. * k * k - 4. * k * k * rmasses);

        // The function f_v is given by equation 36 from Kelner et al. 2006 and
        // is: f_e(x) = g_v*H(x-r) + (h_v1(x) + h_v2(x))*H(r-x), with H(y) the
        // Heaviside function.
        if (k >= rmasses) {    // H(x-r) = 1 and H(r-x) = 0.
            fve = gn;
        } else {    // H(x-r) = 0 and H(r-x) = 1.
            fve = hn1 + hn2;
        }
        Fvespec = fve;
    }
    return Fvespec;
}

double secondary_spectrum(double Ep, double y, std::string_view flavor) {
    double L = std::log(Ep);    // L = ln(Ep/1TeV) as definied in Kelner et al. 2006
                                // for the cross section
    double Fvespec = 0.;        // The spectrum of secondary electrons from pion
                                // decay. Eq 62 from Kelner+06
    if (flavor.compare("muon") == 0) {
        double Betav2, bv2, kv2;      // for the neutrino production (eq. 63-69
                                      // from Kelner et al. 2006)
        double Fi1, Fi2, Fi3, Fv2;    // for the neutrino production (eq. 62
                                      // from Kelner et al. 2006)
        double kv1, bv1, Betav1;      // for the neutrino production (eq. 63-69
                                      // from Kelner et al. 2006)
        double om;                    // energy of neutrino energy over the proton energy
        double F1, F2, F3, F4, F5, F6,
            Fv1;    // for the neutrino production (eq. 66 from Kelner et al.
                    // 2006)

        Betav2 = 1. / (69.5 + 2.65 * L + .3 * L * L);
        bv2 = 1.0 / std::pow((.201 + .062 * L + .00042 * L * L), .25);
        kv2 = (.279 + .141 * L + .0172 * L * L) / (.3 + (2.3 + L) * (2.3 + L));

        Fi1 = 1. + kv2 * y * std::log(10.) * y * std::log(10.);
        Fi2 = 1. + .3 / std::pow(10., (y * bv2));
        Fi3 = -y * std::log(10.);

        Fv2 = Betav2 * Fi1 * Fi1 * Fi1 / (std::pow(10., y) * Fi2) * Fi3 * Fi3 * Fi3 * Fi3 * Fi3;

        kv1 = 1.07 - .086 * L + .002 * L * L;
        bv1 = 1. / (1.67 + .111 * L + .0038 * L * L);
        Betav1 = 1.75 + .204 * L + .010 * L * L;

        om = std::pow(10., y) / .427;

        if (om <= 1.) {
            F1 = 1. - std::pow(om, bv1);
            F2 = 1. + kv1 * std::pow(om, bv1) * F1;
            F3 = 4. * bv1 * std::pow(om, bv1) / F1;
            F5 = 4. * kv1 * bv1 * std::pow(om, bv1) * (1. - 2. * std::pow(om, bv1));
            F4 = F5 / (1. + kv1 * std::pow(om, bv1) * F1);
            F6 = 1. / (std::log(om)) - F3 - F4;
            Fv1 = Betav1 * std::log(om) / (om) * (F1 / F2) * (F1 / F2) * (F1 / F2) * (F1 / F2) * F6;
        } else {
            Fv1 = 0.;
        }
        Fvespec = Fv1 + Fv2;
    } else if (flavor.compare("electron") == 0) {
        double Betae, be, yke;    // The sub-functions that describe the function F_ve(x,E_p)
                                  // that implies the number of elec neutrinos in the interval
                                  // (x,x+dx) per collision. In particular, eqs. 63-65 from
                                  // Kelner et al. 2006.

        // The sub-functions that describe the function F_e(x,E_p) that implies
        // the number of electrons in the interval (x,x+dx) per collision. In
        // particular, eq. 63 from Kelner et al. 2006:
        Betae = 1. / (69.5 + 2.65 * L + .3 * L * L);
        be = 1.0 / std::pow((.201 + .062 * L + .00041 * L * L),
                            .25);    // Eq. 64 from K06.

        yke =
            (.279 + .141 * L + .0172 * L * L) / (.3 + (2.3 + L) * (2.3 + L));    // Eq. 65 from k06

        Fvespec = Betae * (1. + yke * (y * std::log(10.)) * (y * std::log(10.))) *
                  (1. + yke * (y * std::log(10.)) * (y * std::log(10.))) *
                  (1. + yke * (y * std::log(10.)) * (y * std::log(10.))) /
                  (std::pow(10., y) * (1. + .3 / std::pow(10., (y * be)))) * (-y * std::log(10.)) *
                  (-y * std::log(10.)) * (-y * std::log(10.)) * (-y * std::log(10.)) *
                  (-y * std::log(10.));
    }
    return Fvespec;
}

}    // namespace kariba
