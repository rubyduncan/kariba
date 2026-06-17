#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include "kariba/Electrons.hpp"
#include "kariba/Particles.hpp"
#include "kariba/Powerlaw.hpp"
#include "kariba/constants.hpp"

namespace kariba {

//! Class constructor to initialize object
Powerlaw::Powerlaw(size_t size) : Particles(size) {
    plnorm = 1.;

    mass_gr = constants::emgm;
    mass_kev = constants::emgm * constants::gr_to_kev;
}

//! Methods to set momentum/energy arrays
void Powerlaw::set_p(double min, double ucom, double bfield, double betaeff, double r, double fsc) {
    pmin = min;
    pmax = max_p(ucom, bfield, betaeff, r, fsc);

    double pinc = (std::log10(pmax) - std::log10(pmin)) / static_cast<double>(p.size() - 1);

    for (size_t i = 0; i < p.size(); i++) {
        p[i] = std::pow(10., std::log10(pmin) + static_cast<double>(i) * pinc);
        gamma[i] = std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    }
}

void Powerlaw::set_p(double min, double gmax) {
    pmin = min;
    pmax = std::pow(std::pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * constants::cee;

    double pinc = (std::log10(pmax) - std::log10(pmin)) / static_cast<double>(p.size() - 1);

    for (size_t i = 0; i < p.size(); i++) {
        p[i] = std::pow(10., std::log10(pmin) + static_cast<double>(i) * pinc);
        gamma[i] = std::pow(std::pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
    }
}

//! Method to set differential electron number density from known pspec,
//! normalization, and momentum array
void Powerlaw::set_ndens() {
    for (size_t i = 0; i < ndens.size(); i++) {
        ndens[i] = plnorm * std::pow(p[i], -pspec) * std::exp(-p[i] / pmax);
    }
    initialize_gdens();
    differentiate();
}

//! methods to set the slope and normalization
void Powerlaw::set_pspec(double s1) { pspec = s1; }

void Powerlaw::set_norm(double n) {
    plnorm =
        n * (1. - pspec) / (std::pow(p[p.size() - 1], (1. - pspec)) - std::pow(p[0], (1. - pspec)));
}

//! Injection function to be integrated in cooling
double injection_pl_int(double x, void* pars) {
    InjectionPlParams* params = static_cast<InjectionPlParams*>(pars);
    double s = params->s;
    double n = params->n;
    double m = params->m;
    double max = params->max;

    double mom_int = std::pow(std::pow(x, 2.) - 1., 1. / 2.) * m * constants::cee;

    return n * std::pow(mom_int, -s) * std::exp(-mom_int / max);
}

//! Method to solve steady state continuity equation. NOTE: KN cross section not
//! included in IC cooling
void Powerlaw::cooling_steadystate(double ucom, double n0, double bfield, double r,
                                   double betaeff) {
    double Urad = std::pow(bfield, 2.) / (8. * constants::pi) + ucom;
    double pdot_ad = betaeff * constants::cee / r;
    double pdot_rad = (4. * constants::sigtom * constants::cee * Urad) /
                      (3. * mass_gr * std::pow(constants::cee, 2.));
    double tinj = r / (constants::cee);

    double integral, error;
    gsl_function F1;
    auto params = InjectionPlParams{pspec, plnorm, mass_gr, pmax};
    F1.function = &injection_pl_int;
    F1.params = &params;

    for (size_t i = 0; i < gamma.size(); i++) {
        if (i < gamma.size() - 1) {
            gsl_integration_workspace* w1;
            w1 = gsl_integration_workspace_alloc(100);
            gsl_integration_qag(&F1, gamma[i], gamma[i + 1], 1e1, 1e1, 100, 1, w1, &integral,
                                &error);
            gsl_integration_workspace_free(w1);

            ndens[i] =
                (integral / tinj) / (pdot_ad * p[i] / (mass_gr * constants::cee) +
                                     pdot_rad * (gamma[i] * p[i] / (mass_gr * constants::cee)));
        } else {
            ndens[gamma.size() - 1] =
                ndens[gamma.size() - 2] *
                std::pow(p[gamma.size() - 1] / p[gamma.size() - 2], -pspec - 1) * std::exp(-1.);
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
//! acceleration and cooling timescales The estimate is identical to the old
//! agnjet but in momentum space; see Lucchini et al. 2019 for the math of the
//! old version
double Powerlaw::max_p(double ucom, double bfield, double betaeff, double r, double fsc) {
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

// Methods to set energy array for protons
void Powerlaw::set_energy(double gpmin, double fsc, double f_beta, double bfield, double r_g,
                          double z, double r, int infosw, double protdens, double ntargets,
                          double Uradjet, const std::string& outputConfiguration,
                          const std::string& source) {

    double gpmax, logdgp;
    if (fsc < 1.) {
        double pp_targets = ntargets + protdens;

        ProtonTimescales(logdgp, fsc, f_beta, bfield, gpmin, gpmax, r_g, z, r, infosw, pp_targets,
                         Uradjet, outputConfiguration, source);

    } else {
        isEfficient = true;
        gpmax = fsc;
        logdgp =
            std::log10(2. * gpmax / gpmin) /
            static_cast<double>(gamma.size() - 1);    // so as to extend a bit further from γ_p,max
        check_secondary_charged_syn(bfield, gpmax);
    }

    for (size_t i = 0; i < gamma.size(); i++) {
        gamma[i] = std::pow(10., (std::log10(gpmin) + static_cast<double>(i) * logdgp));
        p[i] = sqrt(gamma[i] * gamma[i] - 1.) * mass_gr * constants::cee;
    }
    pmin = p[0];
    pmax = sqrt(gpmax * gpmax - 1.) * mass_gr * constants::cee;
}

void Powerlaw::ProtonTimescales(double& logdgp, double fsc, double f_beta, double bfield,
                                double gpmin, double& gpmax, double r_g, double z, double r,
                                int infosw, double pp_targets, double Uradjet,
                                const std::string& outputConfiguration, const std::string& source) {

    double Tacc0;    // Tacc = Ep*Tacc,0 ==> Tacc,0 = 1/(ηecB)
    double paramG,
        paramH;        // The parameters used to solve Emax equation, i.e.,
                       //  Emax^2 + paramG*Emax - paramH =0 for protons
    double betap;      // the velocity in c of the proton
    double sinel;      // the pp cross section from Kelner et al.2006
    double tescape;    // timescale of the escape especially for the low energy
                       // ones
    double EpTeV;      // energy of the proton in TeV
    double gp;         // Lorentz factor of non-thermal protons

    std::ofstream timescalesFile;

    Tacc0 = 1. / (3. * fsc / 4. * constants::charg * constants::cee * bfield);    // acceleration
    double Tesc = r / (f_beta * constants::cee);                                  // proton escape
    double Tpp = 1. / (constants::sigmapp * pp_targets * constants::cee);         // pp

    // The only photon field here is the companion's/external boosted radiation
    double Tpg0 = 1. / (9.38 * 7.4e-17 * Uradjet / (mass_gr * constants::cee * constants::cee));

    double Tsynp0 = 6. * constants::pi / (constants::sigtom * bfield * bfield) *
                    std::pow(mass_gr * constants::cee, 4) /
                    (constants::emgm * constants::emgm * constants::cee);    // proton synchrotron

    // The parameters used to solve Emax equation, i.e., Emax^2 + paramG*Emax -
    // paramH = 0 for protons:
    paramG = (1. / Tesc + 1. / Tpp) / (1. / Tpg0 + 1. / Tsynp0);
    paramH = (1. / Tacc0) / (1. / Tpg0 + 1. / Tsynp0);

    // Maximun proton energy in erg as calculated by losses:
    double Epmax = (-paramG + sqrt(paramG * paramG + 4. * paramH)) / 2.;

    gpmax = Epmax / (constants::pmgm * constants::cee * constants::cee) + 1.;
    if (Epmax > constants::pmgm * constants::cee * constants::cee) {
        isEfficient = true;
        logdgp =
            std::log10(2. * gpmax / gpmin) /
            static_cast<double>(gamma.size() - 1);    // so as to extend a bit further from γ_p,max
    } else {
        isEfficient = false;
        logdgp = std::log10(10.) /
                 static_cast<double>(gamma.size() - 1);    // make an array between 1 and 10 GeV
    }
    if (infosw >= 2) {
        std::string filepath =
            outputConfiguration + "/Output/Particles/timescales_" + source + ".dat";
        timescalesFile.open(filepath, std::ios::app);
        for (size_t i = 0; i < gamma.size(); i++) {
            gp = std::pow(10., (std::log10(gpmin) + static_cast<double>(i) * logdgp));
            betap = sqrt(gp * gp - 1.) / gp;
            tescape = r / (constants::cee * betap * f_beta);

            EpTeV = constants::pmgm * constants::cee * constants::cee * gp / 1.6;
            sinel = sigma_pp(EpTeV);

            // we set the correct Tpg0 in order to plot the timescales;
            Tpg0 =
                (gp * mass_gr * constants::cee * constants::cee * constants::erg >= 4.8e14)
                    ? 1. / (9.38 * 7.4e-17 * Uradjet / (mass_gr * constants::cee * constants::cee))
                    : 1.e100;
            Tpp = 1. / (constants::Kpp * constants::mbarn * sinel * pp_targets * constants::cee);

            timescalesFile << std::left << std::setw(13) << z / r_g << std::setw(13) << gp
                           << std::setw(13)
                           << Tacc0 * mass_gr * constants::cee * constants::cee * gp
                           << std::setw(13) << tescape << std::setw(13)
                           << Tsynp0 /
                                  (mass_gr * constants::cee * constants::cee * gp * betap * betap)
                           << std::setw(13) << Tpp << std::setw(15)
                           << Tpg0 / (mass_gr * constants::cee * constants::cee) / gp;
            timescalesFile << std::endl;
        }
        timescalesFile.close();
        check_secondary_charged_syn(bfield, gpmax);
    }
}

void Powerlaw::check_secondary_charged_syn(double bfield, double gpmax) {
    if ((bfield * gpmax) >
        7.8e11 * 10.) {    // I *10 because the γ_p,max is in the exponential cutoff
        std::cout << "The pion synchrotron should be taken into account" << std::endl;
        std::cout << "\tB =" << std::setw(15) << bfield << std::setw(10)
                  << " g_p,max=" << std::setw(15) << gpmax << std::endl;
    }
    if ((bfield * gpmax) > 5.6e10 * 10.) {
        std::cout << "The muon synchrotron should be taken into account" << std::endl;
        std::cout << "\tB =" << std::setw(15) << bfield << std::setw(10)
                  << " g_p,max=" << std::setw(15) << gpmax << std::endl;
    }
}

double Powerlaw::sigma_pp(double Ep) {    // cross section of pp in mb (that's
                                          // why I multiply with 1.e-27 at Phie)
    double Ethres = 1.22e-3;    // threshold energy (in TeV) for pp interactions (= 1.22GeV)
    double L = std::log(Ep);    // L = ln(Ep/1TeV) as definied in Kelner et al. 2006
                                // for the cross section
    double sinel;               // σ_inel in mb

    sinel = (Ep >= Ethres) ? (34.3 + 1.88 * L + 0.25 * L * L) * (1. - std::pow((Ethres / Ep), 4)) *
                                 (1. - std::pow((Ethres / Ep), 4))
                           : 1.e-80;
    return sinel;
}

double Powerlaw::set_normprot(double nprot) {
    double Epmin = mass_gr * constants::cee * constants::cee * gamma[0];
    double Epmax = sqrt((pmax * constants::cee) * (pmax * constants::cee) +
                        (mass_gr * constants::cee * constants::cee) *
                            (mass_gr * constants::cee * constants::cee));

    // If break energy is higher than maximum energy, set normalisation for
    // uncooled distribution:
    return nprot * (1. - pspec) / (std::pow(Epmax, 1. - pspec) - std::pow(Epmin, 1. - pspec));
}

//! Method to set differential proton number density per γ from known pspec,
//! normalization and γ array
void Powerlaw::set_gdens(double r, double protdens, double nwind, double bfield, double plfrac,
                         double Uradjet) {
    // if plfrac_p>0, namely a frac of thermal of protons accelerate
    // plnormprot is in #/cm3/erg/sec and gdens is in #/cm3/gamma_p
    if (isEfficient) {
        double plnormprot = set_normprot(protdens) * plfrac * constants::cee / r;
        double Tsynp0, Tchar, Tpp;
        double betap;
        double gpmax =
            sqrt(pmax * pmax / (mass_gr * constants::cee * mass_gr * constants::cee) + 1.);
        Tsynp0 = 6. * constants::pi / (constants::sigtom * bfield * bfield) *
                 std::pow(mass_gr * constants::cee, 4) /
                 (constants::emgm * constants::emgm * constants::cee);    // proton synchrotron
        double Tpg0 = 1. / (9.38 * 7.4e-17 * Uradjet / (mass_gr * constants::cee * constants::cee));
        for (size_t i = 0; i < gamma.size(); i++) {
            betap = sqrt(gamma[i] * gamma[i] - 1.) / gamma[i];
            Tpp =
                1. / (nwind *
                      sigma_pp(constants::pmgm * constants::cee * constants::cee * gamma[i] / 1.6) *
                      betap * constants::cee *
                      constants::mbarn);    // proton energy in TeV for the function
            Tchar = std::pow(constants::cee / r + 1. / Tpp +
                                 gamma[i] * mass_gr * constants::cee * constants::cee *
                                     (1. / Tpg0 + 1. / Tsynp0),
                             -1);
            gdens[i] = plnormprot *
                       std::pow(gamma[i] * mass_gr * constants::cee * constants::cee, -pspec) *
                       std::exp(-gamma[i] / gpmax) * Tchar * mass_gr * constants::cee *
                       constants::cee;
        }
    } else {    // in case Emax<Emin
        for (size_t i = 0; i < gdens.size(); i++) {
            gdens[i] = 1.e-100;
        }
    }
}

// if (Lumsw == 1)
void Powerlaw::set_gdens_pdens(double r, double beta, double Ljet, double ep, double pspec,
                               double& protdens) {
    // Sets the normalization of the accelerated protons assuming that a fraction
    // --ep--  of the jet power goes into the proton acceleration, unless Emax<Emin
    // so set to zero
    if (isEfficient) {
        double G_jet = 1. / sqrt(1. - beta * beta);    // bulk Lorentz factor
        double plnormprot;                             // in #/cm3
        double gpmax =
            sqrt(pmax * pmax / (mass_gr * constants::cee * mass_gr * constants::cee) + 1.);
        double sum = 0;
        double dx = std::log10(gamma[2] / gamma[1]);
        for (size_t i = 0; i < gamma.size(); i++) {
            sum += std::log(10.) * dx * std::pow(gamma[i], -pspec + 2) * exp(-gamma[i] / gpmax);
        }

        plnormprot = ep * Ljet /
                     (mass_gr * constants::cee * constants::cee * constants::pi * r * r * G_jet *
                      beta * constants::cee * sum);
        for (size_t i = 0; i < gdens.size(); i++) {
            gdens[i] = plnormprot * std::pow(gamma[i], -pspec) * exp(-gamma[i] / gpmax);
        }

        sum = 0.;
        for (size_t i = 0; i < gdens.size(); i++) {
            sum += std::log(10.) * gdens[i] * gamma[i] * dx;
        }
        protdens = sum;
    } else {
        for (size_t i = 0; i < gdens.size(); i++) {
            gdens[i] = 1.e-100;
        }
        protdens = 1.e-100;
    }
}

void Powerlaw::set_gdens(double& plfrac_p, double Up, double protdens) {
    // If mass-loading, I use the specific enthalpy to work the normalisation
    if (isEfficient) {
        double gpmax =
            sqrt(pmax * pmax / (mass_gr * constants::cee * mass_gr * constants::cee) + 1.);
        double sum = 0;
        double dx = std::log10(gamma[2] / gamma[1]);
        for (size_t i = 0; i < gamma.size(); i++) {
            sum += std::log(10.) * std::pow(gamma[i], -pspec + 2.) * exp(-gamma[i] / gpmax) * dx;
        }
        double K = std::max(Up / (sum * mass_gr * constants::cee * constants::cee), 0.);

        sum = 0.;
        for (size_t i = 0; i < gamma.size(); i++) {
            sum += std::log(10.) * std::pow(gamma[i], -pspec + 1.) * exp(-gamma[i] / gpmax) * dx;
        }
        double n_nth = K * sum;
        plfrac_p = n_nth / protdens;

        for (size_t i = 0; i < gdens.size(); i++) {
            gdens[i] = K * std::pow(gamma[i], -pspec) * exp(-gamma[i] / gpmax);
        }
    } else {
        plfrac_p = 0;
        for (size_t i = 0; i < gdens.size(); i++) {
            gdens[i] = 1.e-100;
        }
    }
}

//! Function that produces the secondary electrons from pp
void Powerlaw::set_pp_elecs(gsl_interp_accel* acc_Jp, gsl_spline* spline_Jp, double ntot_prot,
                            double nwind, double plfrac, double gammap_min, double Ep_max,
                            double bfield, double r) {

    double ntilde = multiplicity(pspec);    // The number of produced pions for
                                            // a given proton distribution
    double pp_targets = target_protons(ntot_prot, nwind, plfrac);
    double Epcode_max = Ep_max * constants::erg * 1.e-12;    // The proton energy in TeV

    const size_t N = 60;                        // The number of steps of the secondary particle
                                                // (e.g., pion) energy.
    const double gmin = 1.002;                  // the min Lorentz factor of the secondary
    double gmax = Ep_max / constants::emerg;    // the max Lorentz factor of the secondary
    double Ep;                                  // the energy of the proton
    double Bprob;    // Energy distribution of sec electrons for arbitrary pion
                     // distribution
    double xmin = 1.e-3,
           xmax = 1.;    // Min/Max sec particle energy in proton energy
    double ymin, ymax, dy,
        y;                    // The exponent of min and max from above, and the step.
    double Ee;                // energy of secondary electrons in TeV
    double Epimin, Epimax;    // Min/Max pion energy in TeV for the integral
                              // over all pion energies.
    double lEpi;              // std::log10 of pion energy in TeV
    double dw;                // The logarithmic step with which the pion energy increases
    double sum;               // For the integral over all pion energies.
    double sinel;             // The inelastic part of the total cross-section of pp
                              // collisions.
    double qpi;               // The production rate of pions.
    double Jp;                // The number of non-thermal protons/cm3/TeV
    double Fespec;            // Spectrum of sec electrons from pion decay eq62 from
                              // Kelner et al. 2006
    double Fpi;               // eq. 78 from Kelner et al. 2006, the emissivity of electrons
    double Phie;              // Φ_e the energy spectral distribution in #/cm3/TeV/sec for
                              // secondary e
    double fe;                // The energy distribution/probability f_e of electrons
    double tchar;             // characteristic timescale of electron distribution
    double beta_elec;         // beta veloscity of electron
    double tesc;              // escape timescale
    double tsyne;             // electron synchrotron timescale
    double transition;        // The transition between delta approximation and
                              // distributions in TeV

    ymin = std::log10(xmin);         // The exponent of the min energy of the secondary particles.
    ymax = std::log10(xmax);         // The exponent of the max energy of the secondary particles.
    dy = (ymax - ymin) / (N - 1);    // The step of the above

    Bprob = prob();    // The probability for electron production after charged
                       // pion decay.

    transition = 0.16;    // The transition between delta approximation and
                          // distributions.

    // Loop for every electron energy
    for (size_t j = 0; j < gamma.size(); j++) {
        gamma[j] =
            std::pow(10., std::log10(gmin) + static_cast<double>(j) * std::log10(gmax / gmin) /
                                                 static_cast<double>(gamma.size() - 1));
        Ee = gamma[j] * constants::emerg * constants::erg * 1.e-12;    // in TeV
        if (Ee < transition) {
            Epimin = Ee + constants::mpionTeV * constants::mpionTeV / (4. * Ee);
            Epimax = 1.e6;
            dw = (std::log10(Epimax / Epimin)) / (N - 1);
            sum = 0.;
            for (size_t i = 0; i < N; i++) {
                lEpi = std::log10(Epimin) +
                       static_cast<double>(i) * dw;    // The exponent of the pion energy.
                Ep = constants::mprotTeV +
                     std::pow(10., lEpi) / constants::Kpi;    // The mass of the proton in TeV.
                sinel = sigma_pp(Ep);
                Jp = proton_dist(gammap_min, Ep, Epcode_max, spline_Jp, acc_Jp);
                qpi = 2. * ntilde / constants::Kpi * sinel * Jp;
                fe = elec_dist_pp(std::log10(Ee), lEpi);    // eq36 KAB16
                //				Fpi = qpi*std::pow(10.,w)/ sqrt(
                // std::pow(10.,(2.*w))- mpionTeV*mpionTeV)*fe*Bprob; Use the
                // expression below because the above makes a spike at around
                // 1e8eV (disc. wiht Maria)
                Fpi = qpi * std::pow(10., lEpi) / sqrt(std::pow(10., (2. * lEpi))) * fe *
                      Bprob /**1.5*/;
                sum += dw * (Fpi);
            }
            Phie = constants::cee * pp_targets * sum * 1.e-27 *
                   std::log(10.);    // eq 78 in #/cm3/TeV/sec
        } else if ((Ee >= transition) && (Ee <= Epcode_max)) {
            sum = 0.;
            for (size_t i = 0; i < N; i++) {    // The loop over all pion energies.
                y = ymin + static_cast<double>(i) * dy;
                Ep = std::pow(10., (std::log10(Ee) - y));    // Proton enrergy in TeV.
                if (Ep <= Epcode_max) {
                    sinel = sigma_pp(Ep);
                    Jp = proton_dist(gammap_min, Ep, Epcode_max, spline_Jp, acc_Jp);
                    Fespec = elec_spec_pp(Ep, y);
                    sum += dy * (sinel * Jp * Fespec);
                }
                Phie = constants::cee * pp_targets * sum * 1.e-27 * std::log(10.);
            }
        } else {
            Phie = 1.e-50;
        }

        beta_elec = sqrt(gamma[j] * gamma[j] - 1.) / gamma[j];
        tesc = r / (beta_elec * constants::cee);
        tsyne = 6. * constants::pi * constants::emerg /
                (constants::sigtom * constants::cee * bfield * bfield * gamma[j] * beta_elec *
                 beta_elec);
        tchar = std::pow(1. / tsyne + 1. / tesc, -1.);

        gdens[j] = Phie * tchar * Ee / gamma[j];
    }
}

//! The method to set the secondary electrons from pg (I have called the Neutrino
//! object first because I have all the tables from KA09 in this class/file
void Powerlaw::set_pg_electrons(const std::vector<double>& energy,
                                const std::vector<double>& density, double f_beta, double r,
                                double vol, double B) {
    // the density is in erg/s/Hz (because it's a Radiation object)
    double tcool;    // cooling time to account for synchrotron losses
    for (size_t i = 0; i < gamma.size(); i++) {
        gamma[i] = energy[i] / constants::emerg;
        tcool = std::pow(constants::sigtom * B * B * gamma[i] * constants::cee /
                                 (6. * constants::pi * constants::emerg) +
                             f_beta * constants::cee / r,
                         -1);
        gdens[i] = density[i] / (energy[i] * constants::herg * vol) * tcool;    // in #/erg/cm3
    }
}

//! Function that produces the secondary electrons from photon-photon
//! annihilation
void Powerlaw::Qggeefunction(double r, double vol, double bfield, size_t phot_number,
                             const std::vector<double>& en_perseg,
                             const std::vector<double>& lum_perseg, double gmax) {

    double gmin = 1.002;    // the min Lorentz factor of the secondary
    double ng;              // number density of photons with energy 2gammae (#/cm3/erg)
    double Eg;              // std::log of energy(emerg) of photon that collides in order to
                            // produce pairs
    double x, dx;           // the dimensionless energy of target photon and std::log step
                            // of the integration
    double sum;
    double R_gg;    // production rate of pairs from gg (in cm3/sec)
    // double *logx		= new double[phot_number]();	//array of
    // x=std::log10(hv/mec2) double *logNgamma	= new double[phot_number]();
    // //array of std::log10(Ngamma[#/cm3/erg]) double *Ngamma		= new
    // double[phot_number]();	//array of Ngamma[#/cm3/erg]
    std::vector<double> logx;         // array of x=std::log10(hv/mec2)
    std::vector<double> logNgamma;    // array of std::log10(Ngamma[#/cm3/erg])
    std::vector<double> Ngamma;       // array of Ngamma[#/cm3/erg]

    double tchar;        // characteristic timescale
    double beta_elec;    // beta veloscity of electron
    double tsyn;         // synchrotron timescale
    double Rann;         // pair annihilation rate between cold and non-thermal
                         // electrons
    double Ne;           // number density (not per erg) of cold/target electrons
    double Lee_gg;       // losses due to pair annihilation in #/cm3/erg/sec
    // double *tgg_ee 		= new double[phot_number]();	//photon-photon
    // annihilation timescale
    std::vector<double> tgg_ee;    // photon-photon annihilation timescale
    double tcharg;                 // characteristic timescale for photons

    double Qgg_ee;    // in #/cm3/erg

    sum = 1.e-100;
    dx = std::log10(en_perseg[2] / en_perseg[1]);
    for (size_t i = 0; i < phot_number; i++) {
        logx.push_back(std::log10(en_perseg[i] / constants::emerg));
        Ngamma.push_back(lum_perseg[i] * r /
                         (constants::herg * constants::cee * en_perseg[i] * vol));
        if (Ngamma[i] <= 1.e-100) {
            Ngamma[i] = 1.e-100;
        }
        logNgamma.push_back(std::log10(Ngamma[i]));
        sum += en_perseg[i] * std::log(10.) * dx * Ngamma[i];
    }

    for (size_t i = 0; i < phot_number; i++) {
        x = std::pow(10., logx[i]);
        if (x <= 1.) {    // for values higher than one the photons are
                          // considered merely targets
            tgg_ee[i] = 1.e100;
        } else {
            R_gg = constants::sigtom * constants::cee * .652 * (x * x - 1.) / std::pow(x, 3) *
                   std::log(x);
            tgg_ee[i] = r / (constants::cee * .652 * (x * x - 1.) / std::pow(x, 3) * std::log(x));
        }

        // tcharg = std::pow(constants::cee / r + 1. / tgICS + 1. / tgg_ee[i], -1);
        tcharg = std::pow(constants::cee / r + 1. / tgg_ee[i], -1);
        Ngamma[i] *= tcharg / tgg_ee[i];    // photons taken into account for γγ->ee
        logNgamma[i] = std::log10(Ngamma[i]);
    }

    // We interpolate all over the photons added above
    gsl_interp_accel* acc_lNg = gsl_interp_accel_alloc();
    gsl_spline* spline_lNg = gsl_spline_alloc(gsl_interp_steffen, phot_number);
    gsl_spline_init(spline_lNg, logx.data(), logNgamma.data(), phot_number);

    for (size_t i = 0; i < gamma.size(); i++) {
        gamma[i] =
            std::pow(10., std::log10(gmin) + static_cast<double>(i) * std::log10(gmax / gmin) /
                                                 static_cast<double>(gamma.size() - 1));
        Eg = std::log10(2. * gamma[i]);    // 2γ from MK95

        if (Eg >= logx[1] && Eg <= logx[phot_number - 1]) {
            ng = std::pow(10., gsl_spline_eval(spline_lNg, Eg,
                                               acc_lNg));    // n_γ(2γ) from MK95
        } else {
            ng = 1.e-200;
        }

        sum = 0.;
        for (size_t j = 0; j < phot_number; j++) {    // eq.57 from MK95 in units of [n]=#/cm3/erg
            x = std::pow(10., logx[j]);
            R_gg = production_rate(gamma[i], x);
            sum += dx * std::pow(10., logNgamma[j] + logx[j]) * std::log(10.) * R_gg;
            // ln10 and the extra x term because of std::log integration
        }
        Qgg_ee = 4. * ng * sum * constants::emerg;    // #/cm3/erg/sec

        beta_elec = sqrt(gamma[i] * gamma[i] - 1.) / gamma[i];
        tsyn = 6. * constants::pi * constants::emerg /
               (constants::cee * constants::sigtom * bfield * bfield * gamma[i] * beta_elec *
                beta_elec);
        tchar = std::pow(constants::cee / r + 1. / tsyn, -1);

        Qgg_ee *= tchar;    // #/cm3/erg of pairs after photon-photon

        Rann = 3. * constants::sigtom * constants::cee / (8. * gamma[i]) *
               (std::pow(gamma[i], -0.5) + std::log(gamma[i]));
        Ne = Qgg_ee * gamma[i] * constants::emerg;
        Lee_gg = Ne * Rann * Qgg_ee * tchar;    //(non-)cooled pairs collide with cold pairs
        gdens[i] = (Qgg_ee - Lee_gg) * constants::emerg;    // #/cm3/γ
    }

    // we free the space occupied for interpolation
    gsl_spline_free(spline_lNg), gsl_interp_accel_free(acc_lNg);
}

//! simple method to check quantities.
void Powerlaw::test() {
    std::cout << "Power-law distribution;" << std::endl;
    std::cout << "pspec: " << pspec << std::endl;
    std::cout << "Array size: " << p.size() << std::endl;
    std::cout << "Default normalization: " << plnorm << std::endl;
    std::cout << "Particle mass in grams: " << mass_gr << std::endl;
}

}    // namespace kariba
