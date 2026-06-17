#pragma once

#include <vector>

namespace kariba {

//! Structure used for GSL integration
struct PlParams {
    double s;
    double n;
};

//! Structure used for GSL integration
struct BknParams {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
    int cutoff_type; 
};

//! Structure used for GSL integration
struct ThParams {
    double t;
    double n;
    double m;
};

//! Structure used for GSL integration
struct KParams {
    double t;
    double k;
};

//! Structure used for GSL integration
struct InjectionMixedParams {
    double s;
    double t;
    double nth;
    double npl;
    double m;
    double min;
    double max;
    double cutoff;
    int cutoff_type = 0; 
};

//! Structure used for GSL integration
struct InjectionKappaParams {
    double t;
    double k;
    double n;
    double m;
};

//! Structure used for GSL integration
struct InjectionPlParams {
    double s;
    double n;
    double m;
    double max;
};

//! Structure used for GSL integration
struct InjectionBknParams {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
    double n;
    int cutoff_type; 
};

//! Template class for particle distributions
//! This class contains members and methods that are used for thermal,
//! non-thermal and mixed distributions
class Particles {
  protected:
    double mass_gr;     //!< particle mass in grams
    double mass_kev;    //!< same as above but in keV, using electrons as "reference"

    std::vector<double> p;        //!< array of particle momenta
    std::vector<double> ndens;    //!< array of number density per unit volume, per unit momentum
    std::vector<double> gamma;    //!< array of particle kinetic energies for each momentum
    std::vector<double> gdens;    //!< array of number density per unit volume, per unit gamma
    std::vector<double> gdens_diff;    //!< array with differential of number
                                       //!< density for radiation calculation
    int cutoff_type = 0; 

  public:
    Particles(size_t size);

    // cutoff adjustments: 
    int get_cutoff_type() const { return cutoff_type; }
    void set_cutoff_type(int t) { cutoff_type = t; }
    static double cutoff_factor(double x, int cutoff_type = 0); 

    void set_mass(double m);
    void initialize_gdens();
    void initialize_pdens();
    void gdens_differentiate();

    const std::vector<double>& get_p() const { return p; }

    const std::vector<double>& get_pdens() const { return ndens; }

    const std::vector<double>& get_gamma() const { return gamma; }

    const std::vector<double>& get_gdens() const { return gdens; }

    const std::vector<double>& get_gdens_diff() const { return gdens_diff; }

    double count_particles();
    double count_particles_energy();
    double av_p();
    double av_gamma();
    double av_psq();
    double av_gammasq();

    void test_arrays();
};


//need this to be accessible by both mixed and powerlaw, so here goes: 
// this has three options for the behavior of the cutoff in the electron spectrum: 
// 1. the traditional exponential cutoff exp(e/ecut)
// 2. this is from Comisso 2021, exp[(e/ecut)^2] - magnetic turbulence PIC sim 
// 3. also comisso 2021, sech[(e/ecut)^2]

// in principle, meant to act like this: (in set_ndens or injection_mixed_int), in mixed.cpp or powerlaw.cpp: 
// C = cutoff_factor(p[i]/pmax_pl, cutoff_type); then ndens = npl * mom_int thing * C

// Three cutoff types, with x = p/p_cut (momentum)
//cutoff type 0 = exp(-x)
// cutoff type 1 = exp(-x^2)

// cutoff type 2 = sech(x)^2 = 1/cosh(x)^2 
// cosh(x) = ( e^x + e^-x ) / 2
// 
inline double Particles::cutoff_factor(double x, int type) {
    if (type == 0) return std::exp(-x); // classic exponential cutoff 
    if (type == 1) return std::exp(-(x * x)); // 
    if (type == 2) {    // sech^2 cutoff
        const double ax = std::fabs(x); 
        if (ax < 30.0) { //cosh(30) ~ 5e12 
            const double c = std::cosh(ax); //directly calculate sech(x)
            return 1.0 / (c * c);
        } else { // maybe should replace this with taylor exp .
            const double t = std::exp(-2.0 * ax); //approx
            const double denom = 1.0 + t;
            return 4.0 * t / (denom * denom);
        }
    }
    return std::exp(-x); 
}


}    // namespace kariba
