#pragma once
#include <cmath>
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
    int cutoff_type;
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
    int cutoff_type;
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


enum CutoffType { //trying something to make this easier to read, 
    // it will convert directly to int when used 
    Exponential = 0,
    SuperExponential = 1,
    Sech2 = 2

};

class Particles {
  protected:
    double mass_gr;     //!< particle mass in grams
    double mass_kev;    //!< same as above but in keV, using electrons as "reference"

    std::vector<double> p;        //!< array of particle momenta
    std::vector<double> ndens;    //!< array of number density per unit volume, per unit momentum
    std::vector<double> gamma;    //!< array of particle kinetic energies for each momentum
    std::vector<double> gdens;    //!< array of number density per unit volume, per unit gamma
    std::vector<double> pdensp2_diff_logp;    //!< array with differential of number
                                       //!< density for radiation calculation p^-2*dn/dp
    int cutoff_type = Exponential; // setting cutoff type for legacy bhjet without cutoff switch 
  public:
    Particles(size_t size);
    // cutoff adjustments: 
    int get_cutoff_type() const { return cutoff_type; }
    virtual void set_cutoff_type(int t) { cutoff_type = t; }
    static double cutoff_factor(double x, int type = Exponential);

    virtual void set_mass(double m);
    virtual void initialize_gdens();
    virtual void initialize_pdens();
    virtual void differentiate();

    virtual const std::vector<double>& get_p() const { return p; }
    virtual const std::vector<double>& get_pdens() const { return ndens; }
    virtual const std::vector<double>& get_gamma() const { return gamma; }
    virtual const std::vector<double>& get_gdens() const { return gdens; }
    virtual const std::vector<double>& get_pdensp2_diff_logp() const { return pdensp2_diff_logp; }

    virtual double count_particles();
    virtual double count_particles_energy();
    virtual double av_p();
    virtual double av_gamma();
    virtual double av_psq();
    virtual double av_gammasq();

    virtual void test_arrays();
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

inline double Particles::cutoff_factor(double x, int type) 
{
    switch (type) {
    case Exponential: // classic exponential cutoff 
        return std::exp(-x); 
    case SuperExponential:
        return std::exp(-(x * x)); 
    case Sech2: {  // sech^2 cutoff
        const double ax = std::fabs(x);
        if (ax < 30.0) { //cosh(30) ~ 5e12
            const double c = std::cosh(ax); //directly calculate sech(x)
            return 1.0 / (c * c);
        } // this is to avoid issues with cosh, maybe should replace this with taylor exp .
        const double t = std::exp(-2.0 * ax); //approx
        const double denom = 1.0 + t;
        return 4.0 * t / (denom * denom);
    }

    default:
        return std::exp(-x);
    }
}

}    // namespace kariba
