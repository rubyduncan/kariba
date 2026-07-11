// #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <cmath>
#include <gsl/gsl_spline.h>
#include <kariba/Bknpower.hpp>
#include <kariba/Compton.hpp>
#include <kariba/Cyclosyn.hpp>
#include <kariba/Kappa.hpp>
#include <kariba/Mixed.hpp>
#include <kariba/Powerlaw.hpp>
#include <kariba/ShSDisk.hpp>
#include <kariba/Thermal.hpp>
#include <kariba/constants.hpp>

const double EPS = 1.0e-5;

namespace karcst = kariba::constants;

TEST_CASE("Integration tests - Complete workflows") {
    SUBCASE("Single zone jet model workflow") {
        // This test replicates the singlezone example workflow

        size_t nel = 50;    // Smaller arrays for faster testing

        double B = 1.5e-3;
        double n = 9.5e-3;
        double gmin = 4.1e3;
        double gmax = 6.4e4;    // Reduced for faster testing
        double pmin = std::sqrt(gmin * gmin - 1.0) * karcst::emgm * karcst::cee;
        double p = 3.03;

        // Set up electron distribution
        kariba::Powerlaw electrons(nel);
        electrons.set_p(pmin, gmax);
        electrons.set_pspec(p);
        electrons.set_norm(n);
        electrons.set_ndens();

        CHECK(!electrons.get_p().empty());
        CHECK(electrons.count_particles() > 0.0);
        gmin = electrons.get_gamma()[0];
        gmax = electrons.get_gamma()[nel - 1];
        // CHECK(gmin == doctest::Approx(4100.0).epsilon(EPS));
        // CHECK(gmax == doctest::Approx(64000.0).epsilon(EPS));

        // Calculate physical quantities
        double Ue = electrons.av_gamma() * n * karcst::emgm * karcst::cee_cee;
        double Ub = B * B / (8.0 * karcst::pi);
        double equip = Ue / Ub;

        CHECK(Ue > 0.0);
        CHECK(Ub > 0.0);
        CHECK(equip > 0.0);
    }
}
