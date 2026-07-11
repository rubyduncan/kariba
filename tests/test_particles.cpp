#include "doctest.h"

#include <cmath>
#include <kariba/Bknpower.hpp>
#include <kariba/Kappa.hpp>
#include <kariba/Mixed.hpp>
#include <kariba/Powerlaw.hpp>
#include <kariba/Thermal.hpp>
#include <kariba/constants.hpp>

namespace karcst = kariba::constants;

TEST_CASE("Particles subclass functionality") {

    SUBCASE("Cutoff type is actually being set within the code:") {
        kariba::Powerlaw powerlaw(100);
        kariba::Mixed mixed(100);
        kariba::Bknpower bknpower(100);

        // see that it is able to pick up the "hardcoded" backwards compatible type
        CHECK(powerlaw.get_cutoff_type() == 0);
        CHECK(mixed.get_cutoff_type() == 0);
        CHECK(bknpower.get_cutoff_type() == 0);

        //this should check that it can be overwritten from the setting: 
        powerlaw.set_cutoff_type(2);
        mixed.set_cutoff_type(2);
        bknpower.set_cutoff_type(2);

        CHECK(powerlaw.get_cutoff_type() == 2);
        CHECK(mixed.get_cutoff_type() == 2);
        CHECK(bknpower.get_cutoff_type() == 2);

        // check both the functions are actually doing something? 
        
        const double x = 2.0;
        const double c = std::cosh(x); 
        CHECK(kariba::Particles::cutoff_factor(x, 0) == doctest::Approx(std::exp(-x)));
        CHECK(kariba::Particles::cutoff_factor(x, 2) == doctest::Approx(1.0 / (c * c)));
    }


    SUBCASE("Thermal particle distribution") {
        kariba::Thermal thermal(100);

        SUBCASE("Temperature setting and momentum arrays") {
            double temp_kev = 511.0;
            thermal.set_temp_kev(temp_kev);
            thermal.set_p();

            CHECK(!thermal.get_p().empty());
            CHECK(!thermal.get_gamma().empty());

            // Check that gamma values are reasonable for thermal distribution
            const std::vector<double>& gamma = thermal.get_gamma();
            CHECK(gamma[0] >= 1.0);
            CHECK(gamma[99] > gamma[0]);
        }

        SUBCASE("Number density normalization") {
            double temp_kev = 100.0;
            double norm = 1.0;

            thermal.set_temp_kev(temp_kev);
            thermal.set_p();
            thermal.set_norm(norm);
            thermal.set_ndens();

            CHECK(!thermal.get_pdens().empty());
            CHECK(!thermal.get_gdens().empty());

            // Check particle count is positive
            double total_particles = thermal.count_particles();
            CHECK(total_particles > 0.0);
        }

        SUBCASE("Average quantities") {
            double temp_kev = 511.0;
            thermal.set_temp_kev(temp_kev);
            thermal.set_p();
            thermal.set_norm(1.0);
            thermal.set_ndens();

            double avg_gamma = thermal.av_gamma();
            double avg_p = thermal.av_p();

            CHECK(avg_gamma >= 1.0);
            CHECK(avg_p > 0.0);

            // For thermal distribution, average gamma should be reasonable
            CHECK(avg_gamma < 10.0);    // Not too relativistic for 511 keV
        }

        SUBCASE("Gamma-momentum relationship") {
            double temp_kev = 511.0;
            thermal.set_temp_kev(temp_kev);
            thermal.set_p();
            thermal.set_norm(1.0);
            thermal.set_ndens();

            // Test that gamma arrays are consistent
            const std::vector<double>& gamma = thermal.get_gamma();
            const std::vector<double>& p = thermal.get_p();

            // Check relativistic energy-momentum relation for first few points
            for (size_t i = 0; i < 5; i++) {
                double expected_gamma =
                    std::sqrt(1.0 + std::pow(p[i] / (karcst::emgm * karcst::cee), 2.0));
                CHECK(std::abs(gamma[i] - expected_gamma) < 1e-10);
            }
        }
    }

    SUBCASE("Powerlaw particle distribution") {
        kariba::Powerlaw powerlaw(100);

        SUBCASE("Basic powerlaw setup") {
            double pmin = 1e-3 * karcst::emgm * karcst::cee;
            double gmax = 1e3;
            double spec = 2.5;

            powerlaw.set_p(pmin, gmax);
            powerlaw.set_pspec(spec);
            powerlaw.set_norm(1.0);
            powerlaw.set_ndens();

            CHECK(!powerlaw.get_p().empty());
            CHECK(!powerlaw.get_gamma().empty());

            const std::vector<double>& gamma = powerlaw.get_gamma();
            CHECK(gamma[0] >= 1.0);
            // CHECK(gamma[99] <= gmax * 1.001);    // Allow small numerical error

            double total_particles = powerlaw.count_particles();
            CHECK(total_particles > 0.0);
        }

        SUBCASE("Average quantities for powerlaw") {
            double pmin = 1e-3 * karcst::emgm * karcst::cee;
            double gmax = 1e4;

            powerlaw.set_p(pmin, gmax);
            powerlaw.set_pspec(2.0);
            powerlaw.set_norm(1.0);
            powerlaw.set_ndens();

            double avg_gamma = powerlaw.av_gamma();
            CHECK(avg_gamma > 1.0);
            CHECK(avg_gamma < gmax);
        }

        SUBCASE("Cooling steady state") {
            double pmin = 1e-3 * karcst::emgm * karcst::cee;
            double gmax = 1e3;

            powerlaw.set_p(pmin, gmax);
            powerlaw.set_pspec(2.5);
            powerlaw.set_norm(1.0);
            powerlaw.set_ndens();

            // Test cooling calculation doesn't crash
            double ucom = 0.0;    // No external Compton
            double n0 = 1.0;
            double bfield = 1e3;
            double r = 1e10;
            double tshift = 0.1;

            CHECK_NOTHROW(powerlaw.cooling_steadystate(ucom, n0, bfield, r, tshift));
        }
    }

    SUBCASE("Mixed particle distribution") {
        kariba::Mixed mixed(100);

        SUBCASE("Mixed distribution setup") {
            double temp_kev = 511.0;
            double gmax = 1e3;
            double plfrac = 0.1;
            double pspec = 2.0;

            mixed.set_temp_kev(temp_kev);
            mixed.set_p(gmax);
            mixed.set_plfrac(plfrac);
            mixed.set_pspec(pspec);
            mixed.set_norm(1.0);
            mixed.set_ndens();

            CHECK(!mixed.get_p().empty());
            CHECK(!mixed.get_gamma().empty());

            double total_particles = mixed.count_particles();
            CHECK(total_particles > 0.0);
        }

        SUBCASE("Mixed distribution cooling") {
            double temp_kev = 511.0;
            double gmax = 1e3;

            mixed.set_temp_kev(temp_kev);
            mixed.set_p(gmax);
            mixed.set_plfrac(0.1);
            mixed.set_pspec(2.0);
            mixed.set_norm(1.0);
            mixed.set_ndens();

            // Test cooling doesn't crash
            CHECK_NOTHROW(mixed.cooling_steadystate(0.0, 1.0, 1e3, 1e10, 0.1));
        }
    }

    SUBCASE("Kappa particle distribution") {
        kariba::Kappa kappa(100);

        SUBCASE("Kappa distribution setup") {
            double temp_kev = 511.0;
            double gmax = 1e3;
            double kappa_val = 3.0;

            kappa.set_temp_kev(temp_kev);
            kappa.set_p(gmax);
            kappa.set_kappa(kappa_val);
            kappa.set_norm(1.0);
            kappa.set_ndens();

            CHECK(!kappa.get_p().empty());
            CHECK(!kappa.get_gamma().empty());

            double total_particles = kappa.count_particles();
            CHECK(total_particles > 0.0);
        }
    }

    SUBCASE("Broken powerlaw particle distribution") {
        kariba::Bknpower bknpower(100);

        SUBCASE("Broken powerlaw setup") {
            double pmin = 1e-3 * karcst::emgm * karcst::cee;
            double pbrk = 1e-1 * karcst::emgm * karcst::cee;
            double gmax = 1e3;

            bknpower.set_pspec1(-2.0);
            bknpower.set_pspec2(2.5);
            bknpower.set_p(pmin, pbrk, gmax);
            bknpower.set_norm(1.0);
            bknpower.set_ndens();

            CHECK(!bknpower.get_p().empty());
            CHECK(!bknpower.get_gamma().empty());

            double total_particles = bknpower.count_particles();
            CHECK(total_particles > 0.0);
        }
    }
}

TEST_CASE("Particle distribution consistency checks") {
    SUBCASE("Momentum-gamma consistency") {
        kariba::Thermal thermal(50);
        thermal.set_temp_kev(100.0);
        thermal.set_p();
        thermal.set_norm(1.0);
        thermal.set_ndens();

        const std::vector<double>& p = thermal.get_p();
        const std::vector<double>& gamma = thermal.get_gamma();

        // Check relativistic energy-momentum relation: E^2 = (pc)^2 + (mc^2)^2
        for (size_t i = 0; i < 10; i++) {
            double expected_gamma =
                std::sqrt(1.0 + std::pow(p[i] / (karcst::emgm * karcst::cee), 2.0));
            CHECK(std::abs(gamma[i] - expected_gamma) < 1e-10);
        }
    }

    SUBCASE("Array bounds and monotonicity") {
        kariba::Powerlaw powerlaw(50);
        double pmin = 1e-3 * karcst::emgm * karcst::cee;
        double gmax = 1e3;

        powerlaw.set_p(pmin, gmax);
        powerlaw.set_pspec(2.0);
        powerlaw.set_norm(1.0);
        powerlaw.set_ndens();

        const std::vector<double>& p = powerlaw.get_p();
        const std::vector<double>& gamma = powerlaw.get_gamma();

        // Check arrays are monotonically increasing
        for (size_t i = 1; i < 50; i++) {
            CHECK(p[i] > p[i - 1]);
            CHECK(gamma[i] > gamma[i - 1]);
        }

        // // Check maximum gamma is respected (allow small numerical error)
        // CHECK(gamma[49] <= gmax * 1.001);
    }
}
