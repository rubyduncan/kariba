#pragma once

#include "Radiation.hpp"

namespace kariba {

//! Class black body photons, inherited from Radiation.hh
class BBody : public Radiation {
  protected:
    double Tbb;
    double Lbb;
    double normbb;

  public:
    BBody(size_t size = 40);

    void set_temp_kev(double T);
    void set_temp_k(double T);
    void set_temp_hz(double nu);
    void set_lum(double L);
    void bb_spectrum();

    double temp_kev() const;
    double temp_k() const;
    double temp_hz() const;
    double lum() const;
    double norm() const;
    double Urad(double d) const;

    void test();
};

}    // namespace kariba
