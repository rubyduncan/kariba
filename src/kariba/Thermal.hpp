#pragma once

#include "Particles.hpp"

namespace kariba {

//! Class for thermal particles. ndens is number density per unit momentum
class Thermal : public Particles {
  protected:
    double Temp, thnorm, theta;

  public:
    Thermal(size_t size);

    void set_p();
    void set_ndens();
    void set_temp_kev(double T);
    void set_norm(double n);

    double K2(double x);

    void test();
};

}    // namespace kariba
