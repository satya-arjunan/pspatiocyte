#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main(int argc, char *argv[]) {
  const int nx(960); // lattice dimensions
  const int ny(960);
  const int nz(960);
  const bool verbose(true);
  const double duration(1.01); // s
  const double rv(5e-9);
  const int num_A(0);
  const int num_B(64000);
  const int num_C(64000);
  const double D(10e-12); // m^2/s (diffusion coefficient)
  const int nlogs(5000);

  World world(argc, argv, nx, ny, nz, rv);
  world.set_output_numbers_logspace(-4, 1, nlogs);

  Species A("A", D, num_A, world);
  Species B("B", D, num_B, world);
  Species C("C", D, num_C, world);

  const double keff(20e-19); //effective or macroscopic rate
  const double kd(4*M_PI*2*rv*2*D);
  const double ka((keff*kd)/(kd-keff)); //intrinsic rate, used by Spatiocyte
  const double kr_eff(1.35); //reverse rate
  const double kr(ka*kr_eff/keff);

  Reaction forward("B + C -> A", B, C, ka, A, world);
  Reaction reverse("A -> B + C", A, kr, B, C, world);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}


