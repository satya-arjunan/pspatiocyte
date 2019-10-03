#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main(int argc, char *argv[]) {
  const int nx(960); // lattice dimensions
  const int ny(960);
  const int nz(960);
  const bool verbose(true);
  const double duration(0.1); // s
  const double rv(5e-9);
  const int num_A(0);
  const int num_B(64000);
  const int num_C(64000);
  const double D(10e-12); // m^2/s (diffusion coefficient)
  const int nlogs(100);
  const double log_interval(duration/nlogs);

  World world(argc, argv, nx, ny, nz, rv, log_interval);
  Species A("A", D, num_A, world);
  Species B("B", D, num_B, world);
  Species C("C", D, num_C, world);

  const double keff(2.51e-19); //effective or macroscopic rate
  const double kd(4*M_PI*2*rv*2*D);
  const double ka((keff*kd)/(kd-keff)); //intrinsic rate, used by Spatiocyte
  const double kb(1.35); //reverse rate

  Reaction forward("B + C -> A", B, C, ka, A, world);
  Reaction reverse("A -> B + C", A, kb, B, C, world);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}


