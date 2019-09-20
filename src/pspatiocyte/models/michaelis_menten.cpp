#include "world.hpp"
#include "species.hpp"
#include "reaction.hpp"

int main(int argc, char *argv[]) {
  const int nx(64); // lattice dimensions
  const int ny(64);
  const int nz(64);
  const bool verbose(true);
  const double duration(10); // s
  const double volume(909e-18); // m^3
  const double rv(pow(volume/((nz*2.0)*(ny*sqrt(3.0))*
                              (nx*sqrt(8.0/3.0))), 1.0/3.0)); // voxel radius
  std::cout << "voxel radius:" << rv << std::endl;
  const int num_E(9090);
  const int num_S(90910);
  const int num_ES(0);
  const int num_P(0);
  const double D(10e-12); // m^2/s (diffusion coefficient)
  const int nlogs(100);
  const double log_interval(duration/nlogs);

  World world(argc, argv, nx, ny, nz, rv);
  Species E("E", D, num_E, world);
  Species S("S", D, num_S, world);
  Species ES("ES", D, num_ES, world);
  Species P("P", D, num_P, world);

  const double keff(0.01e-18); //effective or macroscopic rate
  const double kd(4*M_PI*2*rv*2*1e-12);
  const double ka((keff*kd)/(kd-keff)); //intrinsic rate, used by Spatiocyte
  const double kr(1); //reverse rate
  const double kcat(1); //production rate

  Reaction forward("E + S -> ES", E, S, ka, ES, world);
  Reaction reverse("ES -> E + S", ES, kr, E, S, world);
  Reaction product("ES -> E + P", ES, kcat, E, P, world);

  world.initialize();
  world.run(log_interval, duration, verbose);
  return 0;
}

