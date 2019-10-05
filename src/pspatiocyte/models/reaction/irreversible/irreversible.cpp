#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main(int argc, char *argv[]) {
  std::string dirname("output");
  double k(0.135); // dissociation rate
  if (argc >= 2) {
    dirname = argv[1];
  }
  if (argc >= 3) {
    k = atof(argv[2]);
  }

  const int nx(960); // lattice dimensions
  const int ny(960);
  const int nz(960);
  const bool verbose(true);
  const double duration(3.01); // s
  const double rv(5e-9);
  const int num_A(64000);
  const int num_B(0);
  const int num_C(0);
  const double D(10e-12); // m^2/s (diffusion coefficient)
  const int nlogs(500);
  const double log_interval(duration/nlogs);

  World world(argc, argv, nx, ny, nz, rv, log_interval, dirname);

  Species A("A", D, num_A, world);
  Species B("B", D, num_B, world);
  Species C("C", D, num_C, world);

  Reaction product("A -> B + C", A, k, B, C, world);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}


