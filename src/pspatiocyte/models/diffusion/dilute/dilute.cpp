#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main(int argc, char *argv[]) {
  std::string dirname("output");
  double D(10e-12); // m^2/s (diffusion coefficient)
  unsigned seed(5);

  if (argc >= 2) {
    dirname = argv[1];
  }
  if (argc >= 3) {
    D = atof(argv[2]);
  }
  if (argc >= 4) {
    seed = atoi(argv[3]);
  }

  const int nx(960); // lattice dimensions
  const int ny(960);
  const int nz(960);
  const bool verbose(true);
  const double duration(0.11); // s
  const double rv(2.5e-9);
  const unsigned nlogs(2000);
  const float end_log(duration);
  const unsigned nmolecules(8);
  const double numbers_log_interval(0); //don't log numbers
  const double coords_log_interval(-1); //log after each diffusion step

  World world(argc, argv, nx, ny, nz, rv, numbers_log_interval, dirname, seed,
              coords_log_interval);

  Species A("A", D, nmolecules, world);

  Vector<float> range(0.05, 0.05, 0.05);
  A.set_populate_range(range);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}

