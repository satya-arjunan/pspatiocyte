#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main(int argc, char *argv[]) {
  std::string dirname("output");
  unsigned seed(5);
  float fraction(0.5);

  if (argc >= 2) {
    dirname = argv[1];
  }
  if (argc >= 3) {
    fraction = atof(argv[2]);
  }
  if (argc >= 4) {
    seed = atoi(argv[3]);
  }

  double D(10e-12); // m^2/s (diffusion coefficient)
  const int nx(960); // lattice dimensions
  const int ny(960);
  const int nz(960);
  const bool verbose(true);
  const double duration(161e-6); // s
  const double rv(2.5e-9);
  const unsigned nlogs(1000);
  const unsigned nmolecules(1);
  const double numbers_log_interval(0); //don't log numbers
  const double coords_log_interval(duration/nlogs);
  const unsigned ncrowders((nx-2)*(ny-2)*(nz-2)*fraction);

  World world(argc, argv, nx, ny, nz, rv, numbers_log_interval, dirname, seed,
              coords_log_interval);

  Species A("A", D, nmolecules, world);
  world.add_output_coords_species(A);

  Species C("C", 0, ncrowders, world);

  Vector<float> range(0.05, 0.05, 0.05);
  A.set_populate_range(range);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}

