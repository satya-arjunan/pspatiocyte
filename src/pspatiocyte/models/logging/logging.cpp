#include <time.h>
#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"
int main(int argc, char *argv[]) {
  std::string dirname("output");
  unsigned seed(5);
  const int nx(32); // lattice dimensions
  const int ny(32);
  const int nz(32);
  const bool verbose(true);
  const double duration(0.001); // s
  const double rv(2.5e-9);
  const unsigned nlogs(2);
  const float end_log(duration);
  const unsigned nmolecules(8);
  const double numbers_log_interval(0); //don't log numbers
  const double coords_log_interval(-1); //log after each diffusion step

  World world(argc, argv, nx, ny, nz, rv, numbers_log_interval, dirname, seed,
              coords_log_interval);

  const unsigned size(world.get_parallel_environment().getsize());
  std::vector<Species*> species_list;
  Vector<float> range(-1, -1, -1);
  for (unsigned i(0); i < size; ++i) {
    std::string rank(std::to_string(i));;
    //go through 8 subvolumes
    for (unsigned j(0); j < 8; ++j) {
      std::string sv(std::to_string(j));;
      Species* ghost(new Species(std::string("Ghost")+rank+sv, 0, 0, world));
      ghost->set_populate_range(range);
      species_list.push_back(ghost);
      /*
      Species* out(new Species(std::string("Out")+rank+sv, 0, 0, world));
      out->set_populate_range(range);
      species_list.push_back(out);
      Species* outg(new Species(std::string("OutGhost")+rank+sv, 0, 0, world));
      outg->set_populate_range(range);
      species_list.push_back(outg);
      Species* vacant(new Species(std::string("Vacant")+rank+sv, 0, 0, world));
      vacant->set_populate_range(range);
      species_list.push_back(vacant);
      Species* ghost(new Species(std::string("Ghost")+rank+sv, 0, 0, world));
      ghost->set_populate_range(range);
      species_list.push_back(ghost);
      */
    }
  }

  world.initialize();
  world.run(duration, verbose);
  return 0;
}

