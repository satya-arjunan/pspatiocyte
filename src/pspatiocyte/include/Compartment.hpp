#ifndef __COMPARTMENT_HPP
#define __COMPARTMENT_HPP

#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "Common.hpp"
#include "Molecule.hpp"
#include "Species.hpp"
#include "Lattice.hpp"
#include "Reaction.hpp"
#include "ParallelEnvironment.hpp"
#include "Vector.hpp"
#include <boost/random.hpp>
#include <limits>
#include <random>

using namespace std;
using namespace boost;

class Compartment {
public:
  Compartment(string name, COMPARTMENT_TYPE type, uint32_t seed,
              const int proc_size, const int proc_id, const int invalid_id,
              const int vacant_id, const int ghost_id):
    is_parallel_(proc_size > 1),
    invalid_id_(invalid_id),
    vacant_id_(vacant_id),
    ghost_id_(ghost_id),
    name_(name),
    type_(type),
    rng_(rd_()),
    gen_(seed),
    totalPropensity_(0),
    reactionTime_(0),
    prev_react_time_(0),
    next_react_time_(std::numeric_limits<double>::infinity()),
    mol_id_min_(std::numeric_limits<unsigned>::max()/proc_size*proc_id),
    mol_id_max_(mol_id_min_+std::numeric_limits<unsigned>::max()/proc_size),
    mol_id_(mol_id_min_) {
      randdbl_ = new boost::variate_generator<
        boost::mt19937&,boost::uniform_real<>>(gen_, uniform_real<>());
      adj_rand_ = new boost::variate_generator<
        boost::mt19937&, boost::uniform_int<>>(gen_, dist_);
  }

  ~Compartment() {
    delete randdbl_;
    delete adj_rand_;
  }

  string getName() {
    return name_;
  }

  COMPARTMENT_TYPE getType() {
    return type_;
  }

  int getNumberOfVoxel() {
    return numberOfVoxel_;
  }

  int getID() {
    return compartmentID_;
  }

  void setID(int ID) {
    compartmentID_ = ID;
  }

  void set_out_id(const int out_id) {
    out_id_ = out_id;
  }
  
  unsigned get_new_mol_id() {
    ++mol_id_;
    if (mol_id_ > mol_id_max_) {
      mol_id_ = mol_id_min_;
    }
    return mol_id_;
  }

  void initialize(Lattice &g, ParallelEnvironment &pe,
                  std::vector<Species*>& species_list);
  double getOldTotalPropensity();

  void throwinMolecules(Species& s, const unsigned size, Lattice &g);
  /*
   *  methods related to independent reaction
   */
  void attachIndependentReaction(Reaction& r);

  bool eliminateReactant(Species &s, Lattice &g, bool binary, int &center,
                         int &neighbor);

  void generateProduct(Species &s, Lattice &g, int coord);

  Molecule& identifyMolecule(Species &s, Lattice &g);

  void calculatePropensity(Reaction &r, ParallelEnvironment &pe);

  double getNewPropensity(double);

  double getNextTime(ParallelEnvironment &pe, double current_time);
  double getNewNextTime(ParallelEnvironment &pe, double current_time);
  /*
   *  methods related to influenced reaction
   */
  void attachInfluencedReaction(Reaction& r);

  void calculateProbability(Reaction &r, Lattice &g);

  void findMaxProbability(Species &s);

  void calculateCollisionTime(Species &s, ParallelEnvironment &pe);

  double executeDirectMethodReaction(Lattice &g, ParallelEnvironment &pe,
                                     double current_time);

  void walk(std::vector<Species*>& species_list,
            Lattice &g, ParallelEnvironment &pe);

  void walk_on_ghost(std::vector<unsigned>& species_ids,
                     std::vector<unsigned>& src_coords,
                     std::vector<unsigned>& tar_coords,
                     Lattice& g, ParallelEnvironment& pe);
  void throwinOneMolecule(Species &s, Lattice &g);

  bool processingReaction(Reaction &r, Lattice &g, ParallelEnvironment &pe);
  
  bool searchAllNeighbors(Lattice &g, const int &coord, int &center,
                          int &neighbor);

  void calculatePropensityNoParallel(Reaction &r);
  void calculateLocalPropensity(Reaction &r, double);

  // independent reaction time
  void setReactionTime(double time) {
      reactionTime_ = time;
  }
  // independent reaction time
  double getReactionTime() {
      return reactionTime_;
  }

  // print coodinates of molecues
  void addCoordinatesSpecies(Species& species);
  void addNumberSpecies(Species& species);
  void outputCoordinatesHeader(Lattice& lattice, ParallelEnvironment& pe);
  void outputCoordinates(const double current_time);
  void outputNumbersHeader(ParallelEnvironment& pe);
  void outputNumbers(const double current_time);

  void react(Lattice& g, Species& s, const unsigned tarID, Voxel& src_voxel,
             Voxel& tar_voxel, const unsigned src_coord,
             const unsigned tar_coord);
  void react_on_ghost(Lattice& g, Species& s, const unsigned tarID,
                      Voxel& src_voxel, Voxel& tar_voxel,
                      const unsigned src_coord, const unsigned tar_coord, 
                      std::vector<SpillMolecule>& spill_coords);
  void do_reaction(Lattice& g, Reaction& r, Voxel& v0, Voxel& v1,
                   const unsigned c0, const unsigned c1);
  void do_reaction_on_ghost(Lattice& g, Reaction& r, Voxel& v0, Voxel& v1,
                            const unsigned c0, const unsigned c1,
                            std::vector<SpillMolecule>& spill_coords);
  void add_molecule(const int species_id, const unsigned coord,
                    const unsigned mol_id, Voxel& voxel);
  unsigned remove_molecule(Lattice& g, const int species_id, Voxel& voxel);
  void populate_outmolecule_in_voxel(const int out_id,
                                     const int species_id,
                                     const unsigned mol_index,
                                     const unsigned coord, Voxel& voxel);
  void vacate_outmolecule_in_voxel(Lattice& g, const int out_id, Voxel& voxel);

  unsigned remove_molecule_from_molecules(Lattice& g, const int species_id,
                                          const unsigned mol_index);
  void add_molecule_on_ghost(const int species_id, const unsigned coord,
                             Voxel& v,
                             std::vector<SpillMolecule>& spill_coords);
  void remove_molecule_on_ghost(Lattice& g, const int species_id,
                                const unsigned coord, Voxel& v,
                                std::vector<SpillMolecule>& spill_coords);
  void check_voxels(Lattice& g, const double id);
  unsigned check_species_size(Lattice& g, ParallelEnvironment& pe,
                              const unsigned sid, const double id);
  void spill_molecules(Lattice& g, const std::vector<SpillMolecule>& coord,
                      const Vector<unsigned>& ghost);
  void jumpin_molecules(Lattice& g, ParallelEnvironment& pe);


private:
  const int is_parallel_;
  const int invalid_id_;
  const int vacant_id_;
  const int ghost_id_;
  unsigned Nx_;
  unsigned Ny_;
  unsigned Nz_;
  unsigned n_out_voxels_ = 0;
  unsigned n_out_ghost_voxels_ = 0;
  unsigned n_ghost_voxels_ = 0;
  unsigned out_cnts_[10];
  unsigned out_ghost_cnts_[10];
  int cntES = 0;
  int out_id_ = vacant_id_; //default value
  std::random_device rd_;
  std::mt19937 rng_;
  std::vector<Species*> species_list_;
  std::vector<OutMolecule> outmolecules_[10];
  std::vector<std::vector<Molecule>> species_molecules_;
  string name_;                        // name of compartment
  COMPARTMENT_TYPE type_;              // volume or surface
  int compartmentID_;                  // compartment ID
  int numberOfVoxel_;                  // number of voxels
//  long long numberOfVoxel_;            // number of voxels
  vector<int> voxelVector_;            // linear coordinates of voxels
  vector<Species> output_coord_species_;
  vector<Species> output_number_species_;
  vector<Reaction*> independent_reactions_;  // list of reactions
  vector<Reaction*> influenced_reactions_;    // list of reactions

  double volume_;
  double local_volume_;
  double next_react_time_;
  double prev_react_time_;

  double totalPropensity_;
  double globalTotalPropensity_;

  // Independent reation time
  double reactionTime_;
  boost::uniform_int<> dist_ = boost::uniform_int<>(0, 11);
  boost::mt19937 gen_;
  boost::variate_generator<boost::mt19937&, boost::uniform_int<>>* adj_rand_;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<>>* randdbl_;

  void writeMolecule(vector<Molecule> &mv);
  std::vector<std::vector<unsigned>> is_reactive_;
  std::vector<std::vector<unsigned>> influenced_reaction_ids_;
  std::vector<std::vector<double>> reaction_probabilities_;
  Vector<unsigned> ghosts_[8];
  Vector<int> mid_span_;
  Vector<int> begin_[8];
  Vector<int> end_[8];
  Vector<int> bit_[8];
  const unsigned mol_id_min_ = 0;
  const unsigned mol_id_max_ = UINT_MAX;
  unsigned mol_id_;
};

#endif /* __COMPARTMENT_HPP */

