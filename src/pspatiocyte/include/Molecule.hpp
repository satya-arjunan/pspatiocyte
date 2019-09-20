#ifndef __MOLECULE_HPP
#define __MOLECULE_HPP

#include "Common.hpp"
using namespace std;

class SpeciesCoord {
public:
  SpeciesCoord(const int s, const unsigned c):
    species_id(s),
    coord(c) {}
  int species_id;
  unsigned coord;
};

class SpillMolecule {
public:
  SpillMolecule(const int s, const unsigned c, const unsigned id):
    species_id(s),
    coord(c),
    mol_id(id) {}
  int species_id;
  unsigned coord;
  unsigned mol_id;
};

class Molecule {
public:
  Molecule(const int c, const unsigned id):
    coord(c),
    mol_id(id) {}
  int coord;
  unsigned mol_id;
};

class OutMolecule {
public:
  OutMolecule(const unsigned c, const unsigned i, const int s):
    coord(c),
    species_id(s),
    mol_index(i) {}
  unsigned coord;
  int species_id;
  unsigned mol_index;
};


#endif /* __MOLECULE_HPP */
