#ifndef __REACTION_HPP
#define __REACTION_HPP

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include "Common.hpp"

using namespace std;

class Reaction {
public:
  Reaction(const Reaction &s) {
    name_ = s.name_;
    type_ = s.type_;
    r0_ = s.r0_;    s0_ = s.s0_;
    r1_ = s.r1_;    s1_ = s.s1_;
    p2_ = s.p2_;    s2_ = s.s2_;
    p3_ = s.p3_;    s3_ = s.s3_;
    k_ = s.k_;
    dt_ = s.dt_;
    reactionID_ = s.reactionID_;
  }

  //A + B -> C
  Reaction(string name, Species& A, Species& B, const double k, Species& C,
           World& world);
  //A + B  -> C + D
  Reaction(string name, Species& A, Species& B, const double k, Species& C,
           Species& D, World& world);
  //A -> C
  Reaction(string name, Species& A, const double k, Species& C, World& world);
  //A -> C + D
  Reaction(string name, Species& A, const double k, Species& C, Species& D,
           World& world);

  ~Reaction() {}
  string get_name() {
    return name_;
  }
  REACTION_TYPE getType() { 
    return type_; 
  }
  Species* getS0() {
    return s0_;
  }
  Species* getS1() {
    return s1_;
  }
  Species* getS2() {
    return s2_;
  }
  Species* getS3() {
    return s3_;
  }
  int getR0() {
    return r0_;
  }
  int getR1() {
    return r1_;
  }
  int getP2() {
    return p2_;
  }
  int getP3() {
    return p3_;
  }
  double getK() {
    return k_;
  }
  double getDt() {
    return dt_;
  }
  int getID() {
    return reactionID_;
  }
  double get_propensity() {
    return propensity_;
  }
  double get_old_propensity() {
    return old_propensity_;
  }
  void set_propensity(double p) {
    propensity_ = p;
  }
  void set_old_propensity(double p) {
    old_propensity_ = p;
  }
  bool operator==(const Reaction &r) const {
    return name_==r.name_ && type_==r.type_ && reactionID_==r.reactionID_;
  }
  double get_probability() {
    return probability_;
  }
  void set_probability(double p) {
    probability_ = p;
  }
  void diagnostics();
private:
  string name_;         // name of reaction
  REACTION_TYPE type_;  // independent or influenced
  Species* s0_;         // species0 in LHS
  Species* s1_;         // species1 in LHS
  Species* s2_;         // species2 in RHS
  Species* s3_;         // species3 in RHS
  int r0_;              // reactant0 in LHS
  int r1_;              // reactant1 in LHS
  int p2_;              // product2 in RHS
  int p3_;              // product3 in RHS
  int nu0_;             // state change0
  int nu1_;             // state change1
  int nu2_;             // state change2
  int nu3_;             // state change3
  double k_;            // diffusion coefficient
  double dt_;           // reaction time
  int reactionID_;      // ID
  double propensity_;    // propensity
  double old_propensity_;    // propensity
  double probability_;
};

#endif /* __REACTION_HPP */
