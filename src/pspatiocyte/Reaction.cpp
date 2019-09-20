#include "Reaction.hpp"
#include "World.hpp"

void Reaction::diagnostics() 
{
    fout << "Reaction:" << endl;
    fout << "   name = " <<  name_ << endl;
    fout << "   type = " <<  type_ << endl;
    fout << "   k    = " <<  k_    << endl;


    fout << "   scheme [ ";

    if(s0_) fout          << s0_->getName(); 
    if(s1_) fout << " + " << s1_->getName();
            fout << " => ";
    if(s2_) fout          << s2_->getName();
    if(s3_) fout << " + " << s3_->getName();
            fout << " ] " << endl;
/*
    if(s0_) fout << "   s0   = " <<  s0_->getName() << endl;
    if(s1_) fout << "   s1   = " <<  s1_->getName() << endl;
    if(s2_) fout << "   s2   = " <<  s2_->getName() << endl;
    if(s3_) fout << "   s3   = " <<  s3_->getName() << endl;
*/
    fout << "   ID   = " <<  reactionID_    << endl;
//    fout << "   dt   = " <<  dt_   << endl;
}


//A + B -> C
Reaction::Reaction(string name, Species& A, Species& B, const double k,
                   Species& C, World& world):
    name_(name),
    type_(INFLUENCED),
    r0_(1),
    s0_(&A),
    r1_(1),
    s1_(&B),
    p2_(1),
    s2_(&C),
    p3_(0),
    s3_(NULL),
    k_(k),
    reactionID_(world.add_influenced_reaction(this)) {}

  //A + B  -> C + D
Reaction::Reaction(string name, Species& A, Species& B, const double k,
                   Species& C, Species& D, World& world):
    name_(name),
    type_(INFLUENCED),
    r0_(1),
    s0_(&A),
    r1_(1),
    s1_(&B),
    p2_(1),
    s2_(&C),
    p3_(1),
    s3_(&D),
    k_(k),
    reactionID_(world.add_influenced_reaction(this)) {}

//A -> C
Reaction::Reaction(string name, Species& A, const double k, Species& C,
                   World& world):
    name_(name),
    type_(INDEPENDENT),
    r0_(1),
    s0_(&A),
    r1_(0),
    s1_(NULL),
    p2_(1),
    s2_(&C),
    p3_(0),
    s3_(NULL),
    k_(k),
    reactionID_(world.add_independent_reaction(this)) {}

//A -> C + D
Reaction::Reaction(string name, Species& A, const double k, Species& C,
                   Species& D, World& world):
    name_(name),
    type_(INDEPENDENT),
    r0_(1),
    s0_(&A),
    r1_(0),
    s1_(NULL),
    p2_(1),
    s2_(&C),
    p3_(1),
    s3_(&D),
    k_(k),
    reactionID_(world.add_independent_reaction(this)) {}
