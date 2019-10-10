//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of pSpatiocyte
//
//        Copyright (C) 2019 Satya N.V. Arjunan
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Motocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Motocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Motocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// and Atsushi Miyauchi
//


#include "Reaction.hpp"
#include "World.hpp"

void Reaction::diagnostics() 
{
    fout << "Reaction:" << endl;
    fout << "   name = " <<  name_ << endl;
    fout << "   type = " <<  type_ << endl;
    fout << "   k    = " <<  k_    << endl;


    fout << "   scheme [ ";

    if(s0_) fout          << s0_->get_name(); 
    if(s1_) fout << " + " << s1_->get_name();
            fout << " => ";
    if(s2_) fout          << s2_->get_name();
    if(s3_) fout << " + " << s3_->get_name();
            fout << " ] " << endl;
/*
    if(s0_) fout << "   s0   = " <<  s0_->get_name() << endl;
    if(s1_) fout << "   s1   = " <<  s1_->get_name() << endl;
    if(s2_) fout << "   s2   = " <<  s2_->get_name() << endl;
    if(s3_) fout << "   s3   = " <<  s3_->get_name() << endl;
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
