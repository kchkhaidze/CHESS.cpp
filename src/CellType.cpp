/*
    A spartial constraint tumour growth model
    Copyright (C) 2018 Timon Heide (timon.heide@icr.ac.uk)
                     & Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef _DEBUG_
#define D(x) x
#else
#define D(x)
#endif

#include "extern_global_variables.hpp"
#include "CellType.hpp"
#include "Cell.hpp"
#include <iostream>
#include <boost/random/uniform_int_distribution.hpp>

// CellType ///////////////////////////////////////////////////////////////////

// Statics:
unsigned int CellType::msNextId = 1;

// Constructors:
CellType::CellType () // Default definition of no death and aggression.
  : mId(0),
    mBirthRate(1.0),
    mAggression(0.0),
    mMu(1.0),
    mPushPower(INT_MAX),
    mNumMembers(0),
    mpMembers(0)
  {
    mColor[0] = 255; mColor[1] = 255; mColor[2] = 255;
  }

CellType::CellType (int clone_id,
                    double birthrate, double deathrate, double aggression,
                    unsigned int push_power, double mu) // Full exp. def.
  : mId(clone_id),
    mBirthRate(birthrate),
    mDeathRate(deathrate),
    mAggression(aggression),
    mMu(mu),
    mPushPower(push_power),
    mNumMembers(0),
    mpMembers(0)
    {
      mColor[0] = 255; mColor[1] = 255; mColor[2] = 255;
    }


CellType::CellType(int clone_id,
                   double birthrate, double deathrate, double aggression,
                   unsigned int push_power, double mu, unsigned char red,
                   unsigned char green, unsigned char blue)
    : mId(clone_id),
      mBirthRate(birthrate),
      mDeathRate(deathrate),
      mAggression(aggression),
      mMu(mu),
      mPushPower(push_power),
      mNumMembers(0),
      mpMembers(0)
    {
      mColor[0] = red; mColor[1] = green; mColor[2] = blue;
    }


// Destructor:
CellType::~CellType() {
  D(std::cout << "CellType::~CellType " << this << std::endl;)
  while(mpMembers.size() != 0) {
    Cell* pCell = mpMembers.back();
    delete pCell;
      D(std::cout << "Deleting cell: " << pCell << std::endl;)
  }
  msNextId = 0;
}


// Getter functions:
unsigned int CellType::Id() const {return mId;}
double CellType::Birthrate() const {return mBirthRate;}
double CellType::Deathrate() const {return mDeathRate;}
double CellType::Aggression() const {return mAggression;};
double CellType::Mu() const {return mMu;};
unsigned int CellType::PushPower() const {return mPushPower;};

unsigned char* CellType::Color() {return mColor;};
unsigned long CellType::NumMembers() const {return mNumMembers;}

Cell* CellType::RandomMember() {
  unsigned long max_i = mNumMembers - 1;
  boost::random::uniform_int_distribution<unsigned long> unif_member(0, max_i);
  unsigned long sel_member = unif_member(rng);

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Sampling ###########" << std::endl;)
  D(std::cout << "  Next reaction member:" << std::endl;)
  D(std::cout << "    Member No: " << mpMembers[sel_member]->Id() << std::endl;)
  D(std::cout << "    Member: " << mpMembers[sel_member] << std::endl;)
  D(std::cout << "    Type: " << this->Id() << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  return mpMembers[sel_member];
}

// Setter functions:
void CellType::RegisterMember(Cell* pCell){
  pCell->TypeIndex(mpMembers.size());
  mpMembers.push_back(pCell); // Append pointer to new cell to member vector
  mNumMembers++;              // Increase the population count
}

void CellType::DeregisterMember(Cell* pCell){
  D(std::cout << "Deregister member " << pCell << std::endl;)
  std::vector<Cell*>::size_type i = pCell->TypeIndex();  // Index of deleted member.
  std::vector<Cell*>::size_type j = mpMembers.size() - 1;  // Index of last member.

  pCell->TypeIndex(-1);

  // Copy last element of vector into the now empty location:
  if (i != j) { // if the deleted element is not the last element.
    mpMembers[j]->TypeIndex(i);
    mpMembers[i] = mpMembers[j];
  }

  mpMembers.pop_back(); // Remove the last element
  mNumMembers--; // Decrease the population count
}

// Print functions:
void CellType::Print() const {
  std::cout << std::endl;
  std::cout << "########## Cell type ##########" << std::endl;
  std::cout << "  ID: " << mId << std::endl;
  std::cout << "  Birth rate: " << mBirthRate << std::endl;
  std::cout << "  Aggression: " << mAggression << std::endl;
  std::cout << "  Members: " << mNumMembers << std::endl;
  std::cout << "###############################" << std::endl;
}

void CellType::PrintAllMembers() const {
  std::cout << std::endl;
  std::cout << "#### Cell type members ########" << std::endl;
  std::cout << "  ID:" << mId << std::endl;
  std::cout << std::endl;
  for(std::vector<Cell*>::size_type i = 0; i != mpMembers.size(); i++) {
    std::cout << "  Member " << i << ": " << std::endl;
    std::cout << "    " << mpMembers[i]->Id() << std::endl;
  }
  std::cout << "###############################" << std::endl;
}
