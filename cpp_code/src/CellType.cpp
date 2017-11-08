//#define DEBUG
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#include "CellType.hpp"
#include "Cell.hpp"
#include <random>
#include <iostream>

// CellType ///////////////////////////////////////////////////////////////////

// Statics:
unsigned int CellType::msNextId = 1;

// Constructors:
CellType::CellType () // Default definition of no death and aggression.
  : mId(msNextId++),
    mBirthRate(1.0),
    mAlpha(0.0),
    mBeta(1.0),
    mAggression(0.0),
    mMu(1.0),
    mPushPower(INT_MAX),
    mNumMembers(0),
    mpMembers(0)
  {
    mColor[0] = 255; mColor[1] = 255; mColor[2] = 255;
  }

CellType::CellType (float birthrate, float alpha, float beta, float aggression,
                    float mu) // Full exp. def.
  : mId(msNextId++),
    mBirthRate(birthrate),
    mAlpha(alpha),
    mBeta(beta),
    mAggression(aggression),
    mMu(mu),
    mPushPower(INT_MAX),
    mNumMembers(0),
    mpMembers(0)
    {
      mColor[0] = 255; mColor[1] = 255; mColor[2] = 255;
    }

CellType::CellType(float birthrate, float alpha, float beta, float aggression,
                   float mu, unsigned char red, unsigned char blue,
                   unsigned char green)
    : mId(msNextId++),
      mBirthRate(birthrate),
      mAlpha(alpha),
      mBeta(beta),
      mAggression(aggression),
      mMu(mu),
      mPushPower(INT_MAX),
      mNumMembers(0),
      mpMembers(0)
    {
      mColor[0] = red; mColor[1] = green; mColor[2] = blue;

    }

// Destructor:

// Getter functions:
unsigned int CellType::Id(){return mId;}
float CellType::Birthrate(){return mBirthRate;}
float CellType::Alpha(){return mAlpha;};
float CellType::Beta(){return mBeta;};
float CellType::Aggression(){return mAggression;};
float CellType::Mu(){return mMu;};
unsigned int CellType::PushPower(){return mPushPower;};

unsigned char* CellType::Color(){return mColor;};
unsigned long CellType::NumMembers() {return mNumMembers;}

Cell* CellType::RandomMember() {
  //:ToDo: Replace with  proper sampling of index
  unsigned long rm = rand() % mNumMembers;

  // Debug messages:
  D(std::cout << std::endl;)
  D(std::cout << "########## Sampling ###########" << std::endl;)
  D(std::cout << "  Next reaction member:" << std::endl;)
  D(std::cout << "    Member No: " << mpMembers[rm]->Id() << std::endl;)
  D(std::cout << "    Member: " << mpMembers[rm] << std::endl;)
  D(std::cout << "    Type: " << this->Id() << std::endl;)
  D(std::cout << "###############################" << std::endl;)

  return mpMembers[rm];
}

// Setter functions:
void CellType::RegisterMember(Cell* pCell){
  mpMembers.push_back(pCell); // Append pointer to new cell to member vector
  mNumMembers++;              // Increase the population count
}

void CellType::DeregisterMember(Cell* pCell){
  D(std::cout << "Deregister member " << pCell << std::endl;)

  std::vector<Cell*>::iterator it;  // iterator for members vector.

  for (it = mpMembers.begin(); it != mpMembers.end(); ++it) { // iter over vector
    if (*it == pCell) {     // if found nth element
      mpMembers.erase(it);  // remove
      break;                // and break loop.
    }
  }
  mNumMembers--; // Decrease the population count
}

// Print functions:
void CellType::Print(){
  std::cout << std::endl;
  std::cout << "########## Cell type ##########" << std::endl;
  std::cout << "  ID: " << mId << std::endl;
  std::cout << "  Birth rate: " << mBirthRate << std::endl;
  std::cout << "  Alpha: " << mAlpha << std::endl;
  std::cout << "  Beta: " << mBeta << std::endl;
  std::cout << "  Aggression: " << mAggression << std::endl;
  std::cout << "  Members: " << mNumMembers << std::endl;
  std::cout << "###############################" << std::endl;
}

void CellType::PrintAllMembers(){
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
