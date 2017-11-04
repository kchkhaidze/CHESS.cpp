#ifndef CELLTYPE_H
#define CELLTYPE_H

// Forward declerations: ///////////////////////////////////////////////////////
class Phylogeny_Node;
class Cell;
class Universe;

// Includes: ///////////////////////////////////////////////////////////////////
#include <vector>

// CellType ///////////////////////////////////////////////////////////////////

class CellType {
    // Identifiers:
    const unsigned int mId;
    static unsigned int msNextId; // Increments from one.

    // Properties:
    float mBirthRate;
    float mAlpha;
    float mBeta;
    float mAggression;
    unsigned int mPushPower;
    unsigned char mColor[3];

    // Member specific:
    unsigned long mNumMembers;
    std::vector<Cell*> mpMembers;

  public:
    // Constructors:
    CellType();
    CellType(float, float, float, float);
    CellType(float, float, float, float, unsigned char, unsigned char,
             unsigned char);

    // Destructor:

    // Getter functions:
    unsigned int Id();
    float Birthrate();
    float Alpha();
    float Beta();
    float Aggression();
    unsigned int PushPower();
    unsigned char* Color();
    unsigned long NumMembers();
    Cell* RandomMember();

    // Setter functions:
    void RegisterMember(Cell*);
    void DeregisterMember(Cell*);

    // Print functions:
    void Print();
    void PrintAllMembers();
};

#endif // CELLTYPE_H
