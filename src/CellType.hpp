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
    int mId;
    static unsigned int msNextId; // Increments from one.

    // Properties:
    double mBirthRate;
    double mDeathRate;
    double mAggression;
    double mMu;
    unsigned int mPushPower;
    unsigned char mColor[3];

    // Member specific:
    unsigned long mNumMembers;
    std::vector<Cell*> mpMembers;

  public:
    // Constructors:
    CellType();
    CellType(int, double, double, double, double, double);
    CellType(int, double, double, double, unsigned int, double);
    CellType(int, double, double, double, double, double,
             unsigned char, unsigned char, unsigned char);
    CellType(int, double, double, double, unsigned int, double,
             unsigned char, unsigned char, unsigned char);

    // Destructor:
    ~CellType();

    // Getter functions:
    unsigned int Id() const;
    double Birthrate() const;
    double Deathrate() const;
    double Aggression() const;
    double Mu() const;
    unsigned int PushPower() const;
    unsigned char* Color();
    unsigned long NumMembers() const;
    Cell* RandomMember();

    // Setter functions:
    void RegisterMember(Cell*);
    void DeregisterMember(Cell*);

    // Print functions:
    void Print() const;
    void PrintAllMembers() const;
};

#endif // CELLTYPE_H
