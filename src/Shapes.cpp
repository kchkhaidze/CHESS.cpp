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

#include "Shapes.hpp"
#include <iostream>

Box::Box(int minx, int miny, int minz, int maxx, int maxy, int maxz)
  : mX(minx), mY(miny), mZ(minz),
    mMinX(minx), mMinY(miny), mMinZ(minz),
    mMaxX(maxx), mMaxY(maxy), mMaxZ(maxz)
  {}

// Other functions:
bool Box::next_coordinate(int& x, int& y, int& z) {

  if (mZ > mMaxZ){
    mZ = mMinZ;
    mY++;
  }

  if (mY > mMaxY){
    mY = mMinY;
    mX++;
  }

  if (mX > mMaxX){
    return false;
  }

  x = mX;
  y = mY;
  z = mZ++;
  return true;
}

void Box::reset() {
  mX = mMinX;
  mY = mMinY;
  mZ = mMinZ;
}


PositionVector::PositionVector() : index(0), elements(0) {};

PositionVector::PositionVector(std::vector <int> x,
                               std::vector <int> y,
                               std::vector <int> z)
{
  if (x.size() == y.size() && x.size() == z.size()) {
    mvX = x;
    mvY = y;
    mvZ = z;
    index = 0;
    elements = x.size();
  } else {
    std::cerr << "Error: Failed to create PositionVector." << std::endl;
    index = 0;
    elements = 0;
  }
}

void PositionVector::insert_index(int x, int y, int z) {
  mvX.push_back(x);
  mvY.push_back(y);
  mvZ.push_back(z);
  elements++;
}

bool PositionVector::next_coordinate(int &x, int&y, int&z) {
  if (index < elements) {
    x = mvX[index];
    y = mvY[index];
    z = mvZ[index];
    index++;
    return(true);
  }
  return false;
}

void PositionVector::reset() {
  index=0;
}
