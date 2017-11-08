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
