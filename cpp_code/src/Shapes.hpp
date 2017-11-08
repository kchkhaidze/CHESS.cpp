#ifndef SHAPE_H
#define SHAPE_H

// abstract class Shape
class Shape {
  public:
    // Other functions:
    virtual bool next_coordinate(int&, int&, int&) =0;
};

class Box: public Shape {
    int mX, mY, mZ;
    int mMinX, mMinY, mMinZ;
    int mMaxX, mMaxY, mMaxZ;

  public:
    // Constructor:
    Box(int, int, int, int, int, int);

    // Other functions:
    bool next_coordinate(int&, int&, int&);
};


#endif // SHAPES_H
