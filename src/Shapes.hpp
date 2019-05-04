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

#ifndef SHAPE_H
#define SHAPE_H

#include <vector>

// abstract class Shape
class Shape {
  public:
    // Other functions:
    virtual bool next_coordinate(int&, int&, int&) =0;
    virtual void reset() =0; // Function reset the shape object.
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
    void reset(); // Function reset the shape object.
};

class PositionVector: public Shape {
    std::vector <int> mvX, mvY, mvZ;
    unsigned long index, elements;

  public:
    // Constructor:
    PositionVector();
    PositionVector(std::vector <int>, std::vector <int>, std::vector <int>);

    // Other functions:
    void insert_index(int, int, int);
    bool next_coordinate(int&, int&, int&);
    void reset(); // Function reset the shape object.
};

#endif // SHAPES_H
