// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id$

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelGenericSparsePixel.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace eutelescope;


EUTelGenericSparsePixel::EUTelGenericSparsePixel()\
  : _xCoord(0),
    _yCoord(0),
    _signal(0)
{
  _xCoord = 0;
  _yCoord = 0;
  _signal = 0;
  _noOfElements = 3;
  _type = kEUTelGenericSparsePixel;
}

EUTelGenericSparsePixel::EUTelGenericSparsePixel(short xCoord, short yCoord, short signal)
  : _xCoord(0),
    _yCoord(0),
    _signal(0)
{
  _xCoord = xCoord;
  _yCoord = yCoord;
  _signal = signal;
  _noOfElements = 3;
  _type = kEUTelGenericSparsePixel;
}


unsigned int EUTelGenericSparsePixel::getNoOfElements() const {
  return _noOfElements;
}

SparsePixelType EUTelGenericSparsePixel::getSparsePixelType() const {
  return _type;
}

void EUTelGenericSparsePixel::print(std::ostream& os) const {
  int bigWidth = 50;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
  os << endl;
  int width = 20;
  os << setw(width) << setiosflags(ios::left) << "Type: "     << _type << endl
     << setw(width) << "Elements: " << _noOfElements << endl
     << setw(width) << "x coord: "  << _xCoord << endl
     << setw(width) << "y coord: "  << _yCoord << endl
     << setw(width) << "signal: "  << _signal << endl;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
}
