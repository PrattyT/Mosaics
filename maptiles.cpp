/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include "maptiles.h"
#include <iostream>
#include <map>
//#include "cs225/RGB_HSL.h"

using namespace std;

Point<3> convertToXYZ(LUVAPixel pixel) {
  return Point<3>(pixel.l, pixel.u, pixel.v);
}

MosaicCanvas *mapTiles(SourceImage const &theSource,
                       vector<TileImage> &theTiles) {
  /**
   * @todo Implement this function!
   */

  vector<Point<3>> points;
  map<Point<3>, int> colorToTile;

  for (unsigned i = 0; i < theTiles.size(); i++) {
    Point<3> point = convertToXYZ(theTiles[i].getAverageColor());
    colorToTile[point] = i;
    points.push_back(point);
  }

  KDTree<3> kdTree(points);

  MosaicCanvas *mosaic =
      new MosaicCanvas(theSource.getRows(), theSource.getColumns());

  for (int i = 0; i < mosaic->getRows(); i++) {
    for (int j = 0; j < mosaic->getColumns(); j++) {

      mosaic->setTile(i, j,
                      &theTiles[colorToTile[kdTree.findNearestNeighbor(
                          convertToXYZ(theSource.getRegionColor(i, j)))]]);
    }
  }

  return mosaic;
}
