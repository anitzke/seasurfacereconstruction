package geometry;

import java.awt.geom.Point2D;
import java.util.Comparator;

import processing.ValueComparator;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Point2DComparator implements Comparator<Point2D> {

  private double tolerance = Math.pow(10, -9);

  public Point2DComparator() {

  }

  public Point2DComparator(double t) {
    this.tolerance = t;
  }

  @Override
  public int compare(Point2D p1, Point2D p2) {
    ValueComparator vc = new ValueComparator();
    if(vc.compare(p1.getX(), p2.getX()) == 0) {
      return vc.compare(p1.getY(), p2.getY());
    }
    return vc.compare(p1.getX(), p2.getX());
  }


}
