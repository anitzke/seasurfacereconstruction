package geometry;

import java.awt.geom.Point2D;
import java.util.Comparator;

import processing.Calculations;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class PointAngleComparator implements Comparator<Point2D> {
  private Point2D p = new Point2D.Double(0, 0);
  public PointAngleComparator() {

  }
  public PointAngleComparator(Point2D p) {
    this.p = p;
  }
  @Override
  public int compare(Point2D p1, Point2D p2) {
    if(p1.equals(p)) {
      return -1;
    }
    if(p2.equals(p)) {
      return 1;
    }
    double a = Calculations.direction(p, p1, p2);
    if(a < 0) {
      return -1;
    }else if(a > 0) {
      return 1;
    }
    double d1 = this.p.distance(p1);
    double d2 = this.p.distance(p2);
    if(d1 < d2) {
      return -1;
    }else if(d1 > d2) {
      return 1;
    }
    return 0;
  }

}
