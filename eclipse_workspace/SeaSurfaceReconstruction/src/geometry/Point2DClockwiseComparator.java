package geometry;

import java.awt.geom.Point2D;
import java.util.Comparator;

import processing.ValueComparator;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Point2DClockwiseComparator implements Comparator<Point2D>{
  private Point2D pivot; 

  public Point2DClockwiseComparator() {
    pivot = new Point2D.Double();
  }

  public Point2DClockwiseComparator(Point2D p) {
    this.pivot = p;
  }

  @Override
  public int compare(Point2D p1, Point2D p2) {
    ValueComparator vc = new ValueComparator();
    double a1 = Math.atan2(p1.getX()-pivot.getX(), p1.getY()-pivot.getY());
    if(a1 < 0) a1 += Math.PI*2;
    
    double a2 = Math.atan2(p2.getX()-pivot.getX(), p2.getY()-pivot.getY());
    if(a2 < 0) a2 += Math.PI*2;
      
    return vc.compare(a1, a2);
  }

}
