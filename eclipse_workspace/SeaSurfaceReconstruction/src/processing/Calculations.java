package processing;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

import shapes3D.Point3D;
import triangulation.Triangle;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Calculations {

  private static double tolerance = Math.pow(10, -9);
  
  public static double[] planeParameters(ElevationPoint[] points) {
    Point2D[] p = new Point2D[points.length];
    double[] h = new double[points.length];
    for(int i = 0; i < points.length; i++) {
      p[i] = points[i].p;
      h[i] = points[i].h;
    }
    return planeParameters(p, h);
  }
  
  public static double[] planeParameters(Point2D[] polygon, double[] heights) {
    
    int p1 = 0;
    int p2 = 1;
    int p3 = 2;
    
    for(int i = 2; i < polygon.length; i++) {
      if(pointToLine(polygon[i], polygon[p1], polygon[p2]) != 0) {
        p3 = i;
        break;
      }
    }
    
    double alpha = (polygon[p2].getY()-polygon[p1].getY())*(heights[p3]-heights[p1])-(polygon[p3].getY()-polygon[p1].getY())*(heights[p2]-heights[p1]);
    double beta = (heights[p2]-heights[p1])*(polygon[p3].getX()-polygon[p1].getX())-(heights[p3]-heights[p1])*(polygon[p2].getX()-polygon[p1].getX());
    double gamma = (polygon[p2].getX()-polygon[p1].getX())*(polygon[p3].getY()-polygon[p1].getY())-(polygon[p3].getX()-polygon[p1].getX())*(polygon[p2].getY()-polygon[p1].getY());
    double delta = -(alpha*polygon[p1].getX()+beta*polygon[p1].getY()+gamma*heights[p1]);
    return new double[]{alpha, beta, gamma, delta};
    
  }

  public static int factorial(int i) {

    if(i <= 1) {
      return 1;
    }

    return i * factorial(i-1);
  }

  /**
   * returns the minimum value of x and y of the given grid and the step size of the grid
   * @param grid
   * @return [min x, min y, step x, step y]
   */
  public static double[] getGridData(Point3D[][] grid) {
    double[] data = {grid[0][0].getX(), 
        grid[0][0].getY(), 
        grid[1][0].getX()-grid[0][0].getX(), 
        grid[0][1].getY()-grid[0][0].getY()};
    return data;
  }

  /**
   * returns the minimum value of x and y of the given grid and the step size of the grid
   * @param grid
   * @return [min x, min y, step x, step y]
   */
  public static double[] getGridData(Point2D[][] grid) {
    double[] data = {grid[0][0].getX(), 
        grid[0][0].getY(), 
        grid[1][0].getX()-grid[0][0].getX(), 
        grid[0][1].getY()-grid[0][0].getY()};
    return data;
  }

  public static Point3D ellip2kart(Point3D ellP, double a, double b) {
    double phi = ellP.getY()*Math.PI/180;
    double lambda = ellP.getX()*Math.PI/180;
    double h = ellP.getZ()*100;

    double e2 = (a*a-b*b)/(a*a);
    double N = a / Math.sqrt(1-e2*Math.sin(phi)*Math.sin(phi));

    double x = (N+h)*Math.cos(phi)*Math.cos(lambda)/10000;
    double y = (N+h)*Math.cos(phi)*Math.sin(lambda)/10000;
    double z = ((1-e2)*N+h)*Math.sin(phi)/10000;
    ellP.setLocation(x, y, z);

    return ellP;
  }

  /**
   * 
   * @param p point to be tested
   * @param l1 first point on line
   * @param l2 second point on line
   * @return 1 if p is left of l1-l2; -1 if p is right of l1-l2; 0 else for collinear
   */
  public static int pointToLine(Point2D p, Point2D l1, Point2D l2) {
    double ldx = l2.getX() - l1.getX();
    double ldy = l2.getY() - l1.getY();
    double pdx = p.getX() - l1.getX();
    double pdy = p.getY() - l1.getY();

    double det = ldx * pdy - ldy * pdx;

    if (det > tolerance) return 1;
    if (det < -tolerance) return -1;
    return 0;
  }

  /**
   * checks whether a list of points represents a polygon with clockwise direction or not
   * @param points
   * @return true if polygon is clockwise
   */
  public static boolean isClockwise(List<Point2D> points) {
    double sum = 0;
    for(int i = 0; i < points.size()-1; i++) {
      sum += (points.get(i+1).getX()-points.get(i).getX()) * (points.get(i+1).getY()+points.get(i).getY());
    }
    sum += (points.get(0).getX()-points.get(points.size()-1).getX()) * (points.get(0).getY()+points.get(points.size()-1).getY());
    return (sum > 0);
  }	

  public static double direction(Point2D a, Point2D b, Point2D c) {
    return (b.getY()-a.getY())*(c.getX()-b.getX()) -
        (b.getX()-a.getX())*(c.getY()-b.getY());
    //		return (b.getX()-a.getX())*(c.getY()-a.getY())-(c.getX()-a.getX())*(b.getY()-a.getY());
  }

  public static Point2D getPolyCenter(Point2D[] poly) {
    double xSum = 0;
    double ySum = 0;
    for(int i = 0; i < poly.length; i++) {
      xSum += poly[i].getX();
      ySum += poly[i].getY();
    }

    return new Point2D.Double(xSum/poly.length, ySum/poly.length);
  }


  /**
   * returns the inner angle of a point setup A-B-C
   * 	A	  C
   * 	 \ /
   * 	  B
   * @param A
   * @param B
   * @param C
   * @return
   */
  public static double getInnerAngle(Point2D A, Point2D B, Point2D C) {
    //calculate slopes
    double beta1 = Math.atan2(A.getX()-B.getX(), A.getY()-B.getY());
    double beta2 = Math.atan2(C.getX()-B.getX(), C.getY()-B.getY());

    //correct slopes
    if(beta1 < 0) beta1 += Math.PI*2;
    if(beta2 < 0) beta2 += Math.PI*2;

    //calculate inner angle
    double beta = beta1-beta2;

    //correct inner angle
    if(beta > Math.PI) beta -= Math.PI*2;
    if(beta < -Math.PI) beta += Math.PI*2;
    return Math.abs(beta);
  }
  
  public static double getAngleClockwise(Point2D A, Point2D B, Point2D C) {
    // calculate arctan2
    double beta1 = getDirectionAngle(B, A);
    double beta2 = getDirectionAngle(B, C);

    // calculate inner angle
    double beta = beta1-beta2;
    
    if(beta > Math.PI) beta -= Math.PI*2;
    if(beta < -Math.PI) beta += Math.PI*2;

    return Math.abs(beta);
  }
  
  /**
   * returns the inner spheric angle of a point setup A-B-C
   *  A   C
   *   \ /
   *    B
   * @param A
   * @param B
   * @param C
   * @return
   */
  public static double getSphericAngle(Point2D A, Point2D B, Point2D C) {
    double c = Calculations.distanceGeogrCoord(A, B);
    double a = Calculations.distanceGeogrCoord(B, C);
    double b = Calculations.distanceGeogrCoord(A, C);
    
    double denom = Math.sin(a) * Math.sin(c);
    double counter = Math.cos(b) - Math.cos(a) * Math.cos(c);
    
    double beta = Math.acos(counter/denom);
    return beta;
  }
  
  public static double getDirectionAngle(Point2D B, Point2D A) {
    // calculate arctan2
    double beta = Math.atan2(A.getX()-B.getX(), A.getY()-B.getY());

    // correct angle by quadrant
    if(beta < 0) beta += Math.PI*2;
    
    return beta;
  }

  public static double getInnerAngles(ArrayList<Point2D> corners) {
    if(corners.size() <= 2) {
      return 0;
    }
    double innerAngle = getInnerAngle(corners.get(corners.size()-1), corners.get(0), corners.get(1));
    for(int i = 1; i < corners.size()-1; i++) {
      innerAngle += getInnerAngle(corners.get(i-1), corners.get(i), corners.get(i+1));
    } 
    innerAngle += getInnerAngle(corners.get(0), corners.get(corners.size()-1), corners.get(corners.size()-2));
    return innerAngle;
  }


  public static boolean pointInTriangle3D(Point3D p, Triangle t) {

    return false;
  }
  
  public static double distanceGeogrCoord(Point2D A, Point2D B) {
    double lat = ((A.getY() + B.getY()) / 2) * (Math.PI/180);
    double dx = 111300 * Math.cos(lat) * (A.getX() - B.getX());
    double dy = 111300 * (A.getY() - B.getY());
    double dist = Math.sqrt(dx*dx + dy*dy);
    
    return dist;
  }

  public static double rounding(double value, double precision) {
    @SuppressWarnings("deprecation")
    Double val = new Double(value);
    if(val.isNaN()) {
      return value;
    }
    return Math.round(value*precision)/precision;
  }


}
