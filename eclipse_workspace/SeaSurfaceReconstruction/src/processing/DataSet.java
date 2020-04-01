package processing;

import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.index.kdtree.KdNode;
import com.vividsolutions.jts.index.kdtree.KdTree;
import com.vividsolutions.jts.index.strtree.STRtree;

import geometry.Point2DClockwiseComparator;
import shapes3D.Point3D;
import triangulation.Triangle;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class DataSet {

  public static List<Point3D> generateCirclePoints(Point3D center, int n, double radius, double[] noiseLimit) {
    List<Point3D> points = new ArrayList<Point3D>();

    double radiantSteps = Math.PI*2/n;

    for(double d = 0; d < Math.PI*2; d += radiantSteps) {
      double nd = noiseItUp(d, noiseLimit[0]);
      double nr = noiseItUp(radius, noiseLimit[1]);
      double nz = noiseItUp(center.getZ(), noiseLimit[2]);
      double x = center.getX() + Math.sin(nd) * nr;
      double y = center.getY() + Math.cos(nd) * nr;
      double z = center.getZ() + nz;
      points.add(new Point3D.Double(x, y, z));
    }

    return points;
  }

  public static double noiseItUp(double value, double noiseLimit) {
    Random r = new Random();
    double rV = -noiseLimit + (noiseLimit - -noiseLimit) * r.nextDouble();
    return value + rV;
  }

  public static ArrayList<Path2D> triangle2polygon(List<Triangle> triangles){
    ArrayList<Path2D> polygons = new ArrayList<Path2D>();
    for(Triangle t:triangles) {
      ArrayList<Point2D> ps = new ArrayList<Point2D>();
      ps.add(t.getA());
      ps.add(t.getB());
      ps.add(t.getC());
      polygons.add(createPolygon(ps));
    }
    return polygons;
  }

  public static Path2D createPolygon(ArrayList<Point2D> points) {
    Path2D path = new Path2D.Double();
    path.moveTo(points.get(0).getX(), points.get(0).getY());
    for(int i = 1; i < points.size(); i++) {
      path.lineTo(points.get(i).getX(), points.get(i).getY());
    }
    path.closePath();		
    return path;
  }

  public static Point3D[] getPointsFromGrid(double[] PoI, double[] gridData, Point3D[][] grid) {
    Point3D[] points = new Point3D[PoI.length/2];
    for(int i = 0; i < PoI.length; i++) {
      double x = PoI[i];
      double y = PoI[++i];
      int indexX = (int)((x - gridData[0]) / gridData[2]);
      int indexY = (int)((y - gridData[1]) / gridData[3]);
      points[(i-1)/2] = new Point3D.Double(x, y, grid[indexX][indexY].getZ());
    }
    return points;
  }

  /**
   * computes every empty triangle formed by given points
   * @param points
   * @return empty triangles
   */
  public static List<Triangle> getEmptyTriangles(Point2D[] points) {

    STRtree pointTree = new STRtree();
    for (Point2D p : points) {
      //	    	System.out.println(p);
      Envelope env = new Envelope(p.getX(), p.getX(), p.getY(), p.getY());
      pointTree.insert(env , p);
    }

    STRtree triangleTree = new STRtree();
    List<Triangle> emptyTriangles = new ArrayList<Triangle>();

    for (int i = 0; i < points.length; i++) {
      for (int j = i + 1; j < points.length; j++) {
        for (int k = j + 1; k < points.length; k++) {
          Triangle t = new Triangle(points[i], points[j], points[k]);

          boolean isEmpty = true;
          for(Object o : pointTree.query(t.getEnvelope())) {
            Point2D p = (Point2D) o;
            if (t.contains(p)) {
              isEmpty = false;
              break;
            }
          }
          if (isEmpty && t.area() > 0) {
            triangleTree.insert(t.getEnvelope(), t);
            emptyTriangles.add(t);
          }
        }
      }  
    }
    System.out.println("There are " + emptyTriangles.size() + " empty triangles.");
    return emptyTriangles;
  }

  /**
   * computes every empty triangle formed by given points
   * !not working reliably!
   * @param points
   * @return empty triangles
   */
  public static List<Triangle> getEmptyTrianglesSweep(Point2D[] points) {
    return getEmptyTrianglesSweep(points, null);
  }

  /**
   * computes every empty triangle formed by given points
   * !not working reliably!
   * @param points
   * @return empty triangles
   */
  public static List<Triangle> getEmptyTrianglesSweep(Point2D[] points, double[] heights) {

    HashMap<Point2D, Double> heightMap = null;
    boolean dim3 = false;
    if(heights != null) {
      dim3 = true;
      heightMap = new HashMap<Point2D, Double>();
      for(int i = 0; i < points.length; i++) {
        heightMap.put(points[i], heights[i]);
      }
    }    

    List<Triangle> emptyTriangles = new ArrayList<Triangle>();

    for(int i = 0; i < points.length; i++) {
      Point2D p = points[i];
      
      // array with points in P\p
      Point2D[] A = new Point2D[points.length-1];
      for(int k = 0, o = 0; k < points.length; k++) {
        if(k == i) {
          continue;
        }
        A[o++] = points[k];
      }

      Arrays.sort(A, new Point2DClockwiseComparator(p)); // sorting clockwise with respect to p
      
      for(int j = 0; j < A.length; j++) {
        Point2D q = A[j];
        double alphaMin = Math.PI;

        for(int k = 1; k < A.length; k++) {

          Point2D r = A[(j + k) % A.length];          

          if(Calculations.pointToLine(r, p, q) == 1) {
            break;
          }
            
          double alpha = Calculations.getInnerAngle(p, q, r);
          if(alpha <= alphaMin && alpha > 0 && Calculations.pointToLine(r, p, q)!= 0) {
            Triangle t;
            if(dim3) {
              t = new Triangle(p, heights[i], q, heightMap.get(q), r, heightMap.get(r));
            } else {
              t = new Triangle(p, q, r);
            }

            emptyTriangles.add(t);
            alphaMin = alpha;
          }
        }
      }
    }   
    Set<Triangle> set = new HashSet<>(emptyTriangles);
    
    return new ArrayList<Triangle>(set);
  }

  public static List<Triangle> getEmptyTriangles3D(Point2D[] points, double[] heights) {

    STRtree pointTree = new STRtree();
    for (Point2D p : points) {
      //	    	System.out.println(p);
      Envelope env = new Envelope(p.getX(), p.getX(), p.getY(), p.getY());
      pointTree.insert(env , p);
    }

    STRtree triangleTree = new STRtree();
    List<Triangle> emptyTriangles = new ArrayList<Triangle>();

    for (int i = 0; i < points.length; i++) {
      for (int j = i + 1; j < points.length; j++) {
        for (int k = j + 1; k < points.length; k++) {
          Triangle t = new Triangle(points[i], heights[i], points[j], heights[j], points[k], heights[k]);

          boolean isEmpty = true;
          for (Object o : pointTree.query(t.getEnvelope())) {
            Point2D p = (Point2D) o;
            if (t.contains(p)) {
              isEmpty = false;
              break;
            }
          }

          if (isEmpty && t.area() > 0) {
            triangleTree.insert(t.getEnvelope(), t);
            emptyTriangles.add(t);
          }
        }
      }  
    }
    //	    System.out.println("There are " + emptyTriangles.size() + " empty triangles.");
    return emptyTriangles;
  }

  public static List<Triangle> getEmptyTrianglesMapped(Point2D[] points, double[] heights, List<int[]> triangleMap) {

    STRtree pointTree = new STRtree();
    for (Point2D p : points) {
      Envelope env = new Envelope(p.getX(), p.getX(), p.getY(), p.getY());
      pointTree.insert(env , p);
    }

    STRtree triangleTree = new STRtree();
    List<Triangle> emptyTriangles = new ArrayList<Triangle>();

    for (int i = 0; i < points.length; i++) {
      for (int j = i + 1; j < points.length; j++) {
        for (int k = j + 1; k < points.length; k++) {
          Triangle t = new Triangle(points[i], heights[i], points[j], heights[j], points[k], heights[k]);

          boolean isEmpty = true;
          for (Object o : pointTree.query(t.getEnvelope())) {
            Point2D p = (Point2D) o;
            if (t.contains(p)) {
              isEmpty = false;
              break;
            }
          }

          if (isEmpty && t.area() > 0) {
            triangleTree.insert(t.getEnvelope(), t);
            emptyTriangles.add(t);
            triangleMap.add(new int[] {i, j, k});
          }
        }
      }  
    }
    //	    System.out.println("There are " + emptyTriangles.size() + " empty triangles.");
    return emptyTriangles;
  }

  public static List<Triangle> getAllTrianglesExceptions(Point2D[] points, double[] heights, List<Point2D> exceptions){
    HashMap<Point2D, Integer> idxMap = new HashMap<Point2D, Integer>();
    for(int i = 0; i < points.length; i++) {
      idxMap.put(points[i], i);
    }

    boolean[] exc = new boolean[points.length];
    for(Point2D p : exceptions) {
      exc[idxMap.get(p)] = true;
    }

    STRtree pointTree = new STRtree();
    for (Point2D p : points) {
      Envelope env = new Envelope(p.getX(), p.getX(), p.getY(), p.getY());
      pointTree.insert(env , p);
    }

    //STRtree triangleTree = new STRtree();
    List<Triangle> emptyTriangles = new ArrayList<Triangle>();

    for (int i = 0; i < points.length; i++) {
      for (int j = i + 1; j < points.length; j++) {
        for (int k = j + 1; k < points.length; k++) {
          Triangle t = new Triangle(points[i], heights[i], points[j], heights[j], points[k], heights[k]);     

          boolean isEmpty = true;

          for (Object o : pointTree.query(t.getEnvelope())) {
            Point2D p = (Point2D) o;
            if(t.contains(p) && !exc[idxMap.get(p)]) {
              isEmpty = false;
              break;
            }
          }

          if (isEmpty && t.area() > 0) {
            emptyTriangles.add(t);
          }
        }
      }  
    }
    //	    System.out.println("There are " + emptyTriangles.size() + " empty triangles.");
    return emptyTriangles;
  }

  public static List<Triangle> getAllTriangles(Point2D[] points, double[] heights){	      
    //STRtree triangleTree = new STRtree();
    List<Triangle> allTriangles = new ArrayList<Triangle>();

    for (int i = 0; i < points.length; i++) {
      for (int j = i + 1; j < points.length; j++) {
        for (int k = j + 1; k < points.length; k++) {
          Triangle t = new Triangle(points[i], heights[i], points[j], heights[j], points[k], heights[k]);       
          if(t.area() > 0) allTriangles.add(t);
        }
      }  
    }
    //	    System.out.println("There are " + emptyTriangles.size() + " empty triangles.");
    return allTriangles;
  }

  /**
   * gets k-order Delaunay triangles (that are empty)
   * @param points
   * @param heights
   * @param numPoints
   * @return
   */
  public static List<Triangle> getOrderKTriangles(Point2D[] points, double[] heights, 
      List<Triangle> candidates, int numPoints) {

    KdTree tree = new KdTree();
    for (Point2D p : points) {
      tree.insert(new Coordinate(p.getX(), p.getY()));
    }
    
    List<Triangle> triangles = new ArrayList<Triangle>();
    
    if(candidates.isEmpty()) {

      for (int i = 0; i < points.length; i++) {
        for (int j = i + 1; j < points.length; j++) {
          for (int k = j + 1; k < points.length; k++) {
            // sets triangle vertices in 2D with counterclockwise ordering
            Triangle t = new Triangle(points[i], heights[i], points[j], heights[j], points[k], heights[k]);;

            Point2D center = t.circumCenter();
            double radius = center.distance(t.getA());

            boolean isEmpty = true;            
            int count = 0;
            for (Object o : tree.query(t.envCircle(center, radius))) {
              KdNode node = (KdNode) o;
              Point2D p = new Point2D.Double(node.getCoordinate().x, 
                  node.getCoordinate().y);

              // if the requested point is one of the triangle vertices
              if(p.equals(t.getA()) || p.equals(t.getB()) || p.equals(t.getC())) {
                continue;
              }

              // if p is in the circumcircle count up 
              if(center.distance(p) < radius) {
                count++;
              } else {
                continue;
              }

              if(count > numPoints || t.contains(p)) {
                isEmpty = false;
                break;
              }
            }

            if(isEmpty && t.area() > 0) {
              triangles.add(t);
            }
          }
        }  
      }
    } else {
      for(Triangle t : candidates) {
        Point2D center = t.circumCenter();
        double radius = center.distance(t.getA());

        boolean isEmpty = true;            
        int count = 0;
        for (Object o : tree.query(t.envCircle(center, radius))) {
          KdNode node = (KdNode) o;
          Point2D p = new Point2D.Double(node.getCoordinate().x, 
              node.getCoordinate().y);

          // if the requested point is one of the triangle vertices
          if(p.equals(t.getA()) || p.equals(t.getB()) || p.equals(t.getC())) {
            continue;
          }

          // if p is in the circumcircle count up 
          if(center.distance(p) < radius) {
            count++;
          } else {
            continue;
          }

          if(count > numPoints || t.contains(p)) {
            isEmpty = false;
            break;
          }
        }

        if(isEmpty && t.area() > 0) {
          triangles.add(t);
        }
      }
    }  
    
    
    return triangles;
  } 

  public static Point2D[][] getGridPoints(double[] x, double[] y){
    double minX = x[0];
    double maxX = x[1];
    int nX = (int)x[2];
    double minY = y[0];
    double maxY = y[1];
    int nY = (int)y[2];

    double rangeX = maxX-minX;
    double deltaX = rangeX/(nX-1);

    double rangeY = maxY-minY;
    double deltaY = rangeY/(nY-1);

    Point2D[][] points = new Point2D[nX][nY];

    for(int i = 0; i < nX; i++) {
      for(int j = 0; j < nY; j++) {
        points[i][j] = new Point2D.Double(minX + i*deltaX, minY + j*deltaY);
      }
    }

    return points;
  }

  /**
   * 
   * @param x [range minimum in x, range maximum in x, number of points in x]
   * @param y [range minimum in y, range maximum in y, number of points in y]
   * @param noise to randomize points a bit
   * @return Point2D[] List of Points
   */
  public static Point2D[] getNoisyGridPoints(double[] x, double[] y, double noise) {
    double minX = x[0];
    double maxX = x[1];
    int nX = (int)x[2];
    double minY = y[0];
    double maxY = y[1];
    int nY = (int)y[2];

    double rangeX = maxX-minX;
    double deltaX = rangeX/(nX-1);

    double rangeY = maxY-minY;
    double deltaY = rangeY/(nY-1);

    Point2D[] points = new Point2D[nX*nY];

    int counter = 0;
    for(int i = 0; i < nX; i++) {
      for(int j = 0; j < nY; j++) {
        double noiseX = Math.random() * noise;
        double noiseY = Math.random() * noise;
        double px;
        if(Math.random() > 0.5) {
          px = minX + i*deltaX + noiseX;
        }else {
          px = minX + i*deltaX - noiseX; 
        }

        double py;
        if(Math.random() > 0.5) {
          py = minY + j*deltaY + noiseY;
        }else {
          py = minY + j*deltaY - noiseY; 
        }
        points[counter] = new Point2D.Double(px, py);
        counter++;
      }
    }
    return points;
  }

  public static double[][] getGridHeights(Point2D[][] points, double[] heightData){
    double[][] gridHeights = new double[points.length][points[0].length];

    double deltaH = (heightData[1]-heightData[0])/(points.length-1);

    for(int i = 0; i < gridHeights.length; i++) {
      for(int j = 0; j < gridHeights[0].length; j++) {
        //				gridHeights[i][j] = (heightData[1]+heightData[0])/2+Math.sin(points[i][j].getX())*((heightData[1]-heightData[0])/2);
        gridHeights[i][j] = heightData[0] + i * deltaH;
      } 
    }

    return gridHeights;
  }

  public static Point2D[] getRandomPoints(int n, double[] dataX, double[] dataY) {
    Point2D[] points = new Point2D[n];
    for (int i = 0; i < n; i++) {
      points[i] = new Point2D.Double(dataX[0] + (Math.random()*(dataX[1]-dataX[0])), dataY[0] + (Math.random()*(dataY[1]-dataY[0])));
    }
    return points;
  }

  public static double[] getRandomHeights(Point2D[] points, double[] heightData) {
    double[] heights = new double[points.length];
    //		double deltaH = (heightData[1]-heightData[0])/(points.length-1);
    double range = heightData[1]-heightData[0];
    for(int i = 0; i < heights.length; i++) {
      heights[i] = (heightData[1]+heightData[0])/2+Math.sin(points[i].getX())*(range/2);
      //			heights[i] = heightData[0] + points[i].getX() * deltaH;
    }
    return heights;
  }

  public static Point3D[] getMeasurements(Point3D[][] data, String monthYear) {

    String[] parts = monthYear.split(", ");

    int year = (Integer.parseInt(parts[1]) - 1993) * 12;
    if(year < 0 || year > 264) {
      System.out.println("invalid input");
      return null;
    }

    int month;

    if(parts[0].equals("January") || parts[0].equals("january")) {
      month = 0;
    }else if(parts[0].equals("February") || parts[0].equals("february")) {
      month = 1;
    }else if(parts[0].equals("March") || parts[0].equals("march")) {
      month = 2;
    }else if(parts[0].equals("April") || parts[0].equals("april")) {
      month = 3;
    }else if(parts[0].equals("May") || parts[0].equals("may")) {
      month = 4;
    }else if(parts[0].equals("June") || parts[0].equals("june")) {
      month = 5;
    }else if(parts[0].equals("July") || parts[0].equals("july")) {
      month = 6;
    }else if(parts[0].equals("August") || parts[0].equals("august")) {
      month = 7;
    }else if(parts[0].equals("September") || parts[0].equals("september")) {
      month = 8;
    }else if(parts[0].equals("October") || parts[0].equals("october")) {
      month = 9;
    }else if(parts[0].equals("November") || parts[0].equals("november")) {
      month = 10;
    }else if(parts[0].equals("December") || parts[0].equals("december")) {
      month = 11;
    }else {
      System.out.println("invalid input");
      return null;
    }

    return data[year + month];
  }

  public static void sortOutList(List ps, int[] ix) {
    for(int i = 0; i < ix.length; i++) {
      ps.set(ix[i], null);
    }
    for(int i = 0; i < ps.size(); i++) {
      if(ps.get(i) == null) {
        ps.remove(i);
        i--;
      }
    }
  }

}
