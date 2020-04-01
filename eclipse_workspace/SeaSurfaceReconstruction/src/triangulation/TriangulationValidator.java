package triangulation;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Stack;
import java.util.TreeSet;

import geometry.PointAngleComparator;
import processing.Calculations;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class TriangulationValidator {

  private List<Triangle> triangles;
  private List<Point2D> points;
  private List<String> validationLog;
  private boolean[] tests;

  public TriangulationValidator() {
    triangles = new ArrayList<Triangle>();
    points = new ArrayList<Point2D>();
    validationLog = new ArrayList<String>();
    tests = new boolean[4];
  }

  public TriangulationValidator(List<Triangle> triangles, List<Point2D> points, List<String> valLog) {
    this.triangles = triangles;
    this.points = points;
    this.validationLog = valLog;
    this.tests = new boolean[4];
  }

  public void setTriangles(List<Triangle> triangles) {
    this.triangles = triangles;
  }

  public void addTriangle(Triangle toAdd) {
    triangles.add(toAdd);
  }

  public void addAllTriangles(List<Triangle> triangles) {
    this.triangles.addAll(triangles);
  }

  public void setPoints(List<Point2D> points) {
    this.points = points;
  }

  public void addPoint(Point2D p) {
    points.add(p);
  }

  public void addAllPoints(List<Point2D> points) {
    this.points.addAll(points);
  }

  public void addValidationLine(String line) {
    validationLog.add(line);
  }

  public void addValidationLines(List<String> lines) {
    validationLog.addAll(lines);
  }

  public void setValidationLog(List<String> log) {
    validationLog = log;
  }

  public void clear() {
    triangles.clear();
    points.clear();
    validationLog.clear();
    tests = new boolean[4];
  }

  /**
   * validates the triangulation T stored in it self
   * @return true if T is valid
   */
  public boolean validate() {

    // check number of triangles
    List<Point2D> convexHull = getConvexHull(this.points);
    int nTrianglesExpected = points.size()*2-(convexHull.size()-1)-2;
    int nTriangles = triangles.size();
    boolean nTest  = nTrianglesExpected == nTriangles;

    tests[0] = nTest;

    if(!nTest) {
      validationLog.add("Expected number of triangles " + nTrianglesExpected + ", actual number of triangles: " + nTriangles);
    }

    // check area of triangles
    double convexArea = 0;
    for(int i = 0; i < convexHull.size()-1; i++) {
      Point2D p1 = convexHull.get(i);
      Point2D p2 = convexHull.get(i+1);
      convexArea += p1.getX()*p2.getY()-p1.getY()*p2.getX();
    }
    convexArea = Math.abs(convexArea/2);

    double triangleArea = 0;
    for(Triangle triangle:triangles) {
      triangleArea += triangle.area();
    }

    boolean areaTest = triangleArea-convexArea < Math.pow(10, -5) && triangleArea-convexArea > -Math.pow(10, -5);

    tests[1] = areaTest;

    if(!areaTest) {
      validationLog.add("Area of convexHull: " + convexArea + ", summarized area of triangles: " + triangleArea);
    }

    // check use of points
    HashMap<Point2D, Integer> pointMap= new HashMap<Point2D, Integer>();
    for(int i = 0; i < points.size(); i++) {
      pointMap.put(points.get(i), i);
    }	
    boolean[] usedPoints = new boolean[points.size()];
    for(Triangle triangle:triangles) {
      Point2D A = triangle.getA();
      Point2D B = triangle.getB();
      Point2D C = triangle.getC();
      usedPoints[pointMap.get(A)] = true;
      usedPoints[pointMap.get(B)] = true;
      usedPoints[pointMap.get(C)] = true;
    }
    boolean pointTest = true;
    int nUsedPoints = 0;
    for(int i = 0; i < usedPoints.length; i++) {
      if(usedPoints[i]) {
        nUsedPoints++;
      }else {
        pointTest = false;
      }
    }

    tests[2] = pointTest;

    if(!pointTest) {
      validationLog.add("Number of points: " + points.size() + ", number of used points: " + nUsedPoints);
    }

    // check overlap of triangles
    List<int[]> overlappingTriangles = new ArrayList<int[]>();
    for(int i = 0; i < triangles.size(); i++) {
      for(int j = i; j < triangles.size(); j++) {
        if(triangles.get(i).intersects(triangles.get(j))) {
          overlappingTriangles.add(new int[] {i, j});
        }
      }
    }
    boolean overlapTest = overlappingTriangles.isEmpty();

    tests[3] = overlapTest;

    if(!overlapTest) {
      String temp =  "following triangles intersect: \n";
      for(int[] oT:overlappingTriangles) {
        temp += triangles.get(oT[0]) + ", " + triangles.get(oT[1]) + "\n";
      }
      validationLog.add(temp);
    }

    return nTest && areaTest && pointTest && overlapTest;
  }

  /**
   * validates the given triangulation T
   * @param triangles List of triangles of T 
   * @param points List of points of T
   * @param validationLog List of Strings (filled by this function)
   * @return true if T is valid
   */
  public static boolean validate(List<Triangle> triangles, List<Point2D> points, List<String> validationLog) {

    boolean log = false;
    if(validationLog != null) {
      log = true;
    }

    // check number of triangles
    List<Point2D> convexHull = getConvexHull(points);
    int nTrianglesExpected = points.size()*2-(convexHull.size()-1)-2;
    int nTriangles = triangles.size();
    boolean nTest  = nTrianglesExpected == nTriangles;

    if(!nTest && log) {
      validationLog.add("Expected number of triangles " + nTrianglesExpected + ", actual number of triangles: " + nTriangles);
    }

    // check area of triangles
    double convexArea = 0;
    for(int i = 0; i < convexHull.size()-1; i++) {
      Point2D p1 = convexHull.get(i);
      Point2D p2 = convexHull.get(i+1);
      convexArea += p1.getX()*p2.getY()-p1.getY()*p2.getX();
    }
    convexArea = Math.abs(convexArea/2);

    double triangleArea = 0;
    for(Triangle triangle:triangles) {
      triangleArea += triangle.area();
    }

    boolean areaTest = triangleArea-convexArea < Math.pow(10, -5) && triangleArea-convexArea > -Math.pow(10, -5);

    if(!areaTest && log) {
      validationLog.add("Area of convexHull: " + convexArea + ", summarized area of triangles: " + triangleArea);
    }

    // check use of points
    HashMap<Point2D, Integer> pointMap= new HashMap<Point2D, Integer>();
    for(int i = 0; i < points.size(); i++) {
      pointMap.put(points.get(i), i);
    }	
    boolean[] usedPoints = new boolean[points.size()];
    for(Triangle triangle:triangles) {
      Point2D A = triangle.getA();
      Point2D B = triangle.getB();
      Point2D C = triangle.getC();
      usedPoints[pointMap.get(A)] = true;
      usedPoints[pointMap.get(B)] = true;
      usedPoints[pointMap.get(C)] = true;
    }
    boolean pointTest = true;
    int nUsedPoints = 0;
    for(int i = 0; i < usedPoints.length; i++) {
      if(usedPoints[i]) {
        nUsedPoints++;
      }else {
        pointTest = false;
      }
    }

    if(!pointTest && log) {
      validationLog.add("Number of points: " + points.size() + ", number of used points: " + nUsedPoints);
    }

    // check overlap of triangles
    List<int[]> overlappingTriangles = new ArrayList<int[]>();
    for(int i = 0; i < triangles.size(); i++) {
      for(int j = i; j < triangles.size(); j++) {
        if(triangles.get(i).intersects(triangles.get(j))) {
          overlappingTriangles.add(new int[] {i, j});
        }
      }
    }
    boolean overlapTest = overlappingTriangles.isEmpty();

    if(!overlapTest && log) {
      String temp =  "following triangles intersect: \n";
      for(int[] oT:overlappingTriangles) {
        temp += triangles.get(oT[0]) + ", " + triangles.get(oT[1]) + "\n";
      }
      validationLog.add(temp);
    }

    return nTest && areaTest && pointTest && overlapTest;
  }

  public static List<Point2D> getConvexHull(List<Point2D> points) {

    Point2D leastPoint = new Point2D.Double(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);

    for(Point2D p:points) {
      if(p.getY() < leastPoint.getY()) {
        leastPoint = p;
      }
      if(p.getY() == leastPoint.getY()) {
        if(p.getX() < leastPoint.getX()) {
          leastPoint = p;
        }
      }
    }

    PointAngleComparator pca = new PointAngleComparator(leastPoint);
    TreeSet<Point2D> angleSortedPointsSet = new TreeSet<Point2D>(pca);
    angleSortedPointsSet.addAll(points);

    List<Point2D> angleSortedPoints = new ArrayList<Point2D>(angleSortedPointsSet);

    int m = angleSortedPoints.size();
    int j = 2;
    while(j < m) {
      while(Calculations.direction(angleSortedPoints.get(0), angleSortedPoints.get(j), angleSortedPoints.get(j-1)) == 0) {
        angleSortedPoints.remove(j-1);
        m--;
        if(j >= m) {
          break;
        }
      }
      j++;
    }

    Stack<Point2D> stack = new Stack<Point2D>();

    stack.push(angleSortedPoints.get(0)); 
    stack.push(angleSortedPoints.get(1)); 
    stack.push(angleSortedPoints.get(2));

    for(int i = 3; i < angleSortedPoints.size(); i++) {
      while(Calculations.direction(nextToPeek(stack), stack.peek(), angleSortedPoints.get(i)) > 0) {
        stack.pop(); 
      }
      stack.push(angleSortedPoints.get(i));
    }
    stack.push(angleSortedPoints.get(0));

    return stack;
  }

  private static Point2D nextToPeek(Stack<Point2D> stack) {
    Point2D i = stack.pop();
    Point2D j = stack.peek();
    stack.push(i);
    return j;		
  }

}
