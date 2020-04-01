package geometry;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeSet;

import processing.Calculations;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class ConvexHull {
  public Point2D[] points;
  public List<Point2D> convexHull;
  
  
  public ConvexHull() {
  }
  
  public ConvexHull(Point2D[] points) {
    this.points = points;
    this.convexHull = computeConvexHull(points);
  }

  public Point2D[] getPoints() {
    return points;
  }

  public List<Point2D> getConvexHull() {
    return convexHull;
  }

  public void setPoints(Point2D[] points) {
    this.points = points;
  }

  public void setConvexHull(List<Point2D> convexHull) {
    this.convexHull = convexHull;
  }

  public List<Point2D> computeConvexHull(Point2D[] points) {

    Point2D leastPoint = new Point2D.Double(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);

    // get least (smallest point)
    for(Point2D p:points) {
      if(p.getY() < leastPoint.getY()) {
        leastPoint = p;
      } else if(p.getY() == leastPoint.getY()) {
        if(p.getX() < leastPoint.getX()) {
          leastPoint = p;
        }
      }
    }

    PointAngleComparator pca = new PointAngleComparator(leastPoint);
    TreeSet<Point2D> angleSortedPointsSet = new TreeSet<Point2D>(pca);
    angleSortedPointsSet.addAll(Arrays.asList(points));

    List<Point2D> angleSortedPoints = new ArrayList<Point2D>(angleSortedPointsSet);

    int m = angleSortedPoints.size();
    int j = 2;
    while(j < m) {
      // find collinear points from least point and remove in angleSortedPoints
      while(Calculations.direction(angleSortedPoints.get(0), angleSortedPoints.get(j), angleSortedPoints.get(j-1)) == 0) {
        angleSortedPoints.remove(j-1);
        m--;
        if(j >= m) {
          break;
        }
      }
      j++;
    }
    angleSortedPoints.size();
    Stack<Point2D> stack = new Stack<Point2D>();

    stack.push(angleSortedPoints.get(0)); 
    stack.push(angleSortedPoints.get(1)); 
    stack.push(angleSortedPoints.get(2));

    // search for points where the orientation of the triangle is counterclockwise
    for(int i = 3; i < angleSortedPoints.size(); i++) {
      while(Calculations.direction(nextToPeek(stack), stack.peek(), angleSortedPoints.get(i)) > 0) {
        stack.pop(); 
      }
      stack.push(angleSortedPoints.get(i));
    }
    stack.push(angleSortedPoints.get(0));
    stack.size();
    return stack;
  }
  
  public boolean[] computeConvexHull(Point2D[] points, Map<Point2D, Integer> pointMap) {
    boolean[] hull = new boolean[points.length];

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
    angleSortedPointsSet.addAll(Arrays.asList(points));

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

    hull[pointMap.get(angleSortedPoints.get(0))] = true; 
    hull[pointMap.get(angleSortedPoints.get(1))] = true; 
    hull[pointMap.get(angleSortedPoints.get(2))] = true;

    Stack<Integer> stack = new Stack<Integer>();

    stack.push(0); stack.push(1); stack.push(2);

    for(int i = 3; i < angleSortedPoints.size(); i++) {
      while(Calculations.direction(angleSortedPoints.get(nextToPeekInteger(stack)), angleSortedPoints.get(stack.peek()), angleSortedPoints.get(i)) > 0) {
        hull[pointMap.get(angleSortedPoints.get(stack.pop()))] = false; 
      }
      stack.push(i);
      hull[pointMap.get(angleSortedPoints.get(i))] = true; 
    }

    List<Point2D> polygon = new ArrayList<Point2D>();
    while(!stack.isEmpty()) {
      polygon.add(angleSortedPoints.get(stack.pop()));
    }
    polygon.add(polygon.get(0));

    return hull;
  }
  
  public List<HashableLine> compConvexHullEdges() {
    List<Point2D> ps = new ArrayList<Point2D>(this.convexHull);
    List<HashableLine> edges = new ArrayList<HashableLine>();
    
    for(int i = 0; i < (ps.size() - 1); i++) {
      HashableLine hl = new HashableLine(ps.get(i), ps.get(i+1));
      checkTurnEdge(hl);
      edges.add(hl);
    }
    
    return edges;
  }
  
  private Point2D nextToPeek(Stack<Point2D> stack) {
    Point2D i = stack.pop();
    Point2D j = stack.peek();
    stack.push(i);
    return j;   
  }

  private int nextToPeekInteger(Stack<Integer> stack) {
    int i = stack.pop();
    int j = stack.peek();
    stack.push(i);
    return j;   
  }
  
  private boolean checkTurnEdge(HashableLine line) {
    if(line.getP1().getX() < line.getP2().getX() || line.getP1().getX() == line.getP2().getX() && line.getP1().getY() < line.getP2().getY()) {
      return false;
    }
    line.turn();
    return true;
  }

  public double area() {
    double convexArea = 0;
    for(int i = 0; i < this.convexHull.size()-1; i++) {
      Point2D p1 = this.convexHull.get(i);
      Point2D p2 = this.convexHull.get(i+1);
      convexArea += p1.getX()*p2.getY()-p1.getY()*p2.getX();
    }
    convexArea = Math.abs(convexArea/2);
    
    return convexArea;
  }

}
