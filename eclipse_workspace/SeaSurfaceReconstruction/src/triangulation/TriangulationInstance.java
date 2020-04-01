package triangulation;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import com.vividsolutions.jts.geom.Envelope;

import geometry.HashableLine;
import processing.Projection;
import shapes3D.Point3D;

/**
 * container-class for triangulation data 
 * 
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class TriangulationInstance implements Cloneable{

  private String name;
  private Point2D[] points;
  private double[] heights;
  private Point2D[][] gridPoints;
  private double[][] gridHeights;
  private List<Triangle> triangles;
  private List<Point2D> convexHull;
  private List<Triangle> candidates;
  private HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap;
  

  public TriangulationInstance() {
    name = "";
    points = new  Point2D[0];
    heights = new double[0];
    gridPoints = new Point2D[0][0];
    gridHeights = new double[0][0];
    triangles = new ArrayList<Triangle>();
    convexHull = new ArrayList<Point2D>();
    candidates = new ArrayList<Triangle>();
    edgeMap = new HashMap<HashableLine, ArrayList<LinkedList<Integer>>>();
  }

  public TriangulationInstance(Point2D[] ps, double[] hs, Point2D[][] gP, double[][] gH, List<Triangle> ts) {
    points = ps;
    heights = hs;
    gridPoints = gP;
    gridHeights = gH;
    triangles = ts;
    convexHull = new ArrayList<Point2D>();
    candidates = new ArrayList<Triangle>();
    edgeMap = new HashMap<HashableLine, ArrayList<LinkedList<Integer>>>();
  }
  
  public TriangulationInstance(Point2D[] ps, double[] hs, Point2D[][] gP, double[][] gH, List<Triangle> ts, 
		  List<Point2D> convexHull) {
	    points = ps;
	    heights = hs;
	    gridPoints = gP;
	    gridHeights = gH;
	    triangles = ts;
	    this.convexHull = convexHull;
	  }

  public TriangulationInstance(String name, Point2D[] points, double[] heights,
      Point2D[][] gridPoints, double[][] gridHeights, List<Triangle> triangles,
      List<Point2D> convexHull, List<Triangle> candidates,
      HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap) {
    this.name = name;
    this.points = points;
    this.heights = heights;
    this.gridPoints = gridPoints;
    this.gridHeights = gridHeights;
    this.triangles = triangles;
    this.convexHull = convexHull;
    this.candidates = candidates;
    this.edgeMap = edgeMap;
  }
 
  public TriangulationInstance(Point2D[] ps, double[] hs, List<Triangle> ts) {
    points = ps;
    heights = hs;
    gridPoints = null;
    gridHeights = null;
    triangles = ts;
  }

  public TriangulationInstance(Point2D[] ps, List<Triangle> ts) {
    points = ps;
    heights = null;
    gridPoints = null;
    gridHeights = null;
    triangles = ts;
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public Point2D[] getPoints() {
    return points;
  }

  public double[] getHeights() {
    return heights;
  }

  public Point2D[][] getGridPoints() {
    return gridPoints;
  }

  public double[][] getGridHeights() {
    return gridHeights;
  }

  public List<Triangle> getTriangles() {
    return triangles;
  }

  public List<Point2D> getConvexHull() {
  return convexHull;
  }
  
  public List<Triangle> getCandidates() {
    return candidates;
  }

  public HashMap<HashableLine, ArrayList<LinkedList<Integer>>> getEdgeMap() {
    return edgeMap;
  }

  public void setPoints(Point2D[] points) {
    this.points = points;
  }

  public void setHeights(double[] heights) {
    this.heights = heights;
  }
  
  public void setHeights(Point3D[] points) {
    double[] h = new double[points.length];
    for (int i=0; i<points.length; i++) {
      h[i] = points[i].getZ();
    }
    this.heights = h;
  }

  public void setGridPoints(Point2D[][] gridPoints) {
    this.gridPoints = gridPoints;
  }

  public void setGridHeights(double[][] gridHeights) {
    this.gridHeights = gridHeights;
  }

  public void setTriangles(List<Triangle> triangles) {
    this.triangles = triangles;
  }
  
  public void setGrid(Point2D[][] gridPoints, double[][] gridHeights) {
    this.gridPoints = gridPoints;
    this.gridHeights = gridHeights;
  }
  
  public void setPoints(Point2D[] points, double[] heights) {
    this.points = points;
    this.heights = heights;
  }
  
  public void setPoints(Point3D[] points) {
    Point2D[] p = new Point2D[points.length];
    double[] h = new double[points.length];
    for (int i=0; i<points.length; i++) {
      p[i] = new Point2D.Double(points[i].getX(), points[i].getY());
      h[i] = points[i].getZ();
    }
    this.points = p;
    this.heights = h;
  }

  public void setConvexHull(List<Point2D> convexHull) {
	this.convexHull = convexHull;
  }
  
  public void setCandidates(List<Triangle> candidates) {
    this.candidates = candidates;
  }

  public void setEdgeMap(
      HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap) {
    this.edgeMap = edgeMap;
  }
  
  

  public Envelope getEnvelope() {
    Envelope env = new Envelope();

    for(Triangle t : triangles) {
      env.expandToInclude(t.getEnvelope());
    }

    return env;
  }
  
  public void project(Projection proj){    
    // stations
    this.points = proj.project(this.points);
    
    // grid
    this.gridPoints = proj.project(this.gridPoints);
    
    // convexHull
    Point2D[] ch = new Point2D[convexHull.size()];
    convexHull.toArray(ch);
    ch = proj.project(ch);
    this.convexHull = Arrays.asList(ch);
    
    // triangles of resulted Triangulation
    List<Triangle> newTriangles = new ArrayList<Triangle>();
    for(int i=0; i < this.triangles.size(); i++) {
      newTriangles.add(proj.project(this.triangles.get(i)));
    }
    this.triangles = newTriangles;
    
    // triangle candidates
    List<Triangle> newCandidates = new ArrayList<Triangle>();
    for(int i=0; i < this.candidates.size(); i++) {
      newCandidates.add(proj.project(this.candidates.get(i)));
    }
    this.candidates = newCandidates;
    
  }
  
  public double area() {
    double triangleArea = 0;
    for(Triangle triangle : this.triangles) {
      //ystem.out.println(triangle.area());
      triangleArea += triangle.area();
    }
    
    return triangleArea;
  }
  
  public TriangulationInstance copy() {
    TriangulationInstance newT = new TriangulationInstance(this.name,
        this.points, this.heights, this.gridPoints, this.gridHeights,
        this.triangles, this.convexHull, this.candidates, this.edgeMap);
 
    return newT;
  }

  public Object clone()throws CloneNotSupportedException{  
    return (TriangulationInstance)super.clone();  
 }
}
