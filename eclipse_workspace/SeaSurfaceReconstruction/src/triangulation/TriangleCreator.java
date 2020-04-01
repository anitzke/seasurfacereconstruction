package triangulation;

import java.awt.geom.Point2D;
import java.util.List;

import processing.DataSet;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public interface TriangleCreator {

  public void setPoints(Point2D[] points);

  void setEmptyTriangles(List<Triangle> emptyTriangles);

  public void setHeights(double[] heights);
  
  /**
   * computes every valid Triangle in a set of points depending on instance.
   * Instance Empty: creates every empty triangle in the point set.
   * Instance Constrained: allows a certain amount of point within Triangles
   * Instance All: creates all Triangles in a point set
   * @return List of Triangles
   */
  public List<Triangle> create();

  public class Empty implements TriangleCreator {

    private Point2D[] points;
    private double[] heights;
    private List<Triangle> emptyTriangles; 

    public Empty() {

    }

    public Empty(Point2D[] points, double[] heights) {
      this.points = points;
      this.heights = heights;
    }

    @Override
    public void setPoints(Point2D[] points) {
      this.points = points;
    }

    @Override
    public void setHeights(double[] heights) {
      this.heights = heights;
    }
    
    @Override
    public void setEmptyTriangles(List<Triangle> emptyTriangles) {
      this.emptyTriangles = emptyTriangles;     
    }

    @Override
    public List<Triangle> create() {
//      List<Triangle> tl1 = DataSet.getEmptyTriangles3D(points, heights););
      List<Triangle> tl2 = DataSet.getEmptyTrianglesSweep(points, heights);
      setEmptyTriangles(tl2);
      return tl2;
    }

  }

  public class Constrained implements TriangleCreator {

    private Point2D[] points;
    private double[] heights;
    private int numPoints;
    private List<Triangle> emptyTriangles;

    public Constrained() {

    }

    public Constrained(int n) {
      this.numPoints = n;
    }

    public Constrained(Point2D[] points, double[] heights, int numPoints) {
      this.points = points;
      this.heights = heights;
      this.numPoints = numPoints;
    }
    
    public Constrained(Point2D[] points, double[] heights, int numPoints,
        List<Triangle> emptyTriangles) {
      this.points = points;
      this.heights = heights;
      this.numPoints = numPoints;
      this.emptyTriangles = emptyTriangles;
    }

    @Override
    public void setPoints(Point2D[] points) {
      this.points = points;
    }

    @Override
    public void setHeights(double[] heights) {
      this.heights = heights;
    }
    
    @Override
    public void setEmptyTriangles(List<Triangle> emptyTriangles) {
      this.emptyTriangles = emptyTriangles;     
    }

    public void setN(int n) {
      this.numPoints = n;
    }

    @Override
    public List<Triangle> create() {
      List<Triangle> tl = DataSet.getOrderKTriangles(points, heights, emptyTriangles, numPoints);
      setEmptyTriangles(tl);
      return tl;
    }

  }

  public class All implements TriangleCreator {

    private Point2D[] points;
    private double[] heights;

    public All() {

    }

    public All(Point2D[] points, double[] heights) {
      this.points = points;
      this.heights = heights;
    }

    @Override
    public void setPoints(Point2D[] points) {
      this.points = points;
    }

    @Override
    public void setHeights(double[] heights) {
      this.heights = heights;
    }
    
    @Override
    public void setEmptyTriangles(List<Triangle> emptyTriangles) {    
    }
    
    @Override
    public List<Triangle> create() {
      return DataSet.getAllTriangles(points, heights);
    }

  }


}
