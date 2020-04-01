package processing;

import java.awt.geom.Path2D;
import java.util.ArrayList;
import java.util.List;

import shapes3D.Point3D;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public interface PointSelector {

  /**
   * selects points based on a criteria
   * NoNaN: returns points with a numeric height value
   * Bounds: returns points with a numeric height value inside a polygon
   * BoundsNaN: returns all points inside a polygon
   * All: returns all points
   * AllNaN: converts all height values to NaN
   * Specific: converts the height value of specific points to NaN
   * @return
   */
  public Point3D[] select();

  public void setPoints(Point3D[] points);

  public static class NoNaN implements PointSelector{

    private Point3D[] allPoints;

    public NoNaN() {

    }

    public NoNaN(Point3D[] points) {
      allPoints = points;
    }

    @Override
    public void setPoints(Point3D[] points) {
      allPoints = points;
    }

    @Override
    public Point3D[] select() {
      List<Point3D> selectedPoints = new ArrayList<Point3D>();
      for(Point3D p : allPoints) {
        if(!Double.isNaN(p.getZ())) {
          selectedPoints.add(p);
        }
      }
      return selectedPoints.toArray(new Point3D[selectedPoints.size()]);
    }

  }

  public static class Bounds implements PointSelector{

    private Point3D[] allPoints;
    private Path2D bounds;

    public Bounds() {

    }

    public Bounds(Point3D[] points, Path2D boundingPolygon) {
      this.allPoints = points;
      this.bounds = boundingPolygon;
    }

    @Override
    public void setPoints(Point3D[] points) {
      allPoints = points;
    }

    public void setBounds(Path2D boundingPolygon) {
      this.bounds = boundingPolygon;
    }

    @Override
    public Point3D[] select() {
      List<Point3D> selectedPoints = new ArrayList<Point3D>();
      for(Point3D p : allPoints) {
        if(bounds.contains(p.getX(), p.getY()) && !Double.isNaN(p.getZ())) {
          selectedPoints.add(p);
        }
      }
      return selectedPoints.toArray(new Point3D[selectedPoints.size()]);
    }

  }

  public static class BoundsNaN implements PointSelector{

    private Point3D[] allPoints;
    private Path2D bounds;

    public BoundsNaN() {

    }

    public BoundsNaN(Point3D[] points, Path2D boundingPolygon) {
      this.allPoints = points;
      this.bounds = boundingPolygon;
    }

    @Override
    public void setPoints(Point3D[] points) {
      allPoints = points;
    }

    public void setBounds(Path2D boundingPolygon) {
      this.bounds = boundingPolygon;
    }

    @Override
    public Point3D[] select() {
      List<Point3D> selectedPoints = new ArrayList<Point3D>();
      for(Point3D p : allPoints) {
        if(bounds.contains(p.getX(), p.getY())) {
          selectedPoints.add(p);
        }
      }
      return selectedPoints.toArray(new Point3D[selectedPoints.size()]);
    }

  }

  public static class All implements PointSelector{

    private Point3D[] allPoints;

    public All() {

    }

    public All(Point3D[] points) {
      allPoints = points;
    }

    @Override
    public Point3D[] select() {
      return allPoints;
    }

    @Override
    public void setPoints(Point3D[] points) {
      allPoints = points;
    }

  }

  public static class AllNaN implements PointSelector {
    private Point3D[] allPoints;

    public AllNaN() {

    }

    public AllNaN(Point3D[] points) {
      allPoints = points;
    }

    @Override
    public void setPoints(Point3D[] points) {
      allPoints = points;
    }

    @Override
    public Point3D[] select() {
      for(Point3D p : allPoints) {
        p.setLocation(p.getX(), p.getY(), Double.NaN);
      }
      return allPoints;
    }
  }

  public static class Specific implements PointSelector{

    private Point3D[] allPoints;
    private List<Point3D> variablePoints;

    public Specific() {

    }

    public Specific(Point3D[] points, List<Point3D> variablePoints) {
      this.allPoints = points;
      this.variablePoints = variablePoints;
    }

    @Override
    public void setPoints(Point3D[] points) {
      allPoints = points;
    }

    public void setVariablePoints(List<Point3D> variablePoints) {
      this.variablePoints = variablePoints;
    }

    @Override
    public Point3D[] select() {
      List<Point3D> selectedPoints = new ArrayList<Point3D>();
      for(Point3D p : allPoints) {
        if(variablePoints.contains(p)) {
          p.setLocation(p.getX(), p.getY(), Double.NaN);
          selectedPoints.add(p);
        }else {
          if(!Double.isNaN(p.getZ())) {
            selectedPoints.add(p);
          }
        }

      }
      return selectedPoints.toArray(new Point3D[selectedPoints.size()]);
    }

  }

}
