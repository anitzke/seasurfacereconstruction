package triangulation;

import java.awt.geom.Path2D;
import java.util.ArrayList;
import java.util.List;

import processing.Projection;
import shapes3D.Point3D;


/**
 * container-class for raw input data 
 * 
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public interface RawData {

  public static class Triangulate implements RawData{

    public List<Point3D[]> points = new ArrayList<Point3D[]>();
    public List<String> identifiers = new ArrayList<String>();
    public List<String> epochI = new ArrayList<String>();
    public List<String> epochJ = new ArrayList<String>();
    public List<Point3D[][]> rasters = new ArrayList<Point3D[][]>();
    public Path2D pointBounds = new Path2D.Double();
    public Path2D rasterBounds = new Path2D.Double();
    public TriangleCreator tc = new TriangleCreator.Empty();
    public Projection proj = new Projection();

    public Triangulate() {

    }

    public Triangulate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters,
        TriangleCreator tc) {
      this.points = points;
      this.identifiers = ids;
      this.epochI = idsI;
      this.epochJ = idsJ;
      this.rasters = rasters;
      this.tc = tc;
    }
    
    public Triangulate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI,  List<String> idsJ, List<Point3D[][]> rasters, 
        Path2D pointBounds, TriangleCreator tc) {
      this.points = points;
      this.identifiers = ids;
      this.epochI = idsI;
      this.epochJ = idsJ;
      this.rasters = rasters;
      this.pointBounds = pointBounds;
      this.tc = tc;
    }
    
    public Triangulate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI,  List<String> idsJ, List<Point3D[][]> rasters, 
        Path2D pointBounds, Path2D rasterBounds, TriangleCreator tc) {
      this.points = points;
      this.identifiers = ids;
      this.epochI = idsI;
      this.epochJ = idsJ;
      this.rasters = rasters;
      this.pointBounds = pointBounds;
      this.rasterBounds = rasterBounds;
      this.tc = tc;
    }

    public void setProj(Projection proj) {
      this.proj = proj;
    }
    
  }

  public static class Evaluate extends Triangulate{

    public TriangulationInstance singleEvaluation;

    public Evaluate() {

    }

    public Evaluate(List<Point3D[]> points, List<String> ids, List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, tc);
    }

    public Evaluate(List<Point3D[]> points, List<String> ids, List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, Path2D pointBounds, TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, pointBounds, tc);
    }

    public Evaluate(List<Point3D[]> points, List<String> ids, List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, Path2D pointBounds, Path2D rasterBounds, TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, pointBounds, rasterBounds, tc);
    }

    public void setSingleEvaluation(TriangulationInstance t) {
      this.singleEvaluation = t;
    }
    
    public void setProj(Projection proj) {
      super.setProj(proj);
    }

  }

  public static class Validate implements RawData{

    public List<TriangulationInstance> triangulations;
    public List<String> identifiers;

    public Validate() {

    }

    public Validate(List<TriangulationInstance> triangulations, List<String> ids) {
      this.triangulations = triangulations;
      this.identifiers = ids;
    }

  }

}
