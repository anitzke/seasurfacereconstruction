package processing;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

import org.gicentre.utils.spatial.Ellipsoid;
import org.gicentre.utils.spatial.LambertConformalConic;

import com.vividsolutions.jts.geom.Envelope;

import shapes3D.Point3D;
import triangulation.RawData;
import triangulation.Triangle;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Projection {
  
  public List<Point3D[]> originPoints;
  public List<Point3D[][]> originRasters;
  public List<Point3D[]> projectedPoints;
  public List<Point3D[][]> projectedRasters;
  public LambertConformalConic LLC;
  public String direction;
  
  
  public Projection() {
  }
  
  
  public Projection(RawData.Triangulate mData) {
    this.originPoints = mData.points;
    this.originRasters = mData.rasters;
    initializeLCC();
  }
  
  public Projection(RawData.Triangulate mData, String direction) {
    this.originPoints = mData.points;
    this.originRasters = mData.rasters;
    initializeLCC();
    this.direction = direction;
  }
  
  
	public Projection(List<Point3D[]> points, List<Point3D[][]> rasters) {
    this.originPoints = points;
    this.originRasters = rasters;
    initializeLCC();
  }
  
  
  public List<Point3D[]> getOriginPoints() {
    return originPoints;
  }


  public List<Point3D[][]> getOriginRasters() {
    return originRasters;
  }


  public List<Point3D[]> getProjectedPoints() {
    return projectedPoints;
  }


  public List<Point3D[][]> getProjectedRasters() {
    return projectedRasters;
  }


  public LambertConformalConic getLLC() {
    return LLC;
  }
  

  public String getDirection() {
    return direction;
  }


  public void setDirection(String direction) {
    this.direction = direction;
  }


  public void setOriginPoints(List<Point3D[]> originPoints) {
    this.originPoints = originPoints;
  }


  public void setOriginRasters(List<Point3D[][]> originRasters) {
    this.originRasters = originRasters;
  }


  public void setProjectedPoints(List<Point3D[]> projectedPoints) {
    this.projectedPoints = projectedPoints;
  }


  public void setProjectedRasters(List<Point3D[][]> projectedRasters) {
    this.projectedRasters = projectedRasters;
  }


  public void setLLC(LambertConformalConic lLC) {
    LLC = lLC;
  }


  public Envelope compEnvelope() {   
    Envelope env = new Envelope();
    for(int i = 0; i < this.originPoints.size(); i++) {
      Point3D[] stations = this.originPoints.get(i);
      for(int j = 0; j < stations.length; j++) {
        env.expandToInclude(stations[j].getX(), stations[j].getY());
      }
    }
    
    return env;
  }
  
  public void initializeLCC() {   
    Envelope env = compEnvelope();
    Ellipsoid ell = new Ellipsoid(Ellipsoid.WGS_84);
    double lat1 = env.getMinY();
    double lat2 = env.getMaxY();
    double lon0 = env.getMinX() + (env.getMaxX() - env.getMinX()) / 2;
    double lat0 = env.getMinY() + (env.getMaxY() - env.getMinY()) / 2;
    double falseEast = 0;
    double falseNorth = 0;
    
    this.LLC = new LambertConformalConic(ell, lat1, lat2, lon0, lat0, falseEast, falseNorth);
  }

  public Point2D project(Point2D point) {
    Point2D new_point = new Point2D.Double();     
    if(direction.equals("toLCC")) {
      new_point = this.LLC.transformCoords(point);
    } else {  // to lon lat
      new_point = this.LLC.invTransformCoords(point);
      return new Point2D.Double(Calculations.rounding(new_point.getX(), 1000),
          Calculations.rounding(new_point.getY(), 1000));
    }       
    return new_point;
  }
  
  public Point2D project(Point3D point) {
   return project(new Point2D.Double(point.getX(), point.getY()));
  }
  
  public Point2D[] project(Point2D[] points) {
    Point2D[] stations = points;
    Point2D[] newStations = new Point2D[points.length];
    for(int j = 0; j < stations.length; j++) {  
      newStations[j] = project(stations[j]);        
    }
    
    return newStations;
  }
  
  public Point3D[] project(Point3D[] points) {
    Point3D[] stations = points;
    Point3D[] newStations = new Point3D[points.length];
    for(int j = 0; j < stations.length; j++) {
      Point2D newP = project(stations[j]);      
      newStations[j] = new Point3D.Double((double) newP.getX(), (double) newP.getY(), stations[j].getZ());        
    }
    
    return newStations;
  }
  
  public Point2D[][] project(Point2D[][] grid) {
    Point2D[][] newGrid = new Point2D[grid.length][grid[0].length];
    for(int j = 0; j < grid.length; j++) {     
      for(int k = 0; k < grid[0].length; k++) {
        Point2D newP = project(grid[j][k]);      
        newGrid[j][k] = new Point2D.Double(newP.getX(), newP.getY());     
      }
    }
    return newGrid;
  }
  
	public List<Point3D[]> project(List<Point3D[]> points) {
		// project TG points
	  List<Point3D[]> projectedPoints = new ArrayList<Point3D[]>();
		for(int i = 0; i < points.size(); i++) {  
		  projectedPoints.add(i, project(points.get(i)));        
      }	
		return projectedPoints;
	}
	
	public List<Point2D[][]> projectRaster2D(List<Point2D[][]> rasters) {
    // project raster
    ArrayList<Point2D[][]> projectedRasters = new ArrayList<Point2D[][]>();
    for(int i = 0; i < rasters.size(); i++) {
      projectedRasters.add(i, project(rasters.get(i)));
    }
    
    return projectedRasters;
	}
	
	public List<Point3D[][]> projectRaster(List<Point3D[][]> rasters) {
	  // project raster
	  ArrayList<Point3D[][]> projectedRasters = new ArrayList<Point3D[][]>();
    for(int i = 0; i < rasters.size(); i++) {
      Point3D[][] grids = rasters.get(i);
      Point3D[][] newGrids = new Point3D[grids.length][grids[0].length];
      for(int j = 0; j < grids.length; j++) {     
        for(int k = 0; k < grids[0].length; k++) {
          Point2D newP = project(grids[j][k]);      
          newGrids[j][k] = new Point3D.Double((double) newP.getX(), (double) newP.getY(), grids[j][k].getZ());        
        }
      }
      projectedRasters.add(i, newGrids);
    }
    
    return projectedRasters;
	}
	
	 public Triangle project(Triangle t) {
	    return new Triangle(project(t.getA()), t.gethA(), 
	        project(t.getB()), t.gethB(), project(t.getC()), t.gethC());
	  }

}
