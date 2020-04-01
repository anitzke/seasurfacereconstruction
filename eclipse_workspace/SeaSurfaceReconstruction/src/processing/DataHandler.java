package processing;

import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import shapes3D.Point3D;
import triangulation.TriangulationInstance;
import ucar.netcdf.Netcdf;
import ucar.netcdf.NetcdfFile;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class DataHandler {
  /**
   * computes common points of two triangulations
   * @param t1
   * @param t2
   * @return common points of t1 and t2
   */
  public static Point3D[] getSharedPoints(TriangulationInstance t1, TriangulationInstance t2) {

    HashSet<Point2D> pointSet = new HashSet<Point2D>();
    for(Point2D p : t2.getPoints()) {
      pointSet.add(p);
    }

    List<Point3D> jointPoints = new ArrayList<Point3D>();
    for(int i = 0; i < t1.getPoints().length; i++) {
      if(pointSet.contains(t1.getPoints()[i])) {
        jointPoints.add(new Point3D.Double(t1.getPoints()[i].getX(), 
            t1.getPoints()[i].getY(), t1.getHeights()[i]));
      }
    }

    return jointPoints.toArray(new Point3D[jointPoints.size()]);
  }
  
  
  /**
   * computes common points of two point sets
   * @param p1
   * @param p2
   * @return common points of p1 and p2
   */
  public static Point3D[] getSharedPoints(Point3D[] p1, Point3D[] p2) {

    HashSet<Point2D> pointSet = new HashSet<Point2D>();
    for(Point3D p : p2) {
      pointSet.add(new Point2D.Double(p.getX(), p.getY()));
    }

    List<Point3D> jointPoints = new ArrayList<Point3D>();
    for(int i = 0; i < p1.length; i++) {
      if(pointSet.contains(new Point2D.Double(p1[i].getX(), p1[i].getY()))) {
        jointPoints.add(new Point3D.Double(p1[i].getX(), 
            p1[i].getY(), p1[i].getZ()));
      }
    }

    return jointPoints.toArray(new Point3D[jointPoints.size()]);
  }

  public static String createTimeString(String name, double[] times) {
    return name + ", formulation time: " + times[0] + ", solve time: " + times[1] + " [" +"\u00B5" + "s]";
  }

  public static String createStatString(String name, double[] stats) {
    return name + ", mean: " + stats[0] + ", sum of squared anomalies: "
        + stats[1] + ", variance: " + stats[2] + ", standard deviation: " + stats[3]
        + ", area: " + stats[4] + ", tide gauges: " + stats[5] + ", candidates: " + stats[6];
  }
  
  public static String createRecString(String name, double[] rec_stats) {
    return name + "," + rec_stats[0] + "," + rec_stats[1] + "," + rec_stats[2];
  }

  /**
   * clips a raster to a bounding polygon
   * @param entireRaster
   * @param points
   * @return clipped raster
   */
  public static Point3D[][] clipRaster(Point3D[][] entireRaster, Path2D boundingPolygon) {

    Rectangle2D bb = boundingPolygon.getBounds();

    Point2D[] boundingPoints = rect2points(bb);

    Point3D[][] clippedRaster = clipRaster(entireRaster, boundingPoints);
    for(int i = 0; i < clippedRaster.length; i++) {
      for(int j = 0; j < clippedRaster[0].length; j++) {
        if(!boundingPolygon.contains(clippedRaster[i][j].getX(), clippedRaster[i][j].getY())) {
          clippedRaster[i][j] = new Point3D.Double(clippedRaster[i][j].getX(), clippedRaster[i][j].getY(), Double.NaN);
        }
      }
    }

    return clippedRaster;
  }
  
  public static double[][] addRasters(double[][] raster1, double[][] raster2){
    double[][] addedRaster = new double[raster1.length][raster1[0].length];
    
    if(raster1.length != raster2.length) {
      return null;
    }
    
    for(int i = 0; i < raster1.length; i ++) {
      if(raster1[i].length != raster2[i].length) {
        return null;
      }
      for(int j = 0; j < raster1[i].length; j++) {
        if(!Double.isNaN(raster1[i][j]) && !Double.isNaN(raster2[i][j])) {
          addedRaster[i][j] = raster1[i][j] - raster2[i][j];
        }else if(!Double.isNaN(raster1[i][j])) {
          addedRaster[i][j] = raster1[i][j];
        }else if(!Double.isNaN(raster2[i][j])){
          addedRaster[i][j] = raster2[i][j];
        }else {
          addedRaster[i][j] = Double.NaN;
        }
        
      }
    }
    
    
    return addedRaster;
  }

  public static Point2D[] rect2points(Rectangle2D rect) {
    return new Point2D[] {new Point2D.Double(rect.getMinX(), rect.getMinY()), 
        new Point2D.Double(rect.getMinX()+rect.getWidth(), rect.getMinY()), 
        new Point2D.Double(rect.getMinX()+rect.getWidth(), rect.getMinY()+rect.getHeight()), 
        new Point2D.Double(rect.getMinX(), rect.getMinY()+rect.getHeight())};
  }

  /**
   * clips a raster to the bounding box of given points
   * @param entireRaster
   * @param points
   * @return clipped raster
   */
  public static Point3D[][] clipRaster(Point3D[][] entireRaster, Point2D[] points) {

    double[] gridData = Calculations.getGridData(entireRaster);

    double minX = Double.MAX_VALUE;
    double maxX = Double.MIN_VALUE;
    double minY = Double.MAX_VALUE;
    double maxY = Double.MIN_VALUE;
    for(Point2D p:points) {
      if(p.getX() < minX) {
        minX = p.getX();
      }
      if(p.getX() > maxX) {
        maxX = p.getX();
      }
      if(p.getY() < minY) {
        minY = p.getY();
      }
      if(p.getY() > maxY) {
        maxY = p.getY();
      }
    }

    int	minXindex = (int) Math.floor((minX - gridData[0])/gridData[2]);
    int	minYindex = (int) Math.floor((minY - gridData[1])/gridData[3]);
    int	maxXindex = (int) Math.ceil((maxX - gridData[0])/gridData[2]);
    int	maxYindex = (int) Math.ceil((maxY - gridData[1])/gridData[3]);

    Point3D[][] clippedRaster = new Point3D[maxXindex-minXindex][maxYindex-minYindex];

    int vx = 0;
    int vy = 0;
    for(int ix = minXindex; ix < maxXindex; ix++) {
      for(int iy = minYindex; iy < maxYindex; iy++) {
        //set all heights of points with the fillValue as height (== point is on land) to NaN
        if(entireRaster[ix][iy].getZ() < 100000) { 
          clippedRaster[vx][vy] = entireRaster[ix][iy];
        }else {
          clippedRaster[vx][vy] = new Point3D.Double(entireRaster[ix][iy].getX(), entireRaster[ix][iy].getY(), Double.NaN);
        }
        vy++;
      }
      vy = 0;
      vx++;
    }

    return clippedRaster;
  }

  /**
   * clips a raster to the bounding box of given points
   * @param entireRaster
   * @param points
   * @return clipped raster
   */
  public static Point3D[][] clipRaster(Point3D[][] entireRaster, Point3D[] points) {

    // get min of lon and lat and resolution of grid
    double[] gridData = Calculations.getGridData(entireRaster);

    double minX = Double.MAX_VALUE;
    double maxX = Double.MIN_VALUE;
    double minY = Double.MAX_VALUE;
    double maxY = Double.MIN_VALUE;

    for(Point3D p:points) {
      if(p.getX() < minX) {
        minX = p.getX();
      }
      if(p.getX() > maxX) {
        maxX = p.getX();
      }
      if(p.getY() < minY) {
        minY = p.getY();
      }
      if(p.getY() > maxY) {
        maxY = p.getY();
      }
    }
    // compute indizes of mins and max of points (stations) in raster
    int	minXindex = (int) Math.floor((minX - gridData[0])/gridData[2]);
    int	minYindex = (int) Math.floor((minY - gridData[1])/gridData[3]);
    int	maxXindex = (int) Math.ceil((maxX - gridData[0])/gridData[2]);
    int	maxYindex = (int) Math.ceil((maxY - gridData[1])/gridData[3]);

    Point3D[][] clippedRaster = new Point3D[maxXindex-minXindex][maxYindex-minYindex];

    int vx = 0;
    int vy = 0;
    for(int ix = minXindex; ix < maxXindex; ix++) {
      for(int iy = minYindex; iy < maxYindex; iy++) {
        //set all heights of points with the fillValue as height (== point is on land) to NaN
        clippedRaster[vx][vy] = entireRaster[ix][iy];
        vy++;
      }
      vy = 0;
      vx++;
    }

    return clippedRaster;
  }	

  /**
   * sets the attributes of a empty triangulation
   * @param triangulation
   * @param points
   * @param raster
   * @return
   */
  public static boolean setTriangulation(TriangulationInstance triangulation, Point3D[] points, Point3D[][] raster) {

    if(points != null) {
      Point2D[] points2D = new Point2D[points.length];
      double[] heights = new double[points.length];
      for(int i = 0; i < points.length; i++) {
        points2D[i] = new Point2D.Double(points[i].getX(), points[i].getY());
        heights[i] = points[i].getZ();
      }

      triangulation.setPoints(points2D);
      triangulation.setHeights(heights);
    }

    if(raster != null) {
      Point2D[][] raster2D = new Point2D[raster.length][raster[0].length];
      double[][] rasterHeights = new double[raster.length][raster[0].length];
      for(int i = 0; i < raster.length; i++) {
        for(int j = 0; j < raster[0].length; j++) {
          raster2D[i][j] = new Point2D.Double(raster[i][j].getX(), raster[i][j].getY());
          rasterHeights[i][j] = raster[i][j].getZ();
        }
      }

      triangulation.setGridPoints(raster2D);
      triangulation.setGridHeights(rasterHeights);
    }		

    return true;
  }

  public static TriangulationInstance normalizeTriangulation(TriangulationInstance triangulation) {

    double[] pointHeights = triangulation.getHeights();
    double[][] gridHeights = triangulation.getGridHeights();

    double maxH = Double.MIN_VALUE;
    double minH = Double.MAX_VALUE;

    for(int i = 0; i < pointHeights.length; i++) {
      if(pointHeights[i] > maxH) {
        maxH = pointHeights[i];
      }
      if(pointHeights[i] < minH) {
        minH = pointHeights[i];
      }
    }

    for(int i = 0; i < gridHeights.length; i++) {
      for(int j = 0; j < gridHeights[i].length; j++) {
        if(gridHeights[i][j] > maxH) {
          maxH = gridHeights[i][j];
        }
        if(gridHeights[i][j] < minH) {
          minH = gridHeights[i][j];
        }
      }
    }

    for(int i = 0; i < pointHeights.length; i++) {
      pointHeights[i] = (pointHeights[i] - minH) / (maxH - minH);
    }

    for(int i = 0; i < gridHeights.length; i++) {
      for(int j = 0; j < gridHeights[i].length; j++) {
        gridHeights[i][j] = (gridHeights[i][j] - minH) / (maxH - minH);
      }
    }

    triangulation.setHeights(pointHeights);
    triangulation.setGridHeights(gridHeights);

    return triangulation;
  }

  public static Point3D[] selectPoints(Point3D[] points, Path2D boundingPolygon) {

    List<Point3D> selectedPoints = new ArrayList<Point3D>();
    for(Point3D p : points) {
      if(boundingPolygon.contains(p.getX(), p.getY()) && !Double.isNaN(p.getZ())) {
        selectedPoints.add(p);
      }
    }

    return selectedPoints.toArray(new Point3D[selectedPoints.size()]);
  }

  public static Point3D[] selectPoints(Point3D[] points) {

    List<Point3D> selectedPoints = new ArrayList<Point3D>();
    for(Point3D p : points) {
      if(!Double.isNaN(p.getZ())) {
        selectedPoints.add(p);
      }
    }

    return selectedPoints.toArray(new Point3D[selectedPoints.size()]);
  }

  public static Point3D[] selectPointsNaN(Point3D[] points, Path2D boundingPolygon) {

    List<Point3D> selectedPoints = new ArrayList<Point3D>();
    for(Point3D p : points) {
      if(boundingPolygon.contains(p.getX(), p.getY())) {
        selectedPoints.add(p);
      }
    }

    return selectedPoints.toArray(new Point3D[selectedPoints.size()]);
  }

  /**
   * searches for the raster file corresponding to the given pattern
   * @param pattern
   * @param rasterFiles
   * @return
   * @throws IOException
   */
  public static Point3D[][] getRaster(String pattern, File[] rasterFiles) throws IOException {

    for(File file : rasterFiles) {
      String fileName = file.getName();
      if(fileName.contains(pattern)) {
        int fileType = getFileType(file);
        if(fileType == 0) {
          return Reader.readAltiNC(file.getAbsolutePath());
        }else {
          return Reader.readElevationNC(file.getAbsolutePath());
        }

      }
    }

    return null;
  }

  /**
   * returns the file type of a nc-file
   * @param rasterFile
   * @return 0 if its altimeter data, 1 else
   * @throws IOException
   */
  public static int getFileType(File rasterFile) throws IOException {
    Netcdf nc = new NetcdfFile(rasterFile, true);
    if(nc.contains("sla")) {
      return 0;
    }else {
      return 1;
    }
  }


  /** convert times into format: yearmonth e.g. 200106 (June 2001)
   * @param times
   * @return
   */
  public static boolean convertTimes(List<String> times){
    for(int i = 0 ; i < times.size(); i++) {
      String s = times.get(i);
      String[] parts = s.split("\\.");
      int month = (int)Math.ceil((Double.parseDouble(parts[1])/100)*1.2);
      if(month < 10) {
        times.set(i, parts[0] + "0" + month);
      }else {
        times.set(i, parts[0] + month);
      }
    }
    return true;
  }

}
