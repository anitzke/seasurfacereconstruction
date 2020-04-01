package processing;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import geometry.HashableLine;
import shapes3D.Point3D;
import triangulation.Triangle;
import triangulation.TriangulationInstance;
import ucar.ma2.Array;
import ucar.ma2.DataType;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.Variable;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Writer {

  /**
   * @param location
   * @param grid
   * @param longitudes
   * @param latitudes
   * @return
   */
  public static boolean writeAltiNC(String location, double[][] grid, double[] longitudes, double[] latitudes) {
    if(grid.length != longitudes.length) {
      System.out.println("wrong dimensions (longitude)");
    }
    if(grid[0].length != latitudes.length) {
      System.out.println("wrong dimensions (latitudes)");
    }
    try {
      NetcdfFileWriter writer = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf3, location, null);
      writer.setFill(true);

      Dimension lonDim = writer.addDimension(null, "lonDim", longitudes.length);
      Dimension latDim = writer.addDimension(null, "latDim", latitudes.length);
      Dimension xysize = writer.addDimension(null, "xysize", grid.length*grid[0].length);

      List<Dimension> dims = new ArrayList<Dimension>();
      dims.add(lonDim);
      dims.add(latDim);
      dims.add(xysize);

      Variable lon = writer.addVariable(null, "lon", DataType.DOUBLE, "lonDim");
      lon.addAttribute(new Attribute("units", "degrees"));

      Variable lat = writer.addVariable(null, "lat", DataType.DOUBLE, "latDim");
      lat.addAttribute(new Attribute("units", "degrees"));

      Variable sla = writer.addVariable(null, "sla", DataType.DOUBLE, "xysize");
      sla.addAttribute(new Attribute("scale_factor", 1.0));
      sla.addAttribute(new Attribute("add_offset", 0.0));
      sla.addAttribute(new Attribute("node_offset", 1.0));
      //			z.addAttribute(new Attribute("valid_range", Array.factory(new double[] {minAnom, maxAnom})));

      writer.addGroupAttribute(null, new Attribute("title", ""));
      writer.addGroupAttribute(null, new Attribute("source", ""));

      writer.create();

      writer.write(lon, Array.factory(longitudes));
      writer.write(lat, Array.factory(latitudes));

      double[] zData = new double[grid.length*grid[0].length];

      int counter = 0;
      for(int j = 0; j < grid[0].length; j++) {
        for(int i = 0; i < grid.length; i++) {
          zData[counter] = grid[i][grid[0].length-1-j];
          counter++;
        }
      }
      writer.write(sla, Array.factory(zData));

      writer.close();
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    } catch (InvalidRangeException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }

  /**
   * @param location
   * @param lines
   * @return
   */
  public static boolean writeLog(String location, List<String> lines) {
    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }

  /**
   * @param location
   * @param points
   * @return
   */
  public static boolean writePointsToCSV(String location, List<Point2D> points) {
    List<String> lines = new ArrayList<String>();
    String firstLine = "\"shapeid\",\"x\",\"y\"";
    lines.add(firstLine);

    int i = 0;
    for(Point2D p : points) {
      String shapeID = "\"" + i + "\"";
      String A = shapeID + ",\"" + String.format("%.3f", p.getX()) + "\",\"" + String.format("%.3f", p.getY()) + "\"";
      lines.add(A);
      i++;
    }

    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }
  
  /**
   * @param location
   * @param points 3D
   * @return
   */
  public static boolean writePointsToCSV(String location, Point3D[] points) {
    List<String> lines = new ArrayList<String>();
    String firstLine = "\"shapeid\",\"x\",\"y\"";
    lines.add(firstLine);

    int i = 0;
    for(Point3D p : points) {
      String shapeID = "\"" + i + "\"";
      String A = shapeID + ",\"" + p.getX() + "\",\"" + p.getY() + "\",\"" + p.getZ() + "\"";
      lines.add(A);
      i++;
    }

    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }
  
  /**
   * @param location
   * @param points 3D
   * @return
   */
  public static boolean writePointsToCSV(String location, Point2D[] points, double[] heights) {
    List<String> lines = new ArrayList<String>();
    String firstLine = "shapeid,x,y,z";
    lines.add(firstLine);

    int i = 0;
    for(Point2D p : points) {
      String shapeID = Integer.toString(i);
      String A = shapeID + "," + p.getX() + "," + p.getY() + "," + heights[i];
      lines.add(A);
      i++;
    }

    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }

  /**
   * @param location
   * @param triangulation
   * @return
   */
  public static boolean writeTrianglesToCSV(String location, TriangulationInstance triangulation) {
    List<String> lines = new ArrayList<String>();
    String firstLine = "\"shapeid\",\"x\",\"y\",\"z\"";
    lines.add(firstLine);

    int i = 0;
    for(Triangle t:triangulation.getTriangles()) {
      String shapeID = "\"" + i + "\"";
      String A = shapeID + ",\"" + t.getA().getX() + "\",\"" + t.getA().getY() + "\",\"" + t.gethA() + "\"";
      String B = shapeID + ",\"" + t.getB().getX() + "\",\"" + t.getB().getY() + "\",\"" + t.gethB() + "\"";
      String C = shapeID + ",\"" + t.getC().getX() + "\",\"" + t.getC().getY() + "\",\"" + t.gethC() + "\"";
      lines.add(A);
      lines.add(B);
      lines.add(C);
      i++;
    }

    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;

  }
  
  /**
   * @param location
   * @param list
   * @return
   */
  public static boolean writeTrianglesToWKT(String location, List<Triangle> list) {
    List<String> lines = new ArrayList<String>();
    String firstLine = "id, triangle";
    lines.add(firstLine);

    int i = 0;
    for(Triangle t: list) {
      String prefix = i + ", \"POLYGON((";
      String A = String.format("%.3f %.3f ", t.getA().getX(), t.getA().getY()) + t.gethA() + ",";
      String B = String.format("%.3f %.3f ", t.getB().getX(), t.getB().getY()) + t.gethB() + ",";
      String C = String.format("%.3f %.3f ", t.getC().getX(), t.getC().getY()) + t.gethC() + ",";
      String A2 = String.format("%.3f %.3f ", t.getA().getX(), t.getA().getY()) + t.gethA();
//      String B = t.getB().getX() + " " + t.getB().getY() + " " + t.gethB() + ",";
//      String C = t.getC().getX() + " " + t.getC().getY() + " " + t.gethC() + ",";
//      String A2 = t.getA().getX() + " " + t.getA().getY() + " " + t.gethA();
      String suffix = "))\"";
      String line = prefix + A + B + C + A2 + suffix;
      lines.add(line);
      i++;
    }

    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }
  
  /**
   * @param location
   * @param triangulation
   * @return
   */
  public static boolean writeTrianglesToWKT(String location, TriangulationInstance triangulation) {
	  return writeTrianglesToWKT(location, triangulation.getTriangles());	  
  }

  /**
   * @param location
   * @param points
   * @return
   */
  public static boolean writePointsToWKT(String location, List<Point3D> points) {
    List<String> lines = new ArrayList<String>();
    String firstLine = "id, points";
    lines.add(firstLine);

    int i = 0;
    for(Point3D p:points) {
      String prefix = i + ", \"POINT(";
      String A = p.getX() + " " + p.getY() + " " + p.getZ();
      String suffix = ")\"";
      String line = prefix + A + suffix;
      lines.add(line);
      i++;
    }

    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }		
    return true;
  }

  /**
   * writes the given Triangulation to the given location as txt-file <br>
   * file-format: <br>
   * l1: number of points n <br>
   * l2: number of triangles m <br>
   * l3 - ln+2: points as a sequence of coordinates seperated by a single , <br>
   * ln+2 - ln+m+2: triangles as a sequenze of integers corresponding to indices in the list of points
   * @param location
   * @param triangulation
   * @return true if writing was successful, false otherwise
   */
  public static boolean writeTriangulation(String location, TriangulationInstance triangulation) {

    Point2D[] points = triangulation.getPoints();
    double[] heights = triangulation.getHeights();


    HashMap<Point2D, Integer> pointMap = new HashMap<Point2D, Integer>();
    for(int i = 0; i < points.length; i++) {
      pointMap.put(points[i], i);
    }

    List<String> lines = new ArrayList<String>();

    String nPoints = ""+points.length;
    lines.add(nPoints);
    String nTriangles = ""+triangulation.getTriangles().size();
    lines.add(nTriangles);

    for(int i = 0; i < points.length; i++) {
//      String pString = points[i].getX() + "," + points[i].getY() + "," + heights[i];
      String pString = String.format("%.3f, %.3f,", points[i].getX(), points[i].getY()) + heights[i];
      lines.add(pString);
    }

    for(Triangle t:triangulation.getTriangles()) {
      String tString = pointMap.get(t.getA()) + "," + pointMap.get(t.getB()) + "," + pointMap.get(t.getC());
      lines.add(tString);
    }

    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }

  /**
   * @param location
   * @param grid
   * @param anomalies
   * @return
   */
  public static boolean writeAnomaliesNC(String location, Point2D[][] grid, double[][] anomalies) {
    try {
      NetcdfFileWriter writer = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf3, location, null);

      double minAnom = Double.MAX_VALUE;
      double maxAnom = Double.MIN_VALUE;
      for(int i = 0; i < anomalies.length; i++) {
        for(int j = 0; j < anomalies[0].length; j++) {
          if(anomalies[i][j] > maxAnom)maxAnom = anomalies[i][j];
          if(anomalies[i][j] < minAnom)minAnom = anomalies[i][j];
        }
      }			
      Dimension side = writer.addDimension(null, "side", 2);
      Dimension xyside = writer.addDimension(null, "xysize", grid.length*grid[0].length);

      List<Dimension> dims = new ArrayList<Dimension>();
      dims.add(side);
      dims.add(xyside);

      Variable x_range = writer.addVariable(null, "x_range", DataType.DOUBLE, "side");
      x_range.addAttribute(new Attribute("units", "meters"));

      Variable y_range = writer.addVariable(null, "y_range", DataType.DOUBLE, "side");
      y_range.addAttribute(new Attribute("units", "meters"));

      Variable z_range = writer.addVariable(null, "z_range", DataType.DOUBLE, "side");
      z_range.addAttribute(new Attribute("units", "meters"));

      Variable spacing = writer.addVariable(null, "spacing", DataType.DOUBLE, "side");

      Variable dimension = writer.addVariable(null, "dimension", DataType.INT, "side");

      Variable z = writer.addVariable(null, "z", DataType.DOUBLE, "xysize");
      z.addAttribute(new Attribute("scale_factor", 1.0));
      z.addAttribute(new Attribute("add_offset", 0.0));
      z.addAttribute(new Attribute("node_offset", 1.0));
      z.addAttribute(new Attribute("valid_range", Array.factory(new double[] {minAnom, maxAnom})));

      writer.addGroupAttribute(null, new Attribute("title", ""));
      writer.addGroupAttribute(null, new Attribute("source", ""));

      writer.create();
      Array space = Array.factory(new double[] {grid[1][0].getX()-grid[0][0].getX(), grid[0][1].getY()-grid[0][0].getY()});
      double sh_x = space.getDouble(0)/2;
      double sh_y = space.getDouble(1)/2;

      writer.write(x_range, Array.factory(new double[] {grid[0][0].getX()-sh_x, grid[grid.length-1][0].getX()+sh_x}));
      writer.write(y_range, Array.factory(new double[] {grid[0][0].getY()-sh_y, grid[0][grid[0].length-1].getY()+sh_y}));


      writer.write(z_range, Array.factory(new double[] {minAnom, maxAnom}));

      writer.write(spacing, space);
      writer.write(dimension, Array.factory(new int[] {grid.length, grid[0].length}));

      double[] zData = new double[anomalies.length*anomalies[0].length];

      int counter = 0;
      for(int j = 0; j < anomalies[0].length; j++) {
        for(int i = 0; i < anomalies.length; i++) {
          zData[counter] = anomalies[i][anomalies[0].length-1-j];
          counter++;
        }
      }
      writer.write(z, Array.factory(zData));

      writer.close();
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    } catch (InvalidRangeException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }

  public static boolean writeScores(String fileName, double... scores) {
    List<String> lines = new ArrayList<String>();
    for(double d:scores) {
      lines.add(String.valueOf(d));
    }
    Path file = Paths.get(fileName);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }
    return true;
  }

  public static void write(String name, Point2D[] points, double[] dataX, double[] dataY, double noise) {
    List<String> lines = new ArrayList<String>();
    String xData = "Parameters for x-direction: min " + dataX[0] + ", max " + dataX[1] + ", n points: " + dataX[2];
    lines.add(xData);
    String yData = "Parameters for y-direction: min " + dataY[0] + ", max " + dataY[1] + ", n points: " + dataY[2];
    lines.add(yData);
    String n = "Noise: " + noise;
    lines.add(n);

    int c = 0;
    for(Point2D p:points) {
      String point = "Point" + c + ": x= " + p.getX() + " y= " + p.getY();
      lines.add(point);
      c++;
    }

    Path file = Paths.get("ErrorData/" + name + ".txt");
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * 
   * @param data[0] = name
   * @param data[1] = dataX
   * @param data[2] = dataY
   * @param data[3] = noise
   * @param data[4] = points
   * @param data[5] = heights
   * @param data[6] = dataGridX
   * @param data[7] = dataGridY
   * @param data[8] = gridNoise
   * @param data[9] = gridPoints
   * @param data[10] = gridHeights				
   */
  public static void write(Object[] data) {
    List<String> lines = new ArrayList<String>();
    String name = (String)data[0];
    double[] dataX = (double[])data[1];
    double[] dataY = (double[])data[2];
    double noise = (double)data[3];
    Point2D[] points = (Point2D[])data[4];
    double[] heights = (double[])data[5];
    double[] dataGridX = (double[])data[6];
    double[] dataGridY = (double[])data[7];
    double gridNoise = (double)data[8];
    Point2D[] gridPoints = (Point2D[])data[9];
    double[] gridHeights = (double[])data[10];


    String xData = "#Parameters for x-direction: min " + dataX[0] + ", max " + dataX[1] + ", n points: " + dataX[2];
    lines.add(xData);
    String yData = "#Parameters for y-direction: min " + dataY[0] + ", max " + dataY[1] + ", n points: " + dataY[2];
    lines.add(yData);
    String n = "#Noise: " + noise;
    lines.add(n);

    for(int c = 0; c < points.length; c++) {
      String point = "Point" + c + ": x= " + points[c].getX() + " y= " + points[c].getY() + " height= " + heights[c];
      lines.add(point);
    }

    String xDataGrid = "#Parameters for x-direction: min " + dataGridX[0] + ", max " + dataGridX[1] + ", n points: " + dataGridX[2];
    lines.add(xDataGrid);
    String yDataGrid = "#Parameters for y-direction: min " + dataGridY[0] + ", max " + dataGridY[1] + ", n points: " + dataGridY[2];
    lines.add(yDataGrid);
    String n2 = "#GridNoise: " + gridNoise;
    lines.add(n2);

    for(int c = 0; c < gridPoints.length; c++) {
      String point = "GridPoint" + c + ": x= " + gridPoints[c].getX() + " y= " + gridPoints[c].getY() + " height= " + gridHeights[c];
      lines.add(point);
    }

    Path file = Paths.get("ErrorData/" + name + ".txt");
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static boolean writeDBMFile(String name, Point2D[] points) {
    List<String> lines = new ArrayList<String>();
    lines.add(String.valueOf(points.length));
    for(Point2D p:points) {
      lines.add(p.getX() + " " + p.getY());
    }
    Path file = Paths.get("ErrorData/" + name + ".txt");
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
      return true;
    } catch (IOException e) {
      e.printStackTrace();
      return false;
    }		
  }

  public static void writeEdgeMap(String location, HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edge,
      List<HashableLine> convexHull) {
    List<String> lines = new ArrayList<String>();
    String firstLine = "edge,left,right,type";
    lines.add(firstLine);

    for(Map.Entry<HashableLine, ArrayList<LinkedList<Integer>>> entry : edge.entrySet()) {
      HashableLine key = entry.getKey();
      String out =  "[" + key.getX1() + " " + key.getY1() + " " + key.getX2() + " " + key.getY2() + 
          "],[" ;
      
      ArrayList<LinkedList<Integer>> val = entry.getValue();
      String left = "";
      String right = "";
      if(!val.get(0).isEmpty()) {
        for(int l : val.get(0)) {
          left = left + l + " ";
        }
        left = left.substring(0, left.length() - 1);
      }
      if(!val.get(1).isEmpty()) {
        for(int r : val.get(1)) {
          right = right + r + " ";
        }
        right = right.substring(0, right.length() - 1);
      }      
      out = out + left + "],[" + right + "],";
      
      if(val.get(0).isEmpty() || val.get(1).isEmpty()) {
        if(convexHull.contains(key)) {
          out = out + 0;
        } else {
          out = out + 1;
        }
      } else {
        out = out + 2;
      }
      
      lines.add(out);
    }


    Path file = Paths.get(location);
    try {
      Files.write(file, lines, Charset.forName("UTF-8"));
    } catch (IOException e) { 
      e.printStackTrace();
    }
  }
}
