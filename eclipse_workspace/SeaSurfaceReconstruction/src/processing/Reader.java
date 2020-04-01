package processing;

import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.imageio.ImageIO;

import shapes3D.Point3D;
import triangulation.Triangle;
import triangulation.TriangulationInstance;
import ucar.netcdf.Netcdf;
import ucar.netcdf.NetcdfFile;
import ucar.netcdf.Variable;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Reader {

  public static List<Point3D[]> readPointWKT(String fileName) {
    Path path = Paths.get(fileName);
    List<String> lines = new ArrayList<String>();
    List<Point3D[]> points = new ArrayList<Point3D[]>();

    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    if(!lines.isEmpty()) {
      Point3D[] pointArray = new Point3D[lines.size()-1];
      for(int i = 1; i < lines.size(); i++) {
        String[] line = lines.get(i).split(", \"POINT\\(|\\)\"| ");
        pointArray[i-1] = new Point3D.Double(Double.parseDouble(line[1]), Double.parseDouble(line[2]), Double.parseDouble(line[3]));
      }
      points.add(pointArray);
    }

    return points;
  }

  public static Path2D readPolygon(String fileName) {
    Path path = Paths.get(fileName);
    List<String> lines = null;

    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    if(lines != null) {
      Path2D polygon = new Path2D.Double();
      boolean first = true;
      for(int i = 1; i < lines.size(); i++) {
        if(!lines.get(i).isEmpty()) {
          String[] parts = lines.get(i).split(",");
          if(first) {
            double x = Double.parseDouble(parts[1].substring(1, parts[1].length()-2));
            double y = Double.parseDouble(parts[2].substring(1, parts[2].length()-2));
            polygon.moveTo(x, y);
            first = false;
          }else {
            double x = Double.parseDouble(parts[1].substring(1, parts[1].length()-2));
            double y = Double.parseDouble(parts[2].substring(1, parts[2].length()-2));
            polygon.lineTo(x, y);
          }
        }
      }
      polygon.closePath();
      return polygon;
    }
    return null;
  }
  
  public static Path2D readPolygon(String fileName, String bla) {
    Path path = Paths.get(fileName);
    List<String> lines = null;

    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      e.printStackTrace();
    }
    if(lines != null) {
      Path2D polygon = new Path2D.Double();
      boolean first = true;
      for(int i = 1; i < lines.size(); i++) {
        if(!lines.get(i).isEmpty()) {
          String[] parts = lines.get(i).split(",");
          if(first) {
            double x = Double.parseDouble(parts[0]);
            double y = Double.parseDouble(parts[1]);
            polygon.moveTo(x, y);
            first = false;
          }else {
            double x = Double.parseDouble(parts[0]);
            double y = Double.parseDouble(parts[1]);
            polygon.lineTo(x, y);
          }
        }
      }
      polygon.closePath();
      return polygon;
    }
    return null;
  }

  public static int[] readPoints2Use(String fileName) {
    Path path = Paths.get(fileName);
    List<String> lines = null;
    int[] pointNumbers = null;
    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    if(lines != null) {
      String[] parts = lines.get(0).split(",");
      pointNumbers = new int[parts.length];
      for(int i = 0; i < parts.length; i++) {
        pointNumbers[i] = Integer.parseInt(parts[i]);
      }
    }
    return pointNumbers;
  }

  public static TriangulationInstance readTriangulation(String fileName){
    Path path = Paths.get(fileName);
    List<String> lines = null;
    TriangulationInstance triangulation = new TriangulationInstance();
    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    if(lines != null) {
      int nPoints = Integer.parseInt(lines.get(0));
      int nTriangles = Integer.parseInt(lines.get(1));

      Point2D[] points = new Point2D[nPoints];
      double[] heights = new double[nPoints];
      int c = 0;
      for(int i = 2; i < nPoints+2; i++) {
        String[] parts = lines.get(i).split(",");
        points[c] = new Point2D.Double(Double.parseDouble(parts[0]), Double.parseDouble(parts[1]));
        heights[c] = Double.parseDouble(parts[2]);
        c++;
      }
      triangulation.setPoints(points);
      triangulation.setHeights(heights);

      List<Triangle> triangles = new ArrayList<Triangle>();
      for(int i = nPoints+2; i < nPoints+2+nTriangles; i++) {
        String[] parts = lines.get(i).split(",");
        int A = Integer.parseInt(parts[0]);
        int B = Integer.parseInt(parts[1]);
        int C = Integer.parseInt(parts[2]);
        Triangle t = new Triangle(points[A], heights[A], points[B], heights[B], points[C], heights[C]);
        triangles.add(t);
      }
      triangulation.setTriangles(triangles);
      String name = path.getFileName().toString();
      int pos = name.lastIndexOf(".");
      if (pos > 0) {
        name = name.substring(0, pos);
      }
      triangulation.setName(name);
      System.out.println("File " + path.getFileName() + " successfully read from " + path.getParent().toAbsolutePath());
      return triangulation;
    }
    return null;
  }

  public static Point3D[][] readTiff(String fileName, double[] gridData){

    try {
      File file = new File(fileName);
      BufferedImage img = ImageIO.read(file);
      int width = img.getWidth();
      int height = img.getHeight();
      //			Raster tile = img.getTile(0, 0);
      WritableRaster raster = img.getRaster();
      //			double[] data = new double[3];
      //			double [] pixel = tile.getPixel(0, 0, data);
      Point3D[][] points = new Point3D[width][height];
      double originX = gridData[0];
      double originY = gridData[1];
      double deltaX = gridData[2];
      double deltaY = gridData[3];
      for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
          //					if(i == 3600 && j == 3600) {
          //						System.out.println("hi");
          //					}
          double[] point = raster.getPixel(i, height-j-1, new double[3]);
          if(point[0] < -1000) point[0] = -500;
          points[i][j] = new Point3D.Double(originX+i*deltaX, originY+j*deltaY, point[0]);
        }
      }
      return points;
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }


    return null;
  }

  public static Point3D[][] readElevationNC(String fileName) throws IOException{


    Netcdf nc = new NetcdfFile(fileName, true);
    int[] index = new int[1];

    Variable lon = nc.get("lon");
    int nlons = lon.getLengths()[0];
    double[] lons = new double[nlons];
    for(int ilon = 0; ilon < nlons; ilon++) {
      index[0] = ilon;
      double l = lon.getDouble(index);
      if(l > 180) {
        lons[ilon] = l - 360;
      }else {
        lons[ilon] = l;
      }
    }

    Variable lat = nc.get("lat");
    int nlats = lat.getLengths()[0];
    double[] lats = new double[nlats];
    for(int ilat = 0; ilat < nlats; ilat++) {
      index[0] = ilat;
      lats[ilat] = lat.getDouble(index);
    }

    Variable band = nc.get("Band1");
    float fill = band.getAttribute("_FillValue").getNumericValue().floatValue();
    int[] ix = new int[2];
    double[][] z = new double[nlons][nlats];
    for(int ilat = 0; ilat < nlats; ilat++) {
      ix[0] = ilat;
      for(int ilon = 0; ilon < nlons; ilon++) {
        ix[1] = ilon;
        z[ilon][ilat] = band.getFloat(ix);
      }
    }

    Point3D[][] ps = new Point3D[nlons][nlats];

    for(int i = 0; i < nlons; i++) {
      for(int j = 0; j < nlats; j++) {
        if(z[i][j] != fill) {
          ps[i][j] = new Point3D.Double(lons[i], lats[j], z[i][j]);
        }else {
          ps[i][j] = new Point3D.Double(lons[i], lats[j], Double.NaN);
        }
      }
    }			
    return ps;

  }

  public static Point3D[][] readAltiNC(String fileName) throws IOException {

    Netcdf nc = new NetcdfFile(fileName, true);
    int[] index = new int[1];

    Variable lon = nc.get("lon");
    int nlons = lon.getLengths()[0];
    double[] lons = new double[nlons];
    for(int ilon = 0; ilon < nlons; ilon++) {
      index[0] = ilon;
      double l = lon.getDouble(index);
      if(l > 180) {
        lons[ilon] = l - 360;
      }else {
        lons[ilon] = l;
      }
    }

    Variable lat = nc.get("lat");
    int nlats = lat.getLengths()[0];
    double[] lats = new double[nlats];
    for(int ilat = 0; ilat < nlats; ilat++) {
      index[0] = ilat;
      lats[ilat] = lat.getDouble(index);
    }

    Variable sla = nc.get("sla"); // sea level average
    float fill = sla.getAttribute("_FillValue").getNumericValue().floatValue();
    int[] ix = new int[2];
    double[][] slas = new double[nlons][nlats];
    for(int ilat = 0; ilat < nlats; ilat++) {
      ix[0] = ilat;
      for(int ilon = 0; ilon < nlons; ilon++) {
        ix[1] = ilon;
        slas[ilon][ilat] = sla.getFloat(ix); // create raster, set sea level average to its lon lat
      }
    }

    List<double[]> slaCopy = Arrays.asList(slas);
    List<double[]> sortedSlas = new ArrayList<double[]>(slaCopy);
    Collections.sort(sortedSlas, Comparator.comparing(s -> lons[slaCopy.indexOf(s)]));

    Arrays.sort(lons); // now from -180 to 180° instead of 0 to 360°

    Point3D[][] ps = new Point3D[nlons][nlats];

    for(int i = 0; i < nlons; i++) {
      for(int j = 0; j < nlats; j++) {
        if(sortedSlas.get(i)[j] != fill) {          
          ps[i][j] = new Point3D.Double(lons[i], lats[j], Calculations.rounding(sortedSlas.get(i)[j]*100, 1000));
        }else {
          ps[i][j] = new Point3D.Double(lons[i], lats[j], Double.NaN);
        }
      }
    }	

    return ps;		

  }

  // changed return type from Point3D[][] to List<Point3D[]>
  public static List<Point3D[]> readD35(String fileName, List<Integer> pointNumbers, List<String> months) {
    Path path = Paths.get(fileName);
    List<String> lines = null;
    List<Point3D[]> ps = new ArrayList<Point3D[]>();
    List<Double> longitude = new ArrayList<Double>();
    List<Double> latitude = new ArrayList<Double>();

    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch(IOException e){
      e.printStackTrace();
    }

    if(lines != null) {
      // first line in file 
      String[] parts = lines.get(2).split(" ");
      if(pointNumbers != null) {
        for(int j = 1; j < parts.length; j++) {
          if(!parts[j].isEmpty()) {
            pointNumbers.add(Integer.parseInt(parts[j]));
          }
        }
      }

      // second line in file contains the longitude of the tide gauge stations
      parts = lines.get(1).split(" ");
      for(int j = 1; j < parts.length; j++) {
        if(!parts[j].isEmpty()) {
          longitude.add(Double.parseDouble(parts[j]));
        }
      }

      // third line in file contains the latitude of the tide gauge stations
      parts = lines.get(2).split(" ");
      for(int j = 1; j < parts.length; j++) {
        if(!parts[j].isEmpty()) {
          latitude.add(Double.parseDouble(parts[j]));
        }
      }

      // from the sixth line the measurements per month starts
      for(int i = 8; i < lines.size(); i++) {
        Point3D[] psTemp = new Point3D[longitude.size()];
        parts = lines.get(i).split(" ");
        int k = 0;
        if(months != null) {
          months.add(parts[0]);
        }
        for(int j = 1; j < parts.length; j++) {
          if(!parts[j].isEmpty()) {
            double z = Double.parseDouble(parts[j]);
            // z divided by 1000 since height (z) is given in mm  (mm -> m)
            psTemp[k] = new Point3D.Double(longitude.get(k), latitude.get(k), Calculations.rounding(z/10, 1000));
            k++;
          }
        }
        ps.add(psTemp);
      }
      //			return ps.toArray(new Point3D[ps.size()][]);
      return ps;

    }
    return new ArrayList<Point3D[]>();
  }

  public static double[] readRandomPoints(String fileName) {
    Path path = Paths.get(fileName);
    List<String> lines = null;
    double[] coords = new double[0];

    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    if(lines != null) {
      coords = new double[(lines.size()-1)*2];
      for(int i = 0; i < lines.size()-1; i++) {
        String[] parts = lines.get(i+1).split(",");
        coords[i*2] = Double.parseDouble(parts[0]);
        coords[i*2+1] = Double.parseDouble(parts[1]);
      }
      return coords;
    }
    return null;
  }

  public static Point2D[] read(String fileName) {
    Path path = Paths.get(fileName);
    List<String> lines = null;
    List<Point2D> ps = new ArrayList<Point2D>();
    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      e.printStackTrace();
    }
    if(lines != null) {
      for(String s : lines) {
        if(s.startsWith("Point")) {
          String[] parts = s.split(" ");
          Point2D p = new Point2D.Double(Double.parseDouble(parts[2]), Double.parseDouble(parts[4]));
          ps.add(p);
        }
      }
      return ps.toArray(new Point2D[ps.size()]);
    }
    return null;
  }

  public static Point2D[] readBDMFile(String fileName) {
    Path path = Paths.get(fileName);
    List<String> lines = null;
    List<Point2D> ps = new ArrayList<Point2D>();
    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      e.printStackTrace();
    }
    int counter = 0;
    if(lines != null) {
      for(String line : lines) {
        if((line.charAt(0) >= '0') && (line.charAt(0) <= '9')) {
          double num = Double.parseDouble(line);
          for(int i = counter+1; i < counter+1+num; i++) {
            String[] parts = lines.get(i).split(" ");
            Double x = null;
            Double y = null;
            for(String s : parts) {
              if(!s.isEmpty()) {
                if(x == null) {
                  x = Double.parseDouble(s);
                }else {
                  y = Double.parseDouble(s);
                }
              }
            }
            Point2D p = new Point2D.Double(x, y);
            ps.add(p);
          }
          break;
        }
        counter++;
      }
      return ps.toArray(new Point2D[ps.size()]);
    }
    return null;
  }

  public static List<Triangle> readBDMTriangulation(String fileName) {
    Path path = Paths.get(fileName);
    List<String> lines = null;
    List<Point2D> ps = new ArrayList<Point2D>();
    List<Triangle> ts = new ArrayList<Triangle>();
    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    int counter = 0;
    int lastPoint = 0;
    if(lines != null) {
      for(String line : lines) {
        if((line.charAt(0) >= '0') && (line.charAt(0) <= '9')) {
          int num = Integer.parseInt(line);
          lastPoint = counter+num;
          for(int i = counter+1; i < counter+1+num; i++) {
            String[] parts = lines.get(i).split(" ");
            Double x = null;
            Double y = null;
            for(String s : parts) {
              if(!s.isEmpty()) {
                if(x == null) {
                  x = Double.parseDouble(s);
                }else {
                  y = Double.parseDouble(s);
                }
              }
            }
            Point2D p = new Point2D.Double(x, y);
            ps.add(p);
          }
          break;
        }
        counter++;
      }
    }else {
      return null;
    }
    for(int i = lastPoint+1;i<lines.size();i++) {
      if(lines.get(i).startsWith("triangles")) {
        int num = Integer.parseInt(lines.get(i+1));
        for(int j = i+2; j < i+2+num; j++) {
          String[] parts = lines.get(j).split(" ");
          int[] points = {Integer.parseInt(parts[0]),
              Integer.parseInt(parts[1]),
              Integer.parseInt(parts[2])};
          Triangle t = new Triangle(ps.get(points[0]), 
              ps.get(points[1]), 
              ps.get(points[2]));
          ts.add(t);
        }
      }
    }


    return ts;
  }

}
