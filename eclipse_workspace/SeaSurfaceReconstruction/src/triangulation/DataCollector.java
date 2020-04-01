package triangulation;

import java.awt.geom.Path2D;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import com.vividsolutions.jts.geom.Envelope;

import processing.DataHandler;
import processing.PointSelector;
import processing.Projection;
import processing.Reader;
import shapes3D.Point3D;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class DataCollector {

  private RawData mData;
  private static HashMap<String, String> parameterMap = new HashMap<String, String>();
  private static Projection projection = new Projection();

  public DataCollector() {

  }

  public DataCollector(String[] params) {
    try {
      mData = collect(params);
    } catch (IOException | UnrecognizedFileFormatException e) {
      e.printStackTrace();
    }
  }

  public static String getParameter(String param) {
    return parameterMap.get(param);
  }
  
  public static Projection getProjection() {
    return projection;
  }


  /**
   * collects the data specified by command line input
   * @param params
   * @return returns an instance of RawData, depending on mode
   * @throws IOException
   * @throws UnrecognizedFileFormatException
   */
  public RawData collect(String[] params) throws IOException, UnrecognizedFileFormatException {

    List<String> obligatories = sortParameters(params);
    RawData mData = null;
    parameterMap.get("-ev");

    switch(parameterMap.get("-mode")) {
    case "compute":
      RawData.Triangulate tempData = new RawData.Triangulate();
      mData = createComputeData(tempData, obligatories);
      break;
    case "evaluate": case "reconstruct":
      RawData.Evaluate tempData3 = new RawData.Evaluate();
      mData = createEvaluateData(tempData3, obligatories);
      break;
    default:
      break;
    }
//    
//    if(!parameterMap.get("-projection").isEmpty()) {
//    	mData = Projection.project(mData);
//    }

    this.mData = mData;
    return mData;
  }

  public RawData getData() {
    return mData;
  }

  /**
   * Sets the default map for parameters and their values. 
   * Also updates the map to the specified parameters. 
   * @param params
   * @return list of essential inputs
   */
  public List<String> sortParameters(String[] paramsString) {
    
    parameterMap.put("-mode", "evaluate");
    parameterMap.put("-candidates", "kODT");
    parameterMap.put("-tc", "0");
    parameterMap.put("-threads", "1");
    parameterMap.put("-vh", "nonan");
    parameterMap.put("-ev", "all");
    parameterMap.put("-on", "all");
    parameterMap.put("-rt", "12");
    parameterMap.put("-store", "experimentNew");
    parameterMap.put("-ma", "0month");
    parameterMap.put("-optimization", "min");
    parameterMap.put("-projection", "");
    parameterMap.put("-stationTraining", "");
  
    Path path = Paths.get(paramsString[0]);
    List<String> lines = new ArrayList<String>();
    List<String> obligatories = new ArrayList<String>();

    try {
      lines = Files.readAllLines(path, StandardCharsets.UTF_8);
    } catch (IOException e) {
      e.printStackTrace();
    }

    if(!lines.isEmpty()) {

    for(int i = 0; i < lines.size(); i++) {
      String[] params = lines.get(i).split(": ");
      if(params[0].equals("-mode")) {
        parameterMap.put("-mode", params[1]);
      } else if(params[0].equals("-candidates")) {
        parameterMap.put("-candidates", params[1]);
      } else if(params[0].equals("-tc")) {
        parameterMap.put("-tc", params[1]);
      }else if(params[0].equals("-threads")){
        parameterMap.put("-threads", params[1]);
      } else if(params[0].equals("-vh")){
        parameterMap.put("-vh", params[1]);
      } else if (params[0].equals("-ev")) {
        parameterMap.put("-ev", params[1]);				
      } else if (params[0].equals("-on")) {
        parameterMap.put("-on", params[1]);       
      } else if (params[0].equals("-rt")) {
        parameterMap.put("-rt", params[1]);       
      }else if (params[0].equals("-store")) {
        parameterMap.put("-store", params[1]); 
      } else if (params[0].equals("-ma")) {
        parameterMap.put("-ma", params[1]); 
      } else if (params[0].equals("-optimization")) {
        parameterMap.put("-optimization", params[1]); 
      } else if (params[0].equals("-projection")) {
        parameterMap.put("-projection", params[1]); 
      } else if (params[0].equals("-stationTraining")) {
        parameterMap.put("-stationTraining", params[1]); 
      } else {
        obligatories.add(params[1]);
      }
    }
    }

    return obligatories;
  }

  /**
   * creates an instance of RawData.Triangulate (for compute mode)
   * @param mData
   * @param obligatories
   * @return RawData.Triangulate
   * @throws IOException
   * @throws UnrecognizedFileFormatException
   */
  public RawData createComputeData(RawData.Triangulate mData, List<String> obligatories) 
      throws IOException, UnrecognizedFileFormatException {

    //create TriangleCreator
    TriangleCreator tc = createTriangleCreator();
    mData.tc = tc;
    PointSelector ps = createPointSelector();

    String pointDirectory = obligatories.get(0);
    String rasterDirectory = obligatories.get(1);

    // load points and any pattern in their names
    List<String> patterns = new ArrayList<String>();
    List<Point3D[]> points = getPoints(pointDirectory, patterns);
    
    // define training epochI and reconstruction epochJ
    if(parameterMap.get("-ev").equals("all")) {
      mData.epochI = patterns; // year and month in list (199401 for January 1994)
    } else {
      String str[] = parameterMap.get("-ev").split(",");
      mData.epochI = Arrays.asList(str);
    }
    
    if(parameterMap.get("-on").equals("all")) {
      mData.epochJ = patterns; // year and month in list (199401 for January 1994)
    } else {
      String str[] = parameterMap.get("-on").split(",");
      mData.epochJ = Arrays.asList(str);
    }
    
    // in case of reconstruction do
    if(!parameterMap.get("-stationTraining").isEmpty()) {
      List<String> patternsTraining = new ArrayList<String>();
      List<Point3D[]> pointsTraining = getPoints(parameterMap.get("-stationTraining"), patternsTraining);
      for(int i=0; i<mData.epochI.size(); i++) {
        int it = patternsTraining.indexOf(mData.epochI.get(i));
        int ir = patterns.indexOf(mData.epochI.get(i));
        points.set(ir, pointsTraining.get(it));
      }
    }

    // load bounds 
    if(obligatories.size() > 2) {
      mData.pointBounds = Reader.readPolygon(obligatories.get(2),"bla");
      if(ps instanceof PointSelector.BoundsNaN) {
        ps = new PointSelector.BoundsNaN(points.get(0), mData.pointBounds);
      }else {
        ps = new PointSelector.Bounds(points.get(0), mData.pointBounds);
      }
      for(int i = 0; i < points.size(); i++) {
        ps.setPoints(points.get(i));
        Point3D[] selectedPoints = ps.select();
        points.set(i, selectedPoints);
      }
      if(obligatories.size() > 3) {
        mData.rasterBounds = Reader.readPolygon(obligatories.get(3));
      }
      if(obligatories.size() > 4) {
        throw new IOException("whatchu doing?");
      }
    }else {
      // select points of method specified with -vh argument parameter
      // (e.g. no points that have no height/NaN)
      for(int i = 0; i < points.size(); i++) {
        ps.setPoints(points.get(i));
        points.set(i, ps.select());
      }
    }
    mData.points = points;
    mData.identifiers = patterns;
    

    // load raster*s and attempts to match rasters to their points based on patterns 
    // (only if input is a directory) 
    File rasterFile = new File(rasterDirectory);
    List<Point3D[][]> rasterList = new ArrayList<Point3D[][]>();
    

    if(rasterFile.isDirectory()) {
      // all files in directory
      File[] rasterFiles = rasterFile.listFiles();
      if(parameterMap.get("-mode").equals("reconstruct")) {
        for(int i = 0; i < mData.epochI.size(); i++) {
          int ii = mData.identifiers.indexOf(mData.epochI.get(i));
          Point3D[][] entireRaster = DataHandler.getRaster(mData.epochI.get(i), rasterFiles);
          if(entireRaster != null) {
            if(obligatories.size() > 3) {
              rasterList.add(DataHandler.clipRaster(entireRaster, mData.rasterBounds));
            }else {
              rasterList.add(DataHandler.clipRaster(entireRaster, points.get(ii)));
            }
          }else {
            throw new IOException("Could not find a corresponding raster for: " + patterns.get(i));
          }
        }
      } else {
      // for all months/year do
        for(int i = 0; i < patterns.size(); i++) {
          Point3D[][] entireRaster = DataHandler.getRaster(patterns.get(i), rasterFiles);
          if(entireRaster != null) {
            if(obligatories.size() > 3) {
              rasterList.add(DataHandler.clipRaster(entireRaster, mData.rasterBounds));
            }else {
              rasterList.add(DataHandler.clipRaster(entireRaster, points.get(i)));
            }
          }else {
            throw new IOException("Could not find a corresponding raster for: " + patterns.get(i));
          }
        }
      }
    }else {
      int fileType = DataHandler.getFileType(rasterFile);
      if(fileType == 0) {
        rasterList.add(Reader.readAltiNC(rasterFile.getAbsolutePath()));
      }else {
        rasterList.add(Reader.readElevationNC(rasterFile.getAbsolutePath()));
      }

    }

//    if(patterns.size() != rasterList.size() && rasterList.size() > 1) {
//      throw new IOException("Something went very wrong!");
//    }
    mData.rasters = rasterList;
    
    
    if(parameterMap.get("-projection").equals("LCC")) {
      Projection proj = new Projection(mData, "toLCC");
      proj.setProjectedPoints(proj.project(proj.getOriginPoints()));
      proj.setProjectedRasters(proj.projectRaster(proj.getOriginRasters()));
      mData.rasters = proj.getProjectedRasters();
      mData.points = proj.getProjectedPoints();
      projection = proj;
    }
    
    return mData;
  }

  /**
   * creates an instance of RawData.Evaluate (for evaluate mode)
   * @param mData
   * @param obligatories
   * @return RawData.Evaluate
   * @throws IOException
   * @throws UnrecognizedFileFormatException
   */
  public RawData createEvaluateData(RawData.Evaluate mData, List<String> obligatories) throws IOException, UnrecognizedFileFormatException {
    RawData.Evaluate raw = (RawData.Evaluate) createComputeData(mData, obligatories);

    if(raw.epochI.size()==1) {
      Point3D[] singleEvalPoints = new Point3D[0];
      Point3D[][] singleEvalRaster = new Point3D[0][];
      String singleEvalID = "";

      int i = raw.identifiers.indexOf(raw.epochI.get(0));
      singleEvalPoints = raw.points.get(i);
      //raw.points.remove(i);
      if(parameterMap.get("-mode").equals("reconstruct")) {
        singleEvalRaster = raw.rasters.get(0);
      } else {
        singleEvalRaster = raw.rasters.get(i);
      }
      //raw.rasters.remove(i);
      singleEvalID = raw.epochI.get(0);
      //raw.identifiers.remove(i);

      if(singleEvalPoints.length == 0 || singleEvalRaster.length == 0 || singleEvalID.isEmpty()) {
        System.out.println("could not find data for: " + raw.epochI.get(0));
      } else {
        TriangulationInstance t = new TriangulationInstance();
        t.setName(singleEvalID);
        DataHandler.setTriangulation(t, singleEvalPoints, singleEvalRaster);
          raw.setSingleEvaluation(t);
      }

    }

    return raw;
  }

  /**
   * creates an instance of RawData.Validate (for validate mode)
   * @param mData
   * @param obligatories
   * @return RawData.Validate
   * @throws IOException
   * @throws UnrecognizedFileFormatException
   */
  public RawData createValidateData(RawData.Validate mData, List<String> obligatories) {

    String triangulationsDirectory = obligatories.get(0);
    List<TriangulationInstance> triangulations = new ArrayList<TriangulationInstance>();
    List<String> patterns = new ArrayList<String>();
    File triangulationsFile = new File(triangulationsDirectory);

    if(triangulationsFile.isDirectory()) {
      File[] files = triangulationsFile.listFiles();
      for(File file : files) {
        triangulations.add(Reader.readTriangulation(file.getAbsolutePath()));
        patterns.add(getNamePattern(file.getName()));
      }
    }else {
      triangulations.add(Reader.readTriangulation(triangulationsFile.getAbsolutePath()));
      patterns.add(getNamePattern(triangulationsFile.getName()));
    }
    mData.triangulations = triangulations;
    mData.identifiers = patterns;

    return mData;
  }

  public Path2D getTriangulationsBoundingBox(List<TriangulationInstance> triangulations) {
    Path2D BB = new Path2D.Double();
    Envelope env = new Envelope();	

    for(TriangulationInstance t : triangulations) {
      env.expandToInclude(t.getEnvelope());
    }

    BB.moveTo(env.getMinX(), env.getMinY());
    BB.lineTo(env.getMaxX(), env.getMinY());
    BB.lineTo(env.getMaxX(), env.getMaxY());
    BB.lineTo(env.getMinX(), env.getMaxY());
    BB.lineTo(env.getMinX(), env.getMinY());
    BB.closePath();

    return BB;
  }

  /**
   * Reads points in pointDirectory. 
   * Also fills the pattern-list with filenames of point-files or months of a d35-file.
   * @param pointDirectory
   * @param patterns
   * @return List of point-arrays
   * @throws UnrecognizedFileFormatException
   */
  public List<Point3D[]> getPoints(String pointDirectory, List<String> patterns) throws UnrecognizedFileFormatException {
    List<Point3D[]> points = new ArrayList<Point3D[]>();
    if(pointDirectory.endsWith("d35")) {
      //			Point3D[][] pointMatrix = Reader.readD35(pointDirectory, null, patterns);
      //			for(int i = 0; i < pointMatrix.length; i++) {
      //				points.add(pointMatrix[i]);
      //			}
      points = Reader.readD35(pointDirectory, null, patterns);

      //convert times into format: yearmonth e.g. 200106 (June 2001)
      DataHandler.convertTimes(patterns);
    } else if(pointDirectory.endsWith("csv")) {
      points = Reader.readPointWKT(pointDirectory);

      patterns.add(getNamePattern(pointDirectory));
    } else if(new File(pointDirectory).isDirectory()) {
      File[] pointFiles = new File(pointDirectory).listFiles();

      for(int i = 0; i < pointFiles.length; i++) {
        if(pointFiles[i].getName().endsWith("csv")) {
          points.add(Reader.readD35(pointFiles[i].getAbsolutePath(), null, patterns).get(0));
          patterns.add(getNamePattern(pointFiles[i].getName()));
        }else {
          String[] parts = pointDirectory.split("\\.");
          System.out.println("unrecognized point-file format: " + parts[parts.length-1]);
        }
      }		
    }else {
      throw new UnrecognizedFileFormatException(pointDirectory);
    }
    return points;
  }

  /**
   * Finds a pattern or name in a String.
   * Pattern starts after the first "_" or "-" after the last "/", if there are any.
   * Pattern ends one position before the last ".".
   * @param directory
   * @return
   */
  public String getNamePattern(String directory) {
    int nameStart = (directory.lastIndexOf("/") > 0) ? directory.lastIndexOf("/") : 0;
    int patternStart = directory.indexOf("_", nameStart);
    if(patternStart < 0) {
      patternStart = directory.indexOf("-", nameStart);
    }
    int patternEnd = directory.lastIndexOf(".");
    return directory.substring(patternStart+1, patternEnd);
  }

  public static class UnrecognizedFileFormatException extends Exception{

    private static final long serialVersionUID = 1L;
    private String directory;

    public UnrecognizedFileFormatException(String pointDirectory) {
      directory = pointDirectory;
    }

    @Override
    public String toString() {
      String[] parts = directory.split("\\.");
      String message = "unrecognized point-file format: " + parts[parts.length-1] + "\n" + 
          "recognized formats: d35, csv";
      return message;
    }

  }

  /**
   * creates an instance of PointSelector based on the "-vh"-parameter in the parameterMap.
   * @return PointSelector
   */
  public PointSelector createPointSelector() {
    String psString = parameterMap.get("-vh");

    if(psString.equals("nonan")) {
      return new PointSelector.NoNaN();
    }else if(psString.equals("bounds")){
      return new PointSelector.Bounds();
    }else if(psString.equals("all")) {
      return new PointSelector.All();
    }else if(psString.equals("allnan")) {
      return new PointSelector.AllNaN();
    }else if(psString.equals("boundsnan")) {
      return new PointSelector.BoundsNaN();
    }else {
      List<Point3D> vPoints = Arrays.asList(Reader.readPointWKT(psString).get(0));
      PointSelector.Specific ps = new PointSelector.Specific();
      ps.setVariablePoints(vPoints);
      return ps;
    }
  }

  /**
   * creates an instance of TriangleCreator based on the "-tc"-parameter in the parameterMap.
   * @return TriangleCreator
   */
  public TriangleCreator createTriangleCreator() {
    String cMethod = parameterMap.get("-candidates");

    TriangleCreator tc;
    try {
      switch(cMethod) {
      case("all"):
        tc = new TriangleCreator.All();
        break;
      case("kODT"):
        int n = Integer.parseInt(parameterMap.get("-tc"));
        tc = new TriangleCreator.Constrained(n);
        break;
      case("empty"):
        tc = new TriangleCreator.Empty();
      break;
      default:
        tc = new TriangleCreator.Empty();
        break;
      } 
    } catch(Exception e) {
      tc = new TriangleCreator.Empty();
      System.out.println("Valid candidate methods: \"empty\", \"kODT\", \"all\"");
    }

    return tc;
  }

}
