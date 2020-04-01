package triangulation;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import processing.Writer;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class DataDistributor {

  private ResultData results;
  private String topFolder = "experiments/experimentNew";

  public DataDistributor(ResultData results) {
    this.results = results;
  }

  public DataDistributor(ResultData results, String writeDirectory) {
    this.results = results;
    this.topFolder = "experiments/" + writeDirectory;
  }

  public void setResults(ResultData results) {
    this.results = results;
  }

  public void setDirectory(String writeDirectory) {
    this.topFolder = writeDirectory;
  }

  /**
   * writes the results into the folder specified by topFolder.
   */
  public void distribute() {

    String anomalyPath = topFolder + "/simpleAnomalies";
    String squaredAnomalyPath = topFolder + "/squaredAnomalies";
    String interpolatedValuesPath = topFolder + "/interpolations";
    String statisticPath = topFolder + "/statistics";
//    String triangulationPath = topFolder + "/triangulations";
    String triangleShapesPath = topFolder + "/triangleCSVs";
    String validationPath = topFolder + "/validations";
    String timerPath = topFolder + "/times";
//    String candidatesPath = topFolder + "/candidates";
//    String edgePath = topFolder + "/edges";
//    String convexHullPath = topFolder + "/convexHull";
    String sharedStationsPath = topFolder + "/sharedStations";
    String reconstructionPath = topFolder + "/reconstructions";
    String sqaDiffPath = topFolder + "/sqaDiff";
    
    checkCreatePath(topFolder);
    
    if(!results.getAnomalies().isEmpty()) {
      checkCreatePath(anomalyPath);
      for(int i = 0; i < results.getAnomalies().size(); i++) {
        Writer.writeAnomaliesNC(anomalyPath + "/" + results.getNames().get(i) + ".nc", results.getTriangulations().get(i).getGridPoints(), results.getAnomalies().get(i));
      }
    }
    
    if(!results.getSquaredAnomalies().isEmpty()) {
      checkCreatePath(squaredAnomalyPath);
      for(int i = 0; i < results.getSquaredAnomalies().size(); i++) {
        Writer.writeAnomaliesNC(squaredAnomalyPath + "/" + results.getNames().get(i)+ ".nc", results.getTriangulations().get(i).getGridPoints(), results.getSquaredAnomalies().get(i));
      }
    }
    
    if(!results.getSqaDiffs().isEmpty()) {
      checkCreatePath(sqaDiffPath);
      for(int i = 0; i < results.getSqaDiffs().size(); i++) {
        Writer.writeAnomaliesNC(sqaDiffPath + "/" + results.getNames().get(i*2).replace("_ILP", "") + ".nc", results.getTriangulations().get(i*2).getGridPoints(), results.getSqaDiffs().get(i));
      }
    }
    
    if(!results.getInterpolatedValues().isEmpty()) {
      checkCreatePath(interpolatedValuesPath);
      for(int i = 0; i < results.getInterpolatedValues().size(); i++) {
//        List<double[][]> test1 = results.getInterpolatedValues();
//        double[][] test = results.getInterpolatedValues().get(0);
        Writer.writeAnomaliesNC(interpolatedValuesPath + "/" + results.getNames().get(i)+ ".nc", results.getTriangulations().get(i).getGridPoints(), results.getInterpolatedValues().get(i));
//        double[] longs = new double[results.getTriangulations().get(i).getGridPoints().length];
//        double[] lats = new double[results.getTriangulations().get(i).getGridPoints()[0].length];
//        for (int j = 0; j< longs.length; j++) {
//          longs[j] = results.getTriangulations().get(i).getGridPoints()[j][0].getX();
//        }
//        for ( int j = 0; j< lats.length; j++) {
//          lats[j] = results.getTriangulations().get(i).getGridPoints()[0][j].getY();
//        }
//        Writer.writeAltiNC(interpolatedValuesPath + "/" + results.getNames().get(i)+ "_I.nc", results.getInterpolatedValues().get(i), longs, lats);
      }
    }
    
    if(!results.getStatistics().isEmpty()) {
      checkCreatePath(statisticPath);
      Writer.writeLog(statisticPath + "/statisticsLog.txt", results.getStatistics());
    }
    
    if(!results.getInterpolationMeans().isEmpty()) {
      checkCreatePath(reconstructionPath);
      Writer.writeLog(reconstructionPath + "/reconstructionsLog.txt", results.getInterpolationMeans());
    }
    
    if(!results.getValidations().isEmpty()) {
      checkCreatePath(validationPath);
      Writer.writeLog(validationPath + "/validationLog.txt", results.getValidations());
    }
    
    if(!results.getTimes().isEmpty()) {
      checkCreatePath(timerPath);
      Writer.writeLog(timerPath + "/timeLog.txt", results.getTimes());
    }
    
    if(!results.getTriangulations().isEmpty()) {
//      checkCreatePath(triangulationPath);
      checkCreatePath(triangleShapesPath);
//      checkCreatePath(convexHullPath);
      checkCreatePath(sharedStationsPath);
//      checkCreatePath(candidatesPath);
//      checkCreatePath(edgePath);

      int j = 0;
      int n = 2;
      for(int i = 0; i < results.getTriangulations().size(); i++) {       
        Writer.writeTrianglesToWKT(triangleShapesPath + "/" + results.getNames().get(i) + ".csv", results.getTriangulations().get(i));
//        Writer.writeTriangulation(triangulationPath + "/" + results.getNames().get(i) + ".t3", results.getTriangulations().get(i));
     
        if(i == j) {
       // edgeMap of triangles in triangulation
//        HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap = Triangulator.getEdgeMap(results.getTriangulations().get(i).getTriangles());
//        ConvexHull ch = new ConvexHull();
//        ch.setConvexHull(results.getTriangulations().get(j).getConvexHull());
//        Writer.writeEdgeMap(edgePath + "/" + results.getNames().get(j) + "_valid.csv", edgeMap, ch.compConvexHullEdges());      

//          Writer.writePointsToCSV(convexHullPath + "/" + results.getNames().get(j) + ".csv",
//              results.getTriangulations().get(j).getConvexHull());
          Writer.writePointsToCSV(sharedStationsPath + "/" + results.getNames().get(j).replace("_ILP", "") + ".csv", results.getTriangulations().get(j).getPoints(), results.getTriangulations().get(j).getHeights());
          j+=n;
        }
      }
    }
  }

  /**
   * checks if a file path exists, creates it if not
   * @param pathString
   * @return
   */
  public static boolean checkCreatePath(String pathString) {
    Path path = Paths.get(pathString);
    if(!Files.exists(path)) {
      new File(pathString).mkdirs();
      return false;
    }
    return true;
  }

}
