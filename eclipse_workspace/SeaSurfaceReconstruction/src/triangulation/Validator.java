package triangulation;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import geometry.ConvexHull;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Validator {

  public static HashMap<String, Integer> validateTriangulation(TriangulationInstance triangulation, List<String> validationLog, boolean detailed) {

    boolean log = false;
    if(validationLog != null) {
      log = true;
    }

    // check number of triangles
    ConvexHull ch = new ConvexHull(triangulation.getPoints());
    List<Point2D> convexHull = ch.getConvexHull();
    int nTrianglesExpected = triangulation.getPoints().length*2-(convexHull.size()-1)-2;
    int nTriangles = triangulation.getTriangles().size();
    boolean nTest  = nTrianglesExpected == nTriangles;

    if(!nTest && log) {
      validationLog.add(triangulation.getName() + ": Expected number of triangles " + nTrianglesExpected + ", actual number of triangles: " + nTriangles);
    }

    // check area of triangles
    double convexArea = ch.area();
    double triangleArea = triangulation.area();

    boolean areaTest = Math.abs(triangleArea-convexArea) < 0.1;

    if(!areaTest && log) {
      validationLog.add(triangulation.getName() + ": Area of convexHull: " + convexArea + ", summarized area of triangles: " + triangleArea);
    }

    // check use of points
    HashMap<Point2D, Integer> pointMap = new HashMap<Point2D, Integer>();
    for(int i = 0; i < triangulation.getPoints().length; i++) {
      pointMap.put(triangulation.getPoints()[i], i);
    }	
    
    boolean[] usedPoints = new boolean[triangulation.getPoints().length];
    for(Triangle triangle : triangulation.getTriangles()) {
      Point2D A = triangle.getA();
      Point2D B = triangle.getB();
      Point2D C = triangle.getC();
      usedPoints[pointMap.get(A)] = true;
      usedPoints[pointMap.get(B)] = true;
      usedPoints[pointMap.get(C)] = true;
    }

    int nUsedPoints = 0;
    List<Point2D> notUsedPoints = new ArrayList<Point2D>();
    int m = 0;
    for (boolean b : usedPoints) {
      nUsedPoints += b ? 1 : 0;
      if (!b) {
        for (HashMap.Entry<Point2D, Integer> entry : pointMap.entrySet()) {
          if (entry.getValue().equals(m)) {
            notUsedPoints.add(entry.getKey());
            break;
          }
        }
      }
      m++;
    }

    boolean pointTest = nUsedPoints == usedPoints.length;

    if(!pointTest && log) {
      validationLog.add(triangulation.getName() +  ": Number of points: " + triangulation.getPoints().length + ", number of used points: " + nUsedPoints);
      validationLog.add(triangulation.getName() +  ": Not used Points: " + notUsedPoints.toString());
    }

    // check overlap of triangles
    List<int[]> overlappingTriangles = new ArrayList<int[]>();
    for(int i = 0; i < triangulation.getTriangles().size(); i++) {
      for(int j = i; j < triangulation.getTriangles().size(); j++) {
        if(triangulation.getTriangles().get(i).intersects(triangulation.getTriangles().get(j))) {
          overlappingTriangles.add(new int[] {i, j});
        }
      }
    }
    boolean overlapTest = overlappingTriangles.isEmpty();

    if(!overlapTest && log) {
      String temp =  triangulation.getName() + ": following triangles intersect: \n";
      for(int[] oT:overlappingTriangles) {
        temp += triangulation.getTriangles().get(oT[0]) + ", " + triangulation.getTriangles().get(oT[1]) + "\n";
      }
      validationLog.add(temp);
    }
    
    HashMap<String, Integer> results = new HashMap<String, Integer>();
    results.put("tests", nTest && areaTest && pointTest && overlapTest ? 1 : 0);
    if(detailed) {
      int triangleDiff = Math.abs(nTrianglesExpected - nTriangles);
      int areaDiff = (int) Math.abs((triangleArea-convexArea));
      int pointDiff = Math.abs(usedPoints.length - nUsedPoints);

      results.put("triangleDiff", triangleDiff);
      results.put("areaDiff", areaDiff);
      results.put("pointDiff", pointDiff);     
    }
    
    return results;      
  }
}