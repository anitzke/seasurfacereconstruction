package evaluation;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.index.kdtree.KdNode;
import com.vividsolutions.jts.index.kdtree.KdTree;

import triangulation.Triangle;
import triangulation.TriangulationInstance;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class TriangulationEvaluator implements Evaluator{

  private static double[][] squaredAnomalies = new double[0][];
  private static double[][] anomalies = new double[0][];
  private static double[][] interpolatedValues = new double[0][];
  
  
  public static double[][] getSquaredAnomalies(){
    return squaredAnomalies;
  }

  public static double[][] getAnomalies(){
    return anomalies;
  }
  
  public static double[][] getinterpolatedValues(){
    return interpolatedValues;
  }

  /**
   * Evaluates a given triangulation with respect to the goal of estimating the surface of given grid points as well as possible
   * @param triangulation 
   */
  @Override
  public double[] evaluate(Object evaluatee) {
    
    TriangulationInstance triangulation = new TriangulationInstance();
    if(evaluatee instanceof TriangulationInstance) {
      triangulation = (TriangulationInstance)evaluatee;
    }

    Point2D[][] gridPoints = triangulation.getGridPoints();
    double[][] gridHeights = triangulation.getGridHeights();
    
    int k = 0;
    KdTree tree = new KdTree();
    HashMap<Point2D, Double> gridMap = new HashMap<Point2D, Double>();
    HashMap<Point2D, ArrayList<Integer>> gridMapIds = new HashMap<Point2D, ArrayList<Integer>>();
    for (int i=0; i<gridPoints.length; i++, k++) {
      int m = 0;
      for (int j=0; j<gridPoints[0].length; j++, m++) {
        Point2D p = gridPoints[i][j];
        ArrayList<Integer> ids = new ArrayList<Integer>();
        ids.add(i);
        ids.add(j);
        tree.insert(new Coordinate(p.getX(), p.getY())); 
        gridMap.put(p, gridHeights[k][m]);
        gridMapIds.put(p, ids);
      }
    }

//    double[] gridData = Calculations.getGridData(gridPoints);
    squaredAnomalies = new double[gridPoints.length][gridPoints[0].length];
    anomalies = new double[gridPoints.length][gridPoints[0].length];
    interpolatedValues = new double[gridPoints.length][gridPoints[0].length];

    for(int i = 0; i < squaredAnomalies.length; i++) {
      Arrays.fill(squaredAnomalies[i], Double.NaN); 
      Arrays.fill(anomalies[i], Double .NaN);
      Arrays.fill(interpolatedValues[i], Double .NaN);
    }
    
//    Projection proj = DataCollector.getProjection();
//    proj.setDirection("");
//    Point2D[][] gridPoints_= proj.project(gridPoints);

    double varianceSum = 0;
    double anomalySum = 0;
    int counter = 0;
    double triangleArea = 0;
    for(Triangle t : triangulation.getTriangles()) {    
//      Triangle t2 = proj.project(t);
//      if((t2.getA().getX()==7.555 |  t2.getA().getX()==8.441 | t2.getA().getX()==-0.615) &
//          (t2.getB().getX()==7.555 |  t2.getB().getX()==8.441 | t2.getB().getX()==-0.615) &
//         ( t2.getC().getX()==7.555 |  t2.getC().getX()==8.441 | t2.getC().getX()==-0.615)) {
//        int bla = 1;
//      }

      triangleArea += t.area();
      double[] planeParams = t.planeParameters();
      for (Object o : tree.query(t.getEnvelope())) {
        KdNode node = (KdNode) o;
        Point2D p = new Point2D.Double(node.getCoordinate().x, 
            node.getCoordinate().y);
        
//        Point2D p_ = proj.project(p);
        
          // evaluate all grid points lying in t
          if(t.contains(p) && !Double.isNaN(gridMap.get(p))) {
            // height interpolated from triangle t
            double zT = -(planeParams[0]*p.getX()
                +planeParams[1]*p.getY()
                +planeParams[3])/planeParams[2];
            
            // absolute difference between interpolated and reference height
            double diff = zT - gridMap.get(p);
      
            varianceSum += diff*diff; // sum of squared differences over all points in t
            anomalySum += diff;       // sum of differences over all points in t
            counter++;
            int ii = gridMapIds.get(p).get(0);
            int jj = gridMapIds.get(p).get(1);
            squaredAnomalies[ii][jj] = diff*diff;
            anomalies[ii][jj] = diff;
            interpolatedValues[ii][jj] = zT;
          }
        }
      }

    double mean = anomalySum/counter;
    double variance = varianceSum/(counter-1);
    double stdDev = Math.sqrt(variance);
   
    double[] stats = {mean, 
        varianceSum, 
        variance, 
        stdDev, 
        triangleArea, 
        triangulation.getPoints().length,
        triangulation.getCandidates().size(),
        };
    
    return stats;
  }

}
