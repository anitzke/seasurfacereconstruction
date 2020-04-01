package triangulation;

import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.index.kdtree.KdNode;
import com.vividsolutions.jts.index.kdtree.KdTree;

import evaluation.TriangulationEvaluator;
import processing.DataHandler;
import processing.DataSet;
import processing.Projection;
import shapes3D.Point3D;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public interface TriangulationWorker{
  /**
   * 
   * Does its work depending on Instance.
   * 
   * Instance Triangulate: triangulates one set of points on one raster.
   * 
   * Instance Evaluate: triangulates any number of point sets and evaluates the
   * triangulations on the rasters for every point set.
   * 
   * Instance Validate: validates any number of triangulations
   * 
   * @return ResultData
   * 
   */
  public ResultData work();
  /**
   * 
   * Simplest variant of worker: Triangulate the point set of each month
   * 
   * @author jonas
   *
   * 
   * 
   */
  public static class Triangulate implements TriangulationWorker {
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
      super();
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
      super();
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
        Path2D pointBounds, TriangleCreator tc, Projection proj) {
      super();
      this.points = points;
      this.identifiers = ids;
      this.epochI = idsI;
      this.epochJ = idsJ;
      this.rasters = rasters;
      this.pointBounds = pointBounds;
      this.tc = tc;
      this.proj = proj;
    }

    public Triangulate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI,  List<String> idsJ, List<Point3D[][]> rasters, 
        Path2D pointBounds, Path2D rasterBounds, TriangleCreator tc) {
      super();
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
    @Override
    public ResultData work() {
      ResultData results = new ResultData();
      TriangulationTimer timer = new TriangulationTimer();
      List<String> nameList = new ArrayList<String>();
      List<double[][]> anomalyList = new ArrayList<double[][]>();
      List<double[][]> squaredAnomalyList = new ArrayList<double[][]>();
      List<String> statisticList = new ArrayList<String>();
      List<TriangulationInstance> triangulationList = new ArrayList<TriangulationInstance>();
      List<String> validationList = new ArrayList<String>();
      List<String> timeList = new ArrayList<String>();

      for(int i = 0; i < epochI.size(); i++) {

        //create and fill triangulations
        TriangulationInstance objectiveTriangulation = new TriangulationInstance();
        TriangulationInstance delaunayTriangulation = new TriangulationInstance();

        int ii = identifiers.indexOf(epochI.get(i));

        objectiveTriangulation.setName("objectiveTriangulation_" + epochI.get(i));
        delaunayTriangulation.setName("delaunayTriangulation_" + epochI.get(i));

        Point3D[][] raster = new Point3D[0][];
        if(rasters.size() > 1) {
          raster = rasters.get(ii);
        }else {
          raster = rasters.get(0);
        }

        DataHandler.setTriangulation(objectiveTriangulation, points.get(ii), raster);
        DataHandler.setTriangulation(delaunayTriangulation, points.get(ii), raster);

        tc.setPoints(objectiveTriangulation.getPoints());
        tc.setHeights(objectiveTriangulation.getHeights());

        //triangulate
        Triangulator t = new triangulation.Triangulator();
        t.addListener(new TriangulationListener() {
          @Override
          public void onEvent(int status) {
            switch(status) {
            case Triangulator.startFormulation:
              timer.recordFormulationStart();
              break;
            case Triangulator.endFormulation:
              timer.recordFormulationEnd();
              break;
            case Triangulator.startSolve:
              timer.recordSolveStart();
              break;
            case Triangulator.endSolve:
              timer.recordSolveEnd();
              break;
            default:
              break;
            }
          }
        });
        t.triangulateSurface(objectiveTriangulation, tc);
        //        t.triangulateHeights(objectiveTriangulation, null, tc);
        t.reset();
        t.triangulateDelaunay(delaunayTriangulation, tc);

        if(delaunayTriangulation.getTriangles().isEmpty() || objectiveTriangulation.getTriangles().isEmpty()) {
          try {
            throw new IOException("Empty Triangulation! in " + identifiers.get(ii));
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
        }

        //store results
        triangulationList.add(objectiveTriangulation);
        nameList.add(objectiveTriangulation.getName());
        triangulationList.add(delaunayTriangulation);
        nameList.add(delaunayTriangulation.getName());

        //evaluate objective triangulation
        TriangulationEvaluator bla = new TriangulationEvaluator();
        double[] objStats = bla.evaluate(objectiveTriangulation);
        anomalyList.add(TriangulationEvaluator.getAnomalies());
        squaredAnomalyList.add(TriangulationEvaluator.getSquaredAnomalies());

        //set heights of delaunay triangulation for evaluation 
        //Evaluate e = new Evaluate();
        //e.changeTriangleHeights(delaunayTriangulation, objectiveTriangulation);

        //evaluate delaunay triangulation
        double[] delStats = bla.evaluate(delaunayTriangulation);
        anomalyList.add(TriangulationEvaluator.getAnomalies());
        squaredAnomalyList.add(TriangulationEvaluator.getSquaredAnomalies());

        statisticList.add(DataHandler.createStatString(objectiveTriangulation.getName(), objStats));
        statisticList.add(DataHandler.createStatString(delaunayTriangulation.getName(), delStats));

        Validate.validate(objectiveTriangulation, validationList);
        Validate.validate(delaunayTriangulation, validationList); 

        timeList.add(DataHandler.createTimeString(objectiveTriangulation.getName(), timer.getTimes()));
      }
      results.addNames(nameList);
      results.addAnomalies(anomalyList);
      results.addSquaredAnomalies(squaredAnomalyList);
      results.addStatistics(statisticList);
      results.addTriangulations(triangulationList);
      results.addValidations(validationList);
      results.addTimes(timeList);

      return results;
    }

    /**
     * 
     * changes the point heights original to the point heights of
     * onTriangulation
     * 
     * @param original
     * @param onTriangulation
     * 
     */
    public void changeTriangleHeights(TriangulationInstance original,
        TriangulationInstance onTriangulation) {
      HashMap<Point2D, Double> pointMap = new HashMap<Point2D, Double>();
      for (int i = 0; i < onTriangulation.getPoints().length; i++) {
        pointMap.put(onTriangulation.getPoints()[i],
            onTriangulation.getHeights()[i]);
      }
      
      double[] heights = new double[original.getHeights().length];
      for (int i = 0; i < original.getPoints().length; i++) {
        heights[i] = pointMap.get(original.getPoints()[i]);
      }
      
      List<Triangle> newTriangles = new ArrayList<Triangle>();
      for (Triangle t : original.getTriangles()) {
        newTriangles.add(new Triangle(t.getA(), pointMap.get(t.getA()), 
            t.getB(), pointMap.get(t.getB()),
            t.getC(), pointMap.get(t.getC())));
      }
      
      original.setHeights(heights);
      original.setTriangles(newTriangles);
    }
   

  }



  /**
   * 
   * Computes for each pair (A,B) of months the intersection of their tide
   * gauges and computes the following triangulations:
   * 
   * - for month A we are given tide gauges
   * - for month B we are given the altimeter data
   * 
   * The method triangulates A on its tide gauges obtaining a 3-dimensional
   * triangulation T. The geometry of T is transferred to the tide gauges of 
   * B obtaining a new 3-dimensional triangulation T'.
   * Further, a delaunay triangulation D. The triangulations T' and D are
   *  evaluated on the altimeter data of B
   * 
   */
  public static class Evaluate extends Triangulate {
    private TriangulationInstance singleEvalTriangulation;
    private ResultData results;
    private TriangulationTimer timer;
    private List<String> nameList;
    private List<double[][]> anomalyList;
    private List<double[][]> squaredAnomalyList;
    private List<double[][]> interpolationsList;
    private List<String> statisticList;
    private List<TriangulationInstance> triangulationList;
    private List<String> validationList;
    private List<String> timeList;
    private List<double[][]> sqaDiffList;

    public Evaluate() {
    }

    public Evaluate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters,
        TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, tc);
    }

    public Evaluate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters,
        Path2D pointBounds, TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, pointBounds, tc);
    }

    public Evaluate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, 
        Path2D pointBounds, Path2D rasterBounds, TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, pointBounds, rasterBounds, 
          tc);
    }

    public Evaluate(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, 
        Path2D pointBounds, Path2D rasterBounds, TriangleCreator tc,
        TriangulationInstance singleEvalTriangulation) {
      super(points, ids, idsI, idsJ, rasters, pointBounds, rasterBounds, 
          tc);
      this.singleEvalTriangulation = singleEvalTriangulation;
    }

    public void setSingleEvaluationTriangulation(TriangulationInstance t) {
      singleEvalTriangulation = t;
    }

    public void setProj(Projection proj) {
      super.setProj(proj);
    }

    /**
     * 
     * ### DOES THE WORK
     * 
     */
    @Override
    public ResultData work() {
      statisticList = new ArrayList<String>();
      validationList = new ArrayList<String>();
      timeList = new ArrayList<String>();

      // toEval: the triangulation that will be evaluated on multiple rasters
      TriangulationInstance toEval = new TriangulationInstance();
      if (singleEvalTriangulation != null) { // single month is evaluated
        toEval = singleEvalTriangulation;
        evaluateOneTriangulation(toEval, identifiers.indexOf(epochI.get(0)));
      } else {
        // each month is evaluated on all given months (276*276 runs) or a
        // special month
        for (int i = 0; i < epochI.size(); i++) {
          toEval.setName(epochI.get(i));
          int ii = identifiers.indexOf(epochI.get(i));
          DataHandler.setTriangulation(toEval, points.get(ii), rasters.get(ii));
          evaluateOneTriangulation(toEval, ii);
          Timestamp timestamp = new Timestamp(System.currentTimeMillis());
          System.out.println(
              timestamp + "     Evaluation of " + toEval.getName() + " ready");
        }
      }
      return results;
    }
    /**
     * 
     * evaluates the input triangulation for another time stemp
     * 
     * @param toEval
     * 
     */
    private void evaluateOneTriangulation(TriangulationInstance toEval, int idEpochI) {
      results = new ResultData();
      timer = new TriangulationTimer();
      nameList = new ArrayList<String>();
      anomalyList = new ArrayList<double[][]>();
      squaredAnomalyList = new ArrayList<double[][]>();
      interpolationsList = new ArrayList<double[][]>();
      triangulationList = new ArrayList<TriangulationInstance>();
      sqaDiffList = new ArrayList<double[][]>();

      for(int i=0; i<epochJ.size(); i++) {        
        int ii = identifiers.indexOf(epochJ.get(i));
//        if(ii == identifiers.indexOf(DataCollector.getParameter("-ev"))) {
//          continue;
//        }
        TriangulationInstance basis = new TriangulationInstance();
        basis.setName(identifiers.get(ii));
        DataHandler.setTriangulation(basis, points.get(ii), rasters.get(ii));
        evaluateTriangulationPair(toEval, basis);
      }

      results.addNames(nameList);
      results.addAnomalies(anomalyList);
      results.addSquaredAnomalies(squaredAnomalyList);
      results.addInterpolatedValues(interpolationsList);
      results.addStatistics(statisticList);
      results.addTriangulations(triangulationList);
      results.addValidations(validationList);
      results.addTimes(timeList);
      results.addSqaDiffs(sqaDiffList);
      DataDistributor distributor = new DataDistributor(results,
          DataCollector.getParameter("-store"));
      distributor.distribute();
    }

    /**
     * 
     * evaluates toEval on the raster of basisTriangulation
     * 
     * @param toEval 
     * @param basisTriangulation
     * 
     */
    private void evaluateTriangulationPair(TriangulationInstance toEval, TriangulationInstance basisTriangulation) {
      // get tide gauges present in both triangulations
      Point3D[] sharedPoints = DataHandler.getSharedPoints(toEval, basisTriangulation);

      // create Triangulations and prepare for ilp and delaunay triangulation
      TriangulationInstance sharedTriangulationILP = new TriangulationInstance();
      TriangulationInstance sharedTriangulationDel = new TriangulationInstance();
      sharedTriangulationILP.setName(toEval.getName() + "_on_" + basisTriangulation.getName() + "_ILP");
      sharedTriangulationDel.setName(toEval.getName() + "_on_" + basisTriangulation.getName() + "_Delaunay");
      DataHandler.setTriangulation(sharedTriangulationILP, sharedPoints, null); // without raster
      DataHandler.setTriangulation(sharedTriangulationDel, sharedPoints, null); // without raster
      sharedTriangulationILP.setGrid(toEval.getGridPoints(), toEval.getGridHeights()); // set grid points and heights for ILP
      sharedTriangulationDel.setGrid(toEval.getGridPoints(), toEval.getGridHeights()); // set grid points and heights for ILP
      tc.setPoints(sharedTriangulationILP.getPoints());
      tc.setHeights(sharedTriangulationILP.getHeights());

      // triangulate common points on initial raster
      Triangulator t = new Triangulator();
      t.addListener(new TriangulationListener() {
        @Override
        public void onEvent(int status) {
          switch (status) {
          case Triangulator.startFormulation:
            timer.recordFormulationStart();
            break;
          case Triangulator.endFormulation:
            timer.recordFormulationEnd();
            break;
          case Triangulator.startSolve:
            timer.recordSolveStart();
            break;
          case Triangulator.endSolve:
            timer.recordSolveEnd();
            break;
          default:
            break;
          }
        }
      });

      List<Triangle> emptyTriangles = DataSet.getEmptyTrianglesSweep(
          sharedTriangulationDel.getPoints(),
          sharedTriangulationDel.getHeights());

      // Delaunay triangulation
      TriangleCreator tcDel = new TriangleCreator.Constrained(
          sharedTriangulationDel.getPoints(),
          sharedTriangulationDel.getHeights(), 0, emptyTriangles);
      t.triangulateDelaunay(sharedTriangulationDel, tcDel);
      timeList.add(DataHandler.createTimeString(
          sharedTriangulationDel.getName(), timer.getTimes()));

      // Triangulation from time A to time B
      tc.setEmptyTriangles(emptyTriangles);
      t.triangulateSurface(sharedTriangulationILP, tc);
      timeList.add(DataHandler.createTimeString(
          sharedTriangulationILP.getName(), timer.getTimes()));

      // evaluate common triangulations on the raster of basisTriangulation
      sharedTriangulationILP.setGrid(basisTriangulation.getGridPoints(),
          basisTriangulation.getGridHeights());
      changeTriangleHeights(sharedTriangulationILP, basisTriangulation);
      sharedTriangulationDel.setGrid(basisTriangulation.getGridPoints(),
          basisTriangulation.getGridHeights());
      changeTriangleHeights(sharedTriangulationDel, basisTriangulation);

      // local      
      List<TriangulationInstance> triangulations = new ArrayList<TriangulationInstance>();
      triangulations.add(sharedTriangulationILP);
      triangulations.add(sharedTriangulationDel);

      List<double[][]> sqas = new ArrayList<double[][]>();
      // evaluate and validate the three created triangulations
      for (TriangulationInstance T : triangulations) {
        // evaluate
        TriangulationEvaluator te = new TriangulationEvaluator();
        double[] stats = te.evaluate(T);
        anomalyList.add(TriangulationEvaluator.getAnomalies());
        sqas.add(TriangulationEvaluator.getSquaredAnomalies());
        squaredAnomalyList.add(TriangulationEvaluator.getSquaredAnomalies());
        statisticList.add(DataHandler.createStatString(T.getName(), stats));
        interpolationsList.add(TriangulationEvaluator.getinterpolatedValues());

        // validate
        Validate.validate(T, validationList);

        // project back to polar coords
        if (DataCollector.getParameter("-projection").equals("LCC")) {
          Projection proj = DataCollector.getProjection();
          proj.setDirection("");
          T.project(proj);
        }

        // save global
        triangulationList.add(T);
        nameList.add(T.getName());
      }
      
      // squared anomalies difference
      double[][] sqa_diff = new double[sqas.get(0).length][sqas.get(0)[0].length];
      for(int i = 0; i < sqas.get(0).length; i++) {
        for(int j = 0; j < sqas.get(0)[0].length; j++) {
          sqa_diff[i][j] = sqas.get(0)[i][j] - sqas.get(1)[i][j];
        }
      }
      
      sqaDiffList.add(sqa_diff);
      
    }
    
  }
  

  public static class Validate implements TriangulationWorker {
    public List<TriangulationInstance> triangulations;
    public List<String> identifiers;
    public Validate() {
    }
    public Validate(List<TriangulationInstance> triangulations) {
      this.triangulations = triangulations;
    }
    public Validate(List<TriangulationInstance> triangulations, List<String> ids) {
      this.triangulations = triangulations;
      this.identifiers = ids;
    }
    @Override
    public ResultData work() {
      ResultData results = new ResultData();
      List<String> nameList = new ArrayList<String>();
      List<String> validationList = new ArrayList<String>();
      for (int i = 0; i < triangulations.size(); i++) {
        triangulations.get(i).setName(identifiers.get(i));
        nameList.add(triangulations.get(i).getName());
        validate(triangulations.get(i), validationList);
      }
      results.addNames(nameList);
      results.addValidations(validationList);
      return results;
    }
    public static void validate(TriangulationInstance t, List<String> validationLog) {
      HashMap<String, Integer> validation = Validator.validateTriangulation(t,
          validationLog, false);
      if (validation.get("tests") == 1) {
        // validationLog.add(t.getName() + ": OK");
      } else {
        validationLog.add(t.getName() + ": not valid");
      }
    }
  }

  public static class Reconstruct extends Triangulate {
    private ResultData results;
    private TriangulationTimer timer;
    private List<String> nameList;
    private List<double[][]> anomalyList;
    private List<double[][]> squaredAnomalyList;
    private List<double[][]> reconstructionsList;
    private List<String> statisticList;
    private List<TriangulationInstance> triangulationList;
    private List<String> validationList;
    private List<String> timeList;
    private List<String> recMeanList;

    public Reconstruct() {
    }

    public Reconstruct(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, 
        TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, tc);
    }

    public Reconstruct(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, 
        Path2D pointBounds, TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, pointBounds, tc);
    }

    public Reconstruct(List<Point3D[]> points, List<String> ids, 
        List<String> idsI, List<String> idsJ, List<Point3D[][]> rasters, 
        Path2D pointBounds, Path2D rasterBounds, TriangleCreator tc) {
      super(points, ids, idsI, idsJ, rasters, pointBounds, rasterBounds, 
          tc);
    }

    public void setProj(Projection proj) {
      super.setProj(proj);
    }
    /**
     * 
     * ### DOES THE WORK
     * 
     */
    @Override
    public ResultData work() {
      results = new ResultData();
      timer = new TriangulationTimer();
      nameList = new ArrayList<String>();
      anomalyList = new ArrayList<double[][]>();
      squaredAnomalyList = new ArrayList<double[][]>();
      statisticList = new ArrayList<String>();
      triangulationList = new ArrayList<TriangulationInstance>();
      validationList = new ArrayList<String>();
      timeList = new ArrayList<String>();
      reconstructionsList = new ArrayList<double[][]>();
      recMeanList = new ArrayList<String>();
      recMeanList.add("epoch,rec_mean_obj,rec_mean_del,tide_gauge_mean");

      Projection proj = DataCollector.getProjection();
      proj.setDirection("");

      for(int i=0; i<epochI.size(); i++) {
        int ii = identifiers.indexOf(epochI.get(i)); // index of epoch i
        int time_delay = Integer.parseInt(DataCollector.getParameter("-rt"));
        for(int j=ii; j>=0; j=j-time_delay) {

          Point3D[] common = DataHandler.getSharedPoints(this.points.get(ii), this.points.get(j));

          // create and fill optimized and delaunay triangulation for comparison
          TriangulationInstance objectiveTriangulation = new TriangulationInstance();
          TriangulationInstance delaunayTriangulation = new TriangulationInstance();
          objectiveTriangulation.setName(identifiers.get(ii) + "_" + identifiers.get(j) + "_o_i");
          delaunayTriangulation.setName(identifiers.get(ii) + "_" + identifiers.get(j) + "_d_i");
          DataHandler.setTriangulation(objectiveTriangulation, common, this.rasters.get(i));
          DataHandler.setTriangulation(delaunayTriangulation, common, this.rasters.get(i));
          this.tc.setPoints(objectiveTriangulation.getPoints());
          this.tc.setHeights(objectiveTriangulation.getHeights());

          // triangulate
          Triangulator t = new Triangulator();
          t.addListener(new TriangulationListener() {
            @Override
            public void onEvent(int status) {
              switch (status) {
              case Triangulator.startFormulation:
                timer.recordFormulationStart();
                break;
              case Triangulator.endFormulation:
                timer.recordFormulationEnd();
                break;
              case Triangulator.startSolve:
                timer.recordSolveStart();
                break;
              case Triangulator.endSolve:
                timer.recordSolveEnd();
                break;
              default:
                break;
              }
            }
          });

          // Compute all empty Triangles
          List<Triangle> emptyTriangles = DataSet.getEmptyTrianglesSweep(
              delaunayTriangulation.getPoints(),
              delaunayTriangulation.getHeights());

          // Delaunay triangulation
          TriangleCreator tcDel = new TriangleCreator.Constrained(
              delaunayTriangulation.getPoints(),
              delaunayTriangulation.getHeights(), 0, emptyTriangles);
          t.triangulateDelaunay(delaunayTriangulation, tcDel);
          timeList.add(DataHandler.createTimeString(
              delaunayTriangulation.getName(), timer.getTimes()));

          // Optimized triangulation
          this.tc.setEmptyTriangles(emptyTriangles);
          t.triangulateSurface(objectiveTriangulation, this.tc);
          timeList.add(DataHandler.createTimeString(
              objectiveTriangulation.getName(), timer.getTimes()));

          // local      
          List<TriangulationInstance> triangulations = new ArrayList<TriangulationInstance>();
          triangulations.add(objectiveTriangulation);
          triangulations.add(delaunayTriangulation);
          
          int m = 0;
          String name = null;
          double[] rec_stats = new double[3];
          for (TriangulationInstance T : triangulations) {
            // evaluate
            TriangulationEvaluator te = new TriangulationEvaluator();
            double[] stats = te.evaluate(T);
            anomalyList.add(TriangulationEvaluator.getAnomalies());
            squaredAnomalyList.add(TriangulationEvaluator.getSquaredAnomalies());
            statisticList.add(DataHandler.createStatString(T.getName(), stats));
            reconstructionsList.add(TriangulationEvaluator.getinterpolatedValues());

            // validate
            Validate.validate(T, validationList);

            // project back to polar coords 
            TriangulationInstance T_copy = new TriangulationInstance();
            try {
              T_copy = (TriangulationInstance) T.clone();
            } catch (CloneNotSupportedException e) {
              // TODO Auto-generated catch block
              e.printStackTrace();
            }
            T_copy.project(proj);

            // save global
            triangulationList.add(T_copy);
            nameList.add(T.getName());
            
            TriangulationInstance eJ = new TriangulationInstance();
            try {
              eJ = (TriangulationInstance) T.clone();
            } catch (CloneNotSupportedException e) {
              // TODO Auto-generated catch block
              e.printStackTrace();
            }
            eJ.setPoints(points.get(j));            

            // change tide gauge heights to epoch j
            changeTriangleHeights(T, eJ); 

            // reconstruct epoch j based on triangulation of epoch i
            double[][] reconstruction = 
                new double[this.rasters.get(i).length][this.rasters.get(i)[0].length];
            double[] rec_mean = new double[1];
            rec_mean[0] = 0;
            reconstruction = reconstruct(T, rec_mean);
            rec_stats[m] = rec_mean[0];
            reconstructionsList.add(reconstruction);

            // save global
            T.project(proj);
            triangulationList.add(T);
            name  = T.getName().split("_")[1];
            nameList.add(T.getName().replace("i", "j"));
            m++;
          }
          
          Point3D[] common_j = DataHandler.getSharedPoints(this.points.get(j), this.points.get(ii));
          double mean_height = 0;
          for(Point3D p : common_j) {
            mean_height = mean_height + p.getZ();
          }
          rec_stats[m] = mean_height/common_j.length;
          recMeanList.add(DataHandler.createRecString(name, rec_stats));
        }
      }

      results.addNames(nameList);
      results.addAnomalies(anomalyList);
      results.addSquaredAnomalies(squaredAnomalyList);
      results.addStatistics(statisticList);
      results.addTriangulations(triangulationList);
      results.addValidations(validationList);
      results.addTimes(timeList);
      results.addInterpolatedValues(reconstructionsList);
      results.addInterpolationMeans(recMeanList);

      DataDistributor distributor = new DataDistributor(results,
          DataCollector.getParameter("-store"));
      distributor.distribute();

      return results;
    }

    public void changeTriangleHeights(TriangulationInstance original,
        Point3D[] points) {
      HashMap<Point2D, Double> pointMap = new HashMap<Point2D, Double>();
      for (int i = 0; i < points.length; i++) {
        pointMap.put(new Point2D.Double(points[i].getX(), points[i].getY()), 
            points[i].getZ());
      }

      double[] heights = new double[original.getHeights().length];
      Point2D[] op = original.getPoints();                  
      for (int i=0; i<op.length; i++) {
        heights[i] = pointMap.get(op[i]);
      }
      original.setHeights(heights);

      List<Triangle> newTriangles = new ArrayList<Triangle>();
      for (Triangle t : original.getTriangles()) {
        newTriangles.add(new Triangle(t.getA(), pointMap.get(t.getA()), 
            t.getB(), pointMap.get(t.getB()),
            t.getC(), pointMap.get(t.getC())));
      }
      
      original.setTriangles(newTriangles);
    }

    public double[][] reconstruct(TriangulationInstance triangulation, double[] rec_mean) {
      Point2D[][] gridPoints = triangulation.getGridPoints();
      double[][] gridHeights = triangulation.getGridHeights();

      double[][] reconstructedValues = new double[gridPoints.length][gridPoints[0].length];
      for(int i=0; i<gridPoints.length; i++) {
        Arrays.fill(reconstructedValues[i], Double .NaN);
      }

      int k = 0;
      KdTree tree = new KdTree();
      HashMap<Point2D, Double> gridMap = new HashMap<Point2D, Double>();
      HashMap<Point2D, ArrayList<Integer>> gridMapIds = new HashMap<Point2D, ArrayList<Integer>>();
      for (int l=0; l<gridPoints.length; l++, k++) {
        int m = 0;
        for (int o=0; o<gridPoints[0].length; o++, m++) {
          Point2D p = gridPoints[l][o];
          ArrayList<Integer> ids = new ArrayList<Integer>();
          ids.add(l);
          ids.add(o);
          tree.insert(new Coordinate(p.getX(), p.getY())); 
          gridMap.put(p, gridHeights[k][m]);
          gridMapIds.put(p, ids);
        }
      }

      double rv_sum = 0;
      double weight_sum = 0;
      for(Triangle tau : triangulation.getTriangles()) {    

        double[] planeParams = tau.planeParameters();
        for (Object o : tree.query(tau.getEnvelope())) {
          KdNode node = (KdNode) o;
          Point2D p = new Point2D.Double(node.getCoordinate().x, 
              node.getCoordinate().y);

          // reconstruct all grid points lying in tau
          if(tau.contains(p) && !Double.isNaN(gridMap.get(p))) {
            // height interpolated from triangle tau
            double zT = -(planeParams[0]*p.getX()
                +planeParams[1]*p.getY()
                +planeParams[3])/planeParams[2];

            int ii = gridMapIds.get(p).get(0);
            int jj = gridMapIds.get(p).get(1);
            
            Projection pro = DataCollector.getProjection();
            proj.setDirection("");
            Point2D p_pro= pro.project(p);

            reconstructedValues[ii][jj] = zT;
            double weight = Math.cos(Math.toRadians(p_pro.getY()));
            rv_sum = rv_sum + zT*weight;
            weight_sum = weight_sum + weight;
          }
        }
      }
      rec_mean[0] = rv_sum/weight_sum;
      return reconstructedValues;
    }
  }
}
