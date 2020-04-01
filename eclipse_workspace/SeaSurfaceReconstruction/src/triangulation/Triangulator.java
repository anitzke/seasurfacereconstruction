package triangulation;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.index.strtree.STRtree;

import geometry.ConvexHull;
import geometry.HashableLine;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import processing.Calculations;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Triangulator {

  public static final int startFormulation = 0;
  public static final int endFormulation = 1;
  public static final int startSolve = 2;
  public static final int endSolve = 3;
  private int threadCount;
  private List<Triangle> usedTriangles;
  private List<Triangle> candidates;
  private HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edges;
  private Point2D[] points;
  private double[] heights; 
  private Point2D[][] gridPoints; 
  private double[][] gridHeights;
  private List<TriangulationListener> listeners = new ArrayList<>();

  public Triangulator() {
    usedTriangles = null;
    points = null;
    heights = null;
    gridPoints = null;
    gridHeights = null;
    threadCount = 1;
    //		gridData = null;
  }

  public void setThreadCount(int n) {
    threadCount = n;
  }

  public void addListener(TriangulationListener listener) {
    listeners.add(listener);
  }

  public int getThreadCount() {
    return threadCount;
  }

  public List<Triangle> getUsedTriangles() {
    return usedTriangles;
  }	

  public Point2D[] getPoints() {
    return points;
  }

  public double[] getHeights() {
    return heights;
  }

  public Point2D[][] getGridPoints() {
    return gridPoints;
  }

  public double[][] getGridHeights() {
    return gridHeights;
  }

  //	public double[] getGridData() {
  //		return gridData;
  //	}

  public List<Triangle> getCandidates() {
	return candidates;
  }

  public void setCandidates(List<Triangle> candidates) {
	this.candidates = candidates;
  }
  
  

  public HashMap<HashableLine, ArrayList<LinkedList<Integer>>> getEdges() {
    return edges;
  }

  public void setEdges(
      HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edges) {
    this.edges = edges;
  }

  public void reset() {
    usedTriangles = null;
    points = null;
    heights = null;
    gridPoints = null;
    gridHeights = null;
    threadCount = 1;
    edges = null;
    //		gridData = null;
  }

  /**
   * computes a triangulation with the goal of estimating the surface of grid points as well as possible
   * @param triangulation required data: points with corresponding heights, grid points with corresponding heights
   * @return a list of Triangles (also sets the triangle list of input Triangulation)
   */
  public List<Triangle> triangulateSurface(TriangulationInstance triangulation, TriangleCreator tc){

   if(triangulation.getGridPoints() == null) {
      System.out.println("Grid Data needed!");
      return null;
    }

    points = triangulation.getPoints();
    heights = triangulation.getHeights();
    gridPoints = triangulation.getGridPoints();
    gridHeights = triangulation.getGridHeights();

    for(TriangulationListener l : listeners) {
      l.onEvent(startFormulation);
    }
    long startForm = System.currentTimeMillis();

    // point Tree of gridPoints
    STRtree pointTree = getPointTree();

    // HashMap of grid points and the corresponding height
    HashMap<Point2D, Double> gridMap = getGridMap();

    // generate every possible triangle for the given set of points 
    // either 3D if heights are given or 2D if not
    List<Triangle> triangleList = tc.create();
 
    //compute convex hull of given points
    ConvexHull ch = new ConvexHull(points);
    triangulation.setConvexHull(ch.getConvexHull());
    List<HashableLine> hullEdges = ch.compConvexHullEdges();

    //compute 2 lists of adjacent triangles for each edge of every triangle
    HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap = getEdgeMap(triangleList);

    
    boolean cond = false;
    while(cond == false) {
      triangleList = getValidCandidates(triangleList, edgeMap, hullEdges);
      edgeMap = getEdgeMap(triangleList);
      cond = containsInvalidEdges(edgeMap, hullEdges);
    }
    candidates = new ArrayList<Triangle>(triangleList);
    edges = new HashMap<HashableLine, ArrayList<LinkedList<Integer>>>(edgeMap);
    triangulation.setCandidates(candidates);
    triangulation.setEdgeMap(edges);

    // from here: Gurobi Optimizer
    try {
      // set optimization sense
      int sense;
      if (DataCollector.getParameter("-optimization").equals("min")) {
        sense = GRB.MINIMIZE;
      } else {
        sense = GRB.MAXIMIZE;
      }
      
      GRBModel model = buildModel(triangleList, edgeMap, ch.getConvexHull(), hullEdges, pointTree, gridMap, sense);

      model.set("OutputFlag", "0");

      for(TriangulationListener l : listeners) {
        l.onEvent(endFormulation);
      }
      long endForm = System.currentTimeMillis();
      for(TriangulationListener l : listeners) {
        l.onEvent(startSolve);
      }
      long startS = System.currentTimeMillis();
      model.optimize();
      for(TriangulationListener l : listeners) {
        l.onEvent(endSolve);
      }
      long endS = System.currentTimeMillis();

      GRBVar[] vars = model.getVars();
      double[] var = new double[vars.length];
      for(int k = 0; k < vars.length; k++) {
        var[k] = vars[k].get(GRB.DoubleAttr.X);
      }
      // get list of triangles with var = 1 
      List<Triangle> chosenOnes = getChosenTriangles(triangleList, vars);

      model.write("triangulatesurface.lp");

      model.dispose();
      model.getEnv().dispose();
      usedTriangles = chosenOnes;
      triangulation.setTriangles(chosenOnes);
      return chosenOnes;

    } catch (GRBException e) {
      e.printStackTrace();
    }

    return null;
  }

  
  public List<Triangle> triangulateSurface(Point2D[] ps, double[] hs, Point2D[][] gP, double[][] gH, TriangleCreator tc){
    return triangulateSurface(new TriangulationInstance(ps, hs, gP, gH, new ArrayList<Triangle>()), tc);
  }
  
  //	
  public List<Triangle> triangulateDelaunay(TriangulationInstance triangulation, TriangleCreator tc){
    
    if(triangulation.getGridPoints() == null) {
      System.out.println("Grid Data needed!");
      return null;
    }

    points = triangulation.getPoints();
    heights = triangulation.getHeights();
    gridPoints = triangulation.getGridPoints();
    gridHeights = triangulation.getGridHeights();

    for(TriangulationListener l : listeners) {
      l.onEvent(startFormulation);
    }
    long startForm = System.currentTimeMillis();

    // generate every possible triangle for the given set of points 
    // either 3D if heights are given or 2D if not
    List<Triangle> triangleList = tc.create();
 
    //compute convex hull of given points
  //compute convex hull of given points
    ConvexHull ch = new ConvexHull(points);
    triangulation.setConvexHull(ch.getConvexHull());
    List<HashableLine> hullEdges = ch.compConvexHullEdges();

    //compute 2 lists of adjacent triangles for each edge of every triangle
    HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap = getEdgeMap(triangleList);
    
    boolean cond = false;
    while(cond == false) {
      triangleList = getValidCandidates(triangleList, edgeMap, hullEdges);
      edgeMap = getEdgeMap(triangleList);
      cond = containsInvalidEdges(edgeMap, hullEdges);
    }
    candidates = new ArrayList<Triangle>(triangleList);

    // from here: Gurobi Optimizer
    try {
      GRBEnv env = new GRBEnv("triangles.log");
      GRBModel model = new GRBModel(env);

      //create variable for every triangle
      GRBVar[] vars = model.addVars(triangleList.size(), GRB.BINARY);       

      // for each edge select at most one incident triangle on each side 
      addEdgeSumConstraints(model, vars, edgeMap);

      // for each edge of convex hull select exactly 1 incident triangle 
      addConvexHullConstraints(model, vars, edgeMap, ch.getConvexHull());
      
      setWeightDelaunay(model, vars, triangleList);

      model.set("OutputFlag", "0");

      for(TriangulationListener l : listeners) {
        l.onEvent(endFormulation);
      }
      long endForm = System.currentTimeMillis();
      for(TriangulationListener l : listeners) {
        l.onEvent(startSolve);
      }
      long startS = System.currentTimeMillis();
      model.optimize();
      for(TriangulationListener l : listeners) {
        l.onEvent(endSolve);
      }
      long endS = System.currentTimeMillis();

      double[] var = new double[vars.length];
      for(int k = 0; k < vars.length; k++) {
        var[k] = vars[k].get(GRB.DoubleAttr.X);
      }
      // get list of triangles with var = 1 
      List<Triangle> chosenOnes = getChosenTriangles(triangleList, vars);

      model.write("triangulatesurface.lp");

      model.dispose();
      model.getEnv().dispose();
      usedTriangles = chosenOnes;
      triangulation.setTriangles(chosenOnes);
      return chosenOnes;

    } catch (GRBException e) {
      e.printStackTrace();
    }

    return null;
  }

  public GRBModel buildModel(List<Triangle> triangleList, HashMap<HashableLine, 
      ArrayList<LinkedList<Integer>>> edgeMap,
		  List<Point2D> convexHull, 
		  List<HashableLine> hullEdges,
		  STRtree pointTree, 
		  HashMap<Point2D, Double> gridMap,
		  int sense) throws GRBException {
    
    GRBEnv env = new GRBEnv("triangles.log");
    GRBModel model = new GRBModel(env);

    //create variable for every triangle
    GRBVar[] vars = model.addVars(triangleList.size(), GRB.BINARY);

    //======================== new Constraints =========================\\    
    //For each edge, select the same number of incident triangles on each of its sides
    addEdgeEqualityConstraints(model, vars, edgeMap, hullEdges);
    
    // for each edge select at most one incident triangle on each side
    addEdgeSumConstraints(model, vars, edgeMap);

    // for each edge of convex hull select exactly 1 incident triangle 
    addConvexHullConstraints(model, vars, edgeMap, convexHull);

    //minimize the sum of squared distances between triangles and sea level
    setWeightSumObjective(model, vars, triangleList, pointTree, gridMap, sense);
    
    return model;
  }


  public void setWeightDelaunay(GRBModel model, GRBVar[] vars, List<Triangle> triangleList) throws GRBException {
    GRBLinExpr obj = new GRBLinExpr();
    int index = 0;
    for(GRBVar var:vars) {
      obj.addTerm(computeWeight(triangleList.get(index), "leastInnerAngle"), var);
      index++;
    }
    model.setObjective(obj, GRB.MAXIMIZE);
  }
  
  
  public void setWeightSumObjective(GRBModel model, GRBVar[] vars, List<Triangle> triangleList, 
      STRtree pointTree, HashMap<Point2D, Double> gridMap, int sense) throws GRBException {

    int steps = vars.length / threadCount + 1;
    double[][] weights = new double[threadCount][];
    Thread[] threads = new Thread[threadCount];

    for(int i = 0; i < threadCount; i++) {
      final Integer fi = i;
      threads[i] = new Thread(new Runnable() {
        @Override
        public void run() {
          weights[fi] = computeSubSet(fi * steps, (fi+1)*steps, triangleList, pointTree, gridMap);
        }
      });
      threads[i].start();
    }

    //wait for every thread to finish before moving on
    try {
      for(int i = 0; i < threads.length; i++) {
        threads[i].join();
      }
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    //join weights into one array
    List<Double> weightList = new ArrayList<Double>();
    for(int i = 0; i < weights.length; i++) {
      for(int j = 0; j < weights[i].length; j++) {
        weightList.add(weights[i][j]);
      }
    }
    //set objective 
    GRBLinExpr obj = new GRBLinExpr();
    for(int i = 0; i < vars.length; i++) {
      obj.addTerm(weightList.get(i), vars[i]);
    }

    model.setObjective(obj, sense);
  }

  public double[] computeSubSet(int startIdx, int endIdx, List<Triangle> triangleList, STRtree pointTree, HashMap<Point2D, Double> gridMap) {
    if(endIdx > triangleList.size()) {
      endIdx = triangleList.size();
    }
    double[] weights = new double[endIdx-startIdx];

    for(int i = startIdx; i < endIdx; i++) {
      weights[i-startIdx] = computeWeight3D(triangleList.get(i), pointTree, gridMap);
    }

    return weights;
  }

  public void addEdgeSumConstraints(GRBModel model, GRBVar[] vars, HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap) throws GRBException {
    int j = 0;
    for(Map.Entry<HashableLine, ArrayList<LinkedList<Integer>>> e : edgeMap.entrySet()) {
      if(!e.getValue().get(0).isEmpty()) {
        GRBLinExpr leftExpr = new GRBLinExpr();
        for(int k : e.getValue().get(0)) {
          leftExpr.addTerm(1, vars[k]);
        }
        String leftName = "l_" + j + "0";
        model.addConstr(leftExpr, GRB.LESS_EQUAL, 1.0, leftName);
      }
      if(!e.getValue().get(1).isEmpty()) {
        GRBLinExpr rightExpr = new GRBLinExpr();
        for(int k : e.getValue().get(1)) {
          rightExpr.addTerm(1, vars[k]);
        }
        String rightName = "r_" + j + "1";
        model.addConstr(rightExpr, GRB.LESS_EQUAL, 1.0, rightName);
      }
      j++;
    }
  }

  public void addEdgeEqualityConstraints(GRBModel model, GRBVar[] vars, 
      HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap, 
      List<HashableLine> hullEdges) throws GRBException {

    int j = 0;
    for(Map.Entry<HashableLine, ArrayList<LinkedList<Integer>>> e : edgeMap.entrySet()) {
      if(!e.getValue().get(0).isEmpty() && !e.getValue().get(1).isEmpty()) {
        // left side of e
        GRBLinExpr leftExpr = new GRBLinExpr(); // linear expression object
        for(int k : e.getValue().get(0)) { // get triangle ids being on the right side of e
          leftExpr.addTerm(1, vars[k]);
        }
        // right side of e
        GRBLinExpr rightExpr = new GRBLinExpr();
        for(int k : e.getValue().get(1)) { // get triangle ids being on the left side of e
          rightExpr.addTerm(1, vars[k]);
        }
        String name = "eq_" + j;
        model.addConstr(leftExpr, GRB.EQUAL, rightExpr, name);
        //      }
      } else if(e.getValue().get(0).isEmpty() & !hullEdges.contains(e.getKey())) {
        // right side of e
        GRBLinExpr rightExpr = new GRBLinExpr();
        for(int k : e.getValue().get(1)) { // get triangle ids being on the left side of e
          rightExpr.addTerm(1, vars[k]);
        }
        String name = "eq_" + j;
        model.addConstr(rightExpr, GRB.EQUAL, 0.0, name);
      } else if(e.getValue().get(1).isEmpty() & !hullEdges.contains(e.getKey())) {
        // left side of e
        GRBLinExpr leftExpr = new GRBLinExpr(); // linear expression object
        for(int k : e.getValue().get(0)) { // get triangle ids being on the right side of e
          leftExpr.addTerm(1, vars[k]);
        }
        String name = "eq_" + j;
        model.addConstr(leftExpr, GRB.EQUAL, 0.0, name);
      }
      j++;
    }
  }

  public void addConvexHullConstraints(GRBModel model, GRBVar[] vars, HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap, List<Point2D> convexHull) throws GRBException {
    for(int h = 0; h < convexHull.size()-1; h++) {
      GRBLinExpr Expr = new GRBLinExpr();
      for(int f : edgeMap.get(new HashableLine(convexHull.get(h), convexHull.get(h+1))).get(0)) {
        Expr.addTerm(1, vars[f]);
      }
      for(int f : edgeMap.get(new HashableLine(convexHull.get(h), convexHull.get(h+1))).get(1)) {
        Expr.addTerm(1, vars[f]);
      }
      String name = "o_" + h;
      model.addConstr(Expr, GRB.EQUAL, 1.0, name);
    }
  }

  public STRtree getPointTree() {
    STRtree pointTree = new STRtree();
    for(int i = 0; i < gridPoints.length; i++) {
      for(int j = 0; j < gridPoints[0].length; j++) {
        Envelope ev = new Envelope(gridPoints[i][j].getX(), gridPoints[i][j].getX(), gridPoints[i][j].getY(), gridPoints[i][j].getY());
        pointTree.insert(ev, gridPoints[i][j]);
      }
    }
    return pointTree;
  }

  public HashMap<Point2D, Double> getGridMap() {
    HashMap<Point2D, Double> gridMap = new HashMap<Point2D, Double>();
    for(int i = 0; i < gridPoints.length; i++) {
      for(int j = 0; j < gridPoints[0].length; j++) {
        gridMap.put(gridPoints[i][j], gridHeights[i][j]);
      }
    }
    return gridMap;
  }

  public static HashMap<HashableLine, ArrayList<LinkedList<Integer>>> getEdgeMap(List<Triangle> triangleList){
    HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap = 
        new HashMap<HashableLine, ArrayList<LinkedList<Integer>>>();
    int i = 0;
    for(Triangle t : triangleList) {
      HashableLine AB = new HashableLine(t.getA(), t.getB());
      HashableLine BC = new HashableLine(t.getB(), t.getC());
      HashableLine CA = new HashableLine(t.getC(), t.getA());

      //make sure every edge has the same orientation (to avoid duplicates)
      checkTurnEdge(AB);
      checkTurnEdge(BC);
      checkTurnEdge(CA);

      //determine on which side the triangle t lies for each of its edges
      int orientation1 = (Calculations.direction(AB.getP1(), AB.getP2(), t.getC()) < 0) ? 0 : 1; 
      if(edgeMap.get(AB) == null) {
        ArrayList<LinkedList<Integer>> lrTriangles = new ArrayList<LinkedList<Integer>>();
        lrTriangles.add(new LinkedList<Integer>());
        lrTriangles.add(new LinkedList<Integer>());
        lrTriangles.get(orientation1).add(i);
        edgeMap.put(AB, lrTriangles);
      }else {
        edgeMap.get(AB).get(orientation1).add(i);
      }		
      int orientation2 = (Calculations.direction(BC.getP1(), BC.getP2(), t.getA()) < 0) ? 0 : 1; 
      if(edgeMap.get(BC) == null) {
        ArrayList<LinkedList<Integer>> lrTriangles = new ArrayList<LinkedList<Integer>>();
        lrTriangles.add(new LinkedList<Integer>());
        lrTriangles.add(new LinkedList<Integer>());
        lrTriangles.get(orientation2).add(i);
        edgeMap.put(BC, lrTriangles);
      }else {
        edgeMap.get(BC).get(orientation2).add(i);
      }
      int orientation3 = (Calculations.direction(CA.getP1(), CA.getP2(), t.getB()) < 0) ? 0 : 1; 
      if(edgeMap.get(CA) == null) {
        ArrayList<LinkedList<Integer>> lrTriangles = new ArrayList<LinkedList<Integer>>();
        lrTriangles.add(new LinkedList<Integer>());
        lrTriangles.add(new LinkedList<Integer>());
        lrTriangles.get(orientation3).add(i);
        edgeMap.put(CA, lrTriangles);
      }else {
        edgeMap.get(CA).get(orientation3).add(i);
      }
      i++;
    }
    return edgeMap;
  }
  
  public boolean containsInvalidEdges(HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap,
      List<HashableLine> hullEdges) {

    for(Map.Entry<HashableLine, ArrayList<LinkedList<Integer>>> e : edgeMap.entrySet()) {
      ArrayList<LinkedList<Integer>> val = e.getValue();
      if((val.get(0).isEmpty() || val.get(1).isEmpty()) && !hullEdges.contains(e.getKey())) {
        return false;
      }
    }
    return true;
  }
  
  /**
   * returns all valid candidates for the triangulation
   * @param triangleList
   * @param edgeMap
   * @return
   */
  public List<Triangle> getValidCandidates(List<Triangle> triangleList, 
      HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap, 
      List<HashableLine> hullEdges ) {
        
    List<Triangle> validCandidates = new ArrayList<Triangle>(triangleList);
    List<Triangle> collection = new ArrayList<Triangle>();
    Set<Integer> treeSet = new TreeSet<>();

     
    for(Map.Entry<HashableLine, ArrayList<LinkedList<Integer>>> e : edgeMap.entrySet()) {
      if((!e.getValue().get(0).isEmpty() && !e.getValue().get(1).isEmpty()) || hullEdges.contains(e.getKey())) {
        continue;
      } else if(e.getValue().get(0).isEmpty() & !hullEdges.contains(e.getKey())) {
        // right side of e
        for(int k : e.getValue().get(1)) { // get triangle ids being on the left side of e
          treeSet.add(k);
        }
      } else if(e.getValue().get(1).isEmpty() & !hullEdges.contains(e.getKey())) {
        for(int k : e.getValue().get(0)) { // get triangle ids being on the right side of e
          treeSet.add(k);
        }
      }
    }
    
    for(int i : treeSet) {
      collection.add(triangleList.get(i));
    }

    validCandidates.removeAll(collection);
    
    return validCandidates;
  }

  public static Double computeWeight(Object o, String modeName) {
    int mode = -1;
    if(modeName.equals("squaredInnerAngle")) {
      mode = 0;
    }
    if(modeName.equals("leastInnerAngle")) {
      mode = 1;
    }
    switch(mode) {
    case 0: 
      if((o instanceof Triangle)) {
        Triangle t = (Triangle)o; 

        double alpha = Calculations.getInnerAngle(t.getB(), t.getA(), t.getC());

        alpha *= alpha;

        double beta = Calculations.getInnerAngle(t.getA(), t.getB(), t.getC());

        beta *= beta;

        double gamma = Calculations.getInnerAngle(t.getA(), t.getC(), t.getB());

        gamma *= gamma;

        return alpha+beta+gamma;
      }else {
        System.out.println("invalid Object for Mode: squaredInnerAngle");
        break;
      }
    case 1:
      if((o instanceof Triangle)) {
        Triangle t = (Triangle)o; 

        double alpha = Calculations.getInnerAngle(t.getB(), t.getA(), t.getC());

        double beta = Calculations.getInnerAngle(t.getA(), t.getB(), t.getC());

        double gamma = Calculations.getInnerAngle(t.getA(), t.getC(), t.getB());

        return Math.min(alpha, Math.min(beta, gamma));
      }else {
        System.out.println("invalid Object for Mode: leastInnerAngle");
        break;
      }
    default: System.out.println("invalid Mode");
    }
    return null;
  }

  public double computeWeight3D(Triangle t, STRtree pointTree, HashMap<Point2D, Double> gridMap) {
    double weight = 0;
    double[] planeParams = t.planeParameters();

    Envelope env = t.getEnvelope();
    double minX = env.getMinX();
    double minY = env.getMinY();
    double maxX = env.getMaxX();
    double maxY = env.getMaxY();

    double[] gridData = Calculations.getGridData(gridPoints);

    int minXind = (int) Math.ceil((minX - gridData[0])/gridData[2]);
    int maxXind = (int) Math.floor((maxX - gridData[0])/gridData[2]);
    int minYind = (int) Math.ceil((minY - gridData[1])/gridData[3]);
    int maxYind = (int) Math.floor((maxY - gridData[1])/gridData[3]);

    for(int i = minXind; i <= maxXind; i++) {
      for(int j = minYind; j <= maxYind; j++) {
        Point2D p = gridPoints[i][j];
        
     // evaluate all grid points lying in t
        if(t.contains(p) && !Double.isNaN(gridMap.get(p))) {
          // height interpolated from triangle t
          double zT = -(planeParams[0]*p.getX()
              +planeParams[1]*p.getY()
              +planeParams[3])/planeParams[2];
          
          // absolute difference between interpolated and reference height
          double diff = zT - gridMap.get(p);
          
          // form the quadratic sum and sum up all quadratic residuals within t
          weight += diff*diff; // cost function for t
        }
      }
    }

    //return the squared sum of all height differences		
    return weight;
  }

  
  public List<Triangle> getChosenTriangles(List<Triangle> triangleList, GRBVar[] vars) throws GRBException {
    ArrayList<Triangle> chosenOnes = new ArrayList<Triangle>();

    for(int k = 0; k < vars.length; k++) {
      if(vars[k].get(GRB.DoubleAttr.X) > 0.5) {
        chosenOnes.add(triangleList.get(k));
      }
    }
    return chosenOnes;
  }
  
  public static boolean checkTurnEdge(HashableLine line) {
    if(line.getP1().getX() < line.getP2().getX() || line.getP1().getX() == line.getP2().getX() && line.getP1().getY() < line.getP2().getY()) {
      return false;
    }
    line.turn();
    return true;
  }
  
}