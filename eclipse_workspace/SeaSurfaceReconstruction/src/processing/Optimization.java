package processing;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.index.strtree.STRtree;

import geometry.HashableLine;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import triangulation.Triangle;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Optimization {
	private GRBEnv env;
	private GRBModel model;
	private GRBVar[] vars;

	/* ---------- constructors ---------- */
	public Optimization() {
	}

	public Optimization(GRBEnv env) {
		this.env = env;
	}

	public Optimization(GRBModel model) {
		this.model = model;
	}

	public Optimization(GRBEnv env, GRBModel model) {
		this.env = env;
		this.model = model;
	}

	public Optimization(GRBEnv env, GRBModel model, GRBVar[] vars) {
		this.env = env;
		this.model = model;
		this.vars = vars;
	}

	/* ---------- getter and setter ---------- */
	public GRBEnv getEnv() {
		return env;
	}
	public void setEnv(GRBEnv env) {
		this.env = env;
	}
	public GRBModel getModel() {
		return model;
	}
	public void setModel(GRBModel model) {
		this.model = model;
	}

	public GRBVar[] getVars() {
		return vars;
	}

	public void setVars(GRBVar[] vars) {
		this.vars = vars;
	}

	/* ---------- weighting ---------- */
	/**
	 * @param triangleList
	 * @param pointTree
	 * @param gridMap
	 * @param threadCount
	 * @throws GRBException
	 */
	public void setWeightSumObjective(List<Triangle> triangleList, STRtree pointTree, 
			HashMap<Point2D, Double> gridMap, int threadCount, Point2D[][] gridPoints) throws GRBException {

		int steps = vars.length / threadCount + 1;
		double[][] weights = new double[threadCount][];
		Thread[] threads = new Thread[threadCount];

		for(int i = 0; i < threadCount; i++) {
			final Integer fi = i;
			threads[i] = new Thread(new Runnable() {
				@Override
				public void run() {
					weights[fi] = computeSubSet(fi * steps, (fi+1)*steps, triangleList, pointTree, gridMap, gridPoints);
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

		model.setObjective(obj, GRB.MINIMIZE);
	}

	/**
	 * @param o
	 * @param modeName
	 * @return
	 */
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

	/**
	 * @param t
	 * @param pointTree
	 * @param gridMap
	 * @return
	 */
	public double computeWeight3D(Triangle t, STRtree pointTree, HashMap<Point2D, Double> gridMap, Point2D[][] gridPoints) {
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

				// searching for all grid points that are in t and have a height (not NaN)
				if(t.contains(p) && !Double.isNaN(gridMap.get(p))) {        	
					// compute height of grid point by interpolation in t (by using
					// hessian normal form of the triangle plane
					double zT = -(planeParams[0]*p.getX()
							+planeParams[1]*p.getY()
							+planeParams[3])/planeParams[2];

					// determine residual of height interpolated and measured
					double dist = Math.abs(zT - gridMap.get(p));

					// form the quadratic sum and sum up all quadratic residuals within t
					weight += dist*dist; // cost function for t
				}
			}
		}

		//return the squared sum of all height differences		
		return weight;
	}

	public double[] computeSubSet(int startIdx, int endIdx, List<Triangle> triangleList, STRtree pointTree, 
	    HashMap<Point2D, Double> gridMap, Point2D[][] gridPoints) {
		if(endIdx > triangleList.size()) {
			endIdx = triangleList.size();
		}
		double[] weights = new double[endIdx-startIdx];

		for(int i = startIdx; i < endIdx; i++) {
			weights[i-startIdx] = computeWeight3D(triangleList.get(i), pointTree, gridMap, gridPoints);
		}

		return weights;
	}

	public void setSingleDistanceObjective(GRBModel model, GRBVar[] vars, Point2D p2min, List<Triangle> triangleList, HashMap<Point2D, Double> gridMap) throws GRBException {
		GRBLinExpr obj = new GRBLinExpr();
		for(int i = 0; i < vars.length; i++) {
			double[] planeParams = triangleList.get(i).planeParameters();
			if(triangleList.get(i).contains(p2min)) {
				double zT = -(planeParams[0]*p2min.getX()
						+planeParams[1]*p2min.getY()
						+planeParams[3])/planeParams[2];
				double dist = Math.abs(zT - gridMap.get(p2min)); 
				obj.addTerm(dist, vars[i]);
			}
		}
		model.setObjective(obj, GRB.MINIMIZE);
	}

	/* ---------- constraints ---------- */

	/**
	 * @param model
	 * @param vars
	 * @param p
	 * @param triangleList
	 * @throws GRBException
	 */
	public void addOverlapConstraint(GRBModel model, GRBVar[] vars, Point2D p, List<Triangle> triangleList) throws GRBException {
		GRBLinExpr expr = new GRBLinExpr();
		for(int i = 0; i < triangleList.size(); i++) {
			if(triangleList.get(i).contains(p)) {
				expr.addTerm(1, vars[i]);
			}
		}
		model.addConstr(expr, GRB.EQUAL, 1, "");
	}

	public void addOverlapConstraints(GRBModel model, GRBVar[] vars, List<Triangle> triangleList, STRtree pointTree) throws GRBException {
		HashMap<Point2D, LinkedList<Integer>> coveringTrianglesMap = new HashMap<Point2D, LinkedList<Integer>>();
		for(int i = 0; i < triangleList.size(); i++) {
			for(Object o : pointTree.query(triangleList.get(i).getEnvelope())) {
				Point2D p = (Point2D) o;
				if(triangleList.get(i).contains(p)) {
					LinkedList<Integer> inMap = coveringTrianglesMap.get(p);
					if(inMap != null) {
						inMap.add(i);
					}else {
						LinkedList<Integer> triangleIndices = new LinkedList<Integer>();
						triangleIndices.add(i);
						coveringTrianglesMap.put(p, triangleIndices);
					}
				}
			}
		}

		for(Point2D p : coveringTrianglesMap.keySet()) {
			GRBLinExpr expr = new GRBLinExpr();
			for(Integer i : coveringTrianglesMap.get(p)) {
				expr.addTerm(1, vars[i]);
			}
			model.addConstr(expr, GRB.EQUAL, 1, "");
		}
	}

	public void addQuantityConstraint(GRBModel model, GRBVar[] vars, List<Point2D> convexHull, Point2D[] points) throws GRBException {
		GRBLinExpr expr = new GRBLinExpr();
		for(int k = 0; k < vars.length; k++) {
			expr.addTerm(1, vars[k]);
		}
		model.addConstr(expr, GRB.EQUAL, points.length*2-(convexHull.size()-1)-2, "k");
	}

	public void addTriangleConstraints(GRBModel model, GRBVar[] vars, ArrayList<Integer[][]> adjTriangles) throws GRBException {
		for(int k = 0; k < adjTriangles.size(); k++) {
			for(Integer[] ints :adjTriangles.get(k)) {
				if(ints.length > 0) {
					GRBLinExpr Expr = new GRBLinExpr();
					for(int h:ints) {
						Expr.addTerm(1, vars[h]);
					}
					String name = "c_" + k;
					model.addConstr(Expr, GRB.GREATER_EQUAL, vars[k], name);
				}
			}
		}
	}

	public void addEdgeSumConstraints(GRBModel model, GRBVar[] vars, HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap) throws GRBException {
		int j = 0;
		for(Map.Entry<HashableLine, ArrayList<LinkedList<Integer>>> e : edgeMap.entrySet()) {
			if(!e.getValue().get(0).isEmpty()) {
				GRBLinExpr leftExpr = new GRBLinExpr();
				for(int k:e.getValue().get(0)) {
					leftExpr.addTerm(1, vars[k]);
				}
				String leftName = "l_" + j + "0";
				model.addConstr(leftExpr, GRB.LESS_EQUAL, 1.0, leftName);
			}
			if(!e.getValue().get(1).isEmpty()) {
				GRBLinExpr rightExpr = new GRBLinExpr();
				for(int k:e.getValue().get(1)) {
					rightExpr.addTerm(1, vars[k]);
				}
				String rightName = "r_" + j + "1";
				model.addConstr(rightExpr, GRB.LESS_EQUAL, 1.0, rightName);
			}
			j++;
		}
	}

	public void addEdgeEqualityConstraints(GRBModel model, GRBVar[] vars, HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap) throws GRBException {
		int j = 0;
		for(Map.Entry<HashableLine, ArrayList<LinkedList<Integer>>> e : edgeMap.entrySet()) {
			if(!e.getValue().get(0).isEmpty() && !e.getValue().get(1).isEmpty()) {
				GRBLinExpr leftExpr = new GRBLinExpr();
				for(int k:e.getValue().get(0)) {
					leftExpr.addTerm(1, vars[k]);
				}
				GRBLinExpr rightExpr = new GRBLinExpr();
				for(int k:e.getValue().get(1)) {
					rightExpr.addTerm(1, vars[k]);
				}
				String name = "eq_" + j;
				model.addConstr(leftExpr, GRB.EQUAL, rightExpr, name);
			}
			j++;
		}
	}

	public void addConvexHullConstraints(GRBModel model, GRBVar[] vars, HashMap<HashableLine, ArrayList<LinkedList<Integer>>> edgeMap, List<Point2D> convexHull) throws GRBException {
		for(int h = 0; h < convexHull.size()-1; h++) {
			GRBLinExpr Expr = new GRBLinExpr();
			for(int f:edgeMap.get(new HashableLine(convexHull.get(h), convexHull.get(h+1))).get(0)) {
				Expr.addTerm(1, vars[f]);
			}
			for(int f:edgeMap.get(new HashableLine(convexHull.get(h), convexHull.get(h+1))).get(1)) {
				Expr.addTerm(1, vars[f]);
			}
			String name = "o_" + h;
			model.addConstr(Expr, GRB.EQUAL, 1.0, name);
		}
	}

	public Point2D getClosestGridPoint(Point2D p, HashMap<Point2D, Double> gridMap) {
		Point2D[] mappedPoints = gridMap.keySet().toArray(new Point2D[gridMap.size()]);

		double minDist = Double.MAX_VALUE;
		int minIdx = 0;
		for(int i = 0; i < mappedPoints.length; i++) {
			double dist = p.distance(mappedPoints[i]);
			if(dist < minDist && !Double.isNaN(gridMap.get(mappedPoints[i]))) {
				minDist = dist;
				minIdx = i;
			}
		}

		return mappedPoints[minIdx];
	}
}
