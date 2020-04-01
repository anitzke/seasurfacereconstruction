package triangulation;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.TreeSet;

import com.vividsolutions.jts.geom.Envelope;

import geometry.ElevationPointClockwiseComparator;
import geometry.HashableLine;
import geometry.Point2DComparator;
import processing.Calculations;
import processing.ElevationPoint;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Triangle implements Cloneable{
  private Point2D a;
  private Point2D b;
  private Point2D c;
  private Envelope env;
  private double hA = 0;
  private double hB = 0;
  private double hC = 0;


  public Triangle(Point2D p, Point2D q, Point2D r) {
    if(Calculations.pointToLine(r, p, q) > 0) {
      a = p;
      b = q;
      c = r;

    } else {
      a = p;
      b = r;
      c = q;
    }

    env = new Envelope();
    env.expandToInclude(a.getX(), a.getY());
    env.expandToInclude(b.getX(), b.getY());
    env.expandToInclude(c.getX(), b.getY());
  }

  public Triangle(Point2D p, double hA, Point2D q, double hB, Point2D r, double hC) {
    if(Calculations.pointToLine(r, p, q) > 0) {
      a = p;
      b = q;
      c = r;
      this.hA = hA;
      this.hB = hB;
      this.hC = hC;
    } else {
      a = p;
      b = r;
      c = q;
      this.hA = hA;
      this.hB = hC;
      this.hC = hB;
    }

    env = new Envelope();
    env.expandToInclude(a.getX(), a.getY());
    env.expandToInclude(b.getX(), b.getY());
    env.expandToInclude(c.getX(), c.getY());

  }

  /**
   * @return the a
   */
  public Point2D getA() {
    return a;
  }

  /**
   * @param a the a to set
   */
  public void setA(Point2D a) {
    this.a = a;
  }

  /**
   * @return the b
   */
  public Point2D getB() {
    return b;
  }

  /**
   * @param b the b to set
   */
  public void setB(Point2D b) {
    this.b = b;
  }

  /**
   * @return the c
   */
  public Point2D getC() {
    return c;
  }

  /**
   * @param c the c to set
   */
  public void setC(Point2D c) {
    this.c = c;
  }

  /**
   * @return the env
   */
  public Envelope getEnv() {
    return env;
  }

  /**
   * @param env the env to set
   */
  public void setEnv(Envelope env) {
    this.env = env;
  }

  /**
   * @return the hA
   */
  public double gethA() {
    return hA;
  }

  /**
   * @param hA the hA to set
   */
  public void sethA(double hA) {
    this.hA = hA;
  }

  /**
   * @return the hB
   */
  public double gethB() {
    return hB;
  }

  /**
   * @param hB the hB to set
   */
  public void sethB(double hB) {
    this.hB = hB;
  }

  /**
   * @return the hC
   */
  public double gethC() {
    return hC;
  }

  /**
   * @param hC the hC to set
   */
  public void sethC(double hC) {
    this.hC = hC;
  }

  /**
   * for each triangle edge check if p is on the right, if so, then p is in the 
   * triangle (counterclockwise ordering of abc)
   * @param p
   * @return
   */
  public boolean contains(Point2D p) {
    if (Calculations.pointToLine(p, a, b) <= 0) return false;
    if (Calculations.pointToLine(p, b, c) <= 0) return false;
    if (Calculations.pointToLine(p, c, a) <= 0) return false;
    return true;
  }
  
  public boolean on_boundary(Point2D p) {
    if (Calculations.pointToLine(p, a, b) != 0) return false;
    if (Calculations.pointToLine(p, b, c) != 0) return false;
    if (Calculations.pointToLine(p, c, a) != 0) return false;
    return true;
  }
  
  public boolean semicontains(Point2D p) {
    if (Calculations.pointToLine(p, a, b) < 0) return false;
    if (Calculations.pointToLine(p, b, c) < 0) return false;
    if (Calculations.pointToLine(p, c, a) < 0) return false;
    return true;
  }
  
  public boolean[] containsOther(Triangle other) {
    
    boolean[] contains = new boolean[3];
    
    if(this.contains(other.getA())) {
      contains[0] = true;
    }
    if(this.contains(other.getB())) {
      contains[1] = true;
    }
    if(this.contains(other.getC())) {
      contains[2] = true;
    }
    return contains;
  }
  
  public boolean[] otherContains(Triangle other) {
    
    boolean[] contains = new boolean[3];
    
    if(other.contains(this.getA())) {
      contains[0] = true;
    }
    if(other.contains(this.getB())) {
      contains[1] = true;
    }
    if(other.contains(this.getC())) {
      contains[2] = true;
    }
    return contains;
  }
  
  public Envelope getEnvelope() {
    return env;
  }
  public boolean intersects(Triangle t2) {

    if(!this.env.intersects(t2.getEnvelope())) {
      return false;
    }

    if (this.contains(t2.getA())) return true;
    if (this.contains(t2.getB())) return true;
    if (this.contains(t2.getC())) return true;

    Line2D AB = new HashableLine(this.a, this.b);
    Line2D BC = new HashableLine(this.b, this.c);
    Line2D CA = new HashableLine(this.c, this.a);

    Line2D ABo = new HashableLine(t2.a, t2.b);
    Line2D BCo = new HashableLine(t2.b, t2.c);
    Line2D CAo = new HashableLine(t2.c, t2.a);

    if(AB.intersectsLine(ABo) && 
        !AB.getP1().equals(ABo.getP1()) && 
        !AB.getP2().equals(ABo.getP2()) && 
        !AB.equals(ABo) &&
        !AB.getP1().equals(ABo.getP2()) && 
        !AB.getP2().equals(ABo.getP1())) {
      return true;
    }
    if(AB.intersectsLine(BCo) && 
        !AB.getP1().equals(BCo.getP1()) && 
        !AB.getP2().equals(BCo.getP2()) && 
        !AB.equals(BCo) &&
        !AB.getP1().equals(BCo.getP2()) && 
        !AB.getP2().equals(BCo.getP1())) {
      return true;
    }
    if(AB.intersectsLine(CAo) && 
        !AB.getP1().equals(CAo.getP1()) && 
        !AB.getP2().equals(CAo.getP2()) && 
        !AB.equals(CAo) &&
        !AB.getP1().equals(CAo.getP2()) && 
        !AB.getP2().equals(CAo.getP1())) {
      return true;
    }
    if(BC.intersectsLine(ABo) && 
        !BC.getP1().equals(ABo.getP1()) && 
        !BC.getP2().equals(ABo.getP2()) && 
        !BC.equals(ABo) &&
        !BC.getP1().equals(ABo.getP2()) && 
        !BC.getP2().equals(ABo.getP1())) {
      return true;
    }
    if(BC.intersectsLine(BCo) && 
        !BC.getP1().equals(BCo.getP1()) && 
        !BC.getP2().equals(BCo.getP2()) && 
        !BC.equals(BCo) &&
        !BC.getP1().equals(BCo.getP2()) && 
        !BC.getP2().equals(BCo.getP1())) {
      return true;
    }
    if(BC.intersectsLine(CAo) && 
        !BC.getP1().equals(CAo.getP1()) && 
        !BC.getP2().equals(CAo.getP2()) && 
        !BC.equals(CAo) &&
        !BC.getP1().equals(CAo.getP2()) && 
        !BC.getP2().equals(CAo.getP1())) {
      return true;
    }
    if(CA.intersectsLine(ABo) && 
        !CA.getP1().equals(ABo.getP1()) && 
        !CA.getP2().equals(ABo.getP2()) && 
        !CA.equals(ABo) &&
        !CA.getP1().equals(ABo.getP2()) && 
        !CA.getP2().equals(ABo.getP1())) {
      return true;
    }
    if(CA.intersectsLine(BCo) && 
        !CA.getP1().equals(BCo.getP1()) && 
        !CA.getP2().equals(BCo.getP2()) && 
        !CA.equals(BCo) &&
        !CA.getP1().equals(BCo.getP2()) && 
        !CA.getP2().equals(BCo.getP1())) {
      return true;
    }
    if(CA.intersectsLine(CAo) && 
        !CA.getP1().equals(CAo.getP1()) && 
        !CA.getP2().equals(CAo.getP2()) && 
        !CA.equals(CAo) &&
        !CA.getP1().equals(CAo.getP2()) && 
        !CA.getP2().equals(CAo.getP1())) {
      return true;
    }

    return false;
  }
  
  public ElevationPoint[] intersection(Triangle other) {
    
    if(this.equals(other)) {
      return new ElevationPoint[] {new ElevationPoint(this.a, this.hA), 
          new ElevationPoint(this.b, this.hB), 
          new ElevationPoint(this.c, this.hC)};
    }
        
    if(this.getEnvelope().intersects(other.getEnvelope())) {
     if(this.intersects(other)) {
       
       // if this contains other completely: intersection = other
       if(this.contains(other.getA()) && this.contains(other.getB()) && this.contains(other.getC())) {
         return new ElevationPoint[] {new ElevationPoint(other.getA(), other.gethA()), 
             new ElevationPoint(other.getB(), other.gethB()), 
             new ElevationPoint(other.getC(), other.gethC())};
       }
       
       // if other contains this completely: intersection = this
       if(other.contains(this.getA()) && other.contains(this.getB()) && other.contains(this.getC())) {
         return new ElevationPoint[] {new ElevationPoint(this.a, this.hA), 
                new ElevationPoint(this.b, this.hB), 
                new ElevationPoint(this.c, this.hC)};
       }
               
       List<ElevationPoint> pointsOfInterest = new ArrayList<>();
       
       Point2D[] pointsThis = new Point2D[] {this.a, this.b, this.c};
       double[] heightsThis = new double[] {this.hA, this.hB, this.hC};
       Point2D[] pointsOther = new Point2D[] {other.getA(), other.getB(), other.getC()};
       double[] heightsOther = new double[] {other.gethA(), other.gethB(), other.gethC()};
       
       for(int i = 0; i < 3; i++) {
         if(this.contains(pointsOther[i])) {
           pointsOfInterest.add(new ElevationPoint(pointsOther[i], heightsOther[i]));
         }
         if(other.contains(pointsThis[i])) {
           pointsOfInterest.add(new ElevationPoint(pointsThis[i], heightsThis[i]));
         }
       }
       
       ElevationPoint[] intersections = getAllIntersections(other);
       pointsOfInterest.addAll(Arrays.asList(intersections));
       
       Point2D[] simplePoints = new Point2D[pointsOfInterest.size()];
       for(int k = 0; k < pointsOfInterest.size(); k++) {
         simplePoints[k] = pointsOfInterest.get(k).p;
       }
       Point2D c = Calculations.getPolyCenter(simplePoints);
       TreeSet<ElevationPoint> clockwisePoly = new TreeSet<>(new ElevationPointClockwiseComparator(c));
       clockwisePoly.addAll(pointsOfInterest);
       
       return clockwisePoly.toArray(new ElevationPoint[clockwisePoly.size()]);       
     }
   }
    
    return new ElevationPoint[0];
  }
     
  public ElevationPoint[] getAllIntersections(Triangle other) {
    List<ElevationPoint> intersections = new ArrayList<>();
    
    HashableLine AB = new HashableLine(this.a, this.b);
    HashableLine BC = new HashableLine(this.b, this.c);
    HashableLine CA = new HashableLine(this.c, this.a);
    
    HashableLine ABo = new HashableLine(other.a, other.b);
    HashableLine BCo = new HashableLine(other.b, other.c);
    HashableLine CAo = new HashableLine(other.c, other.a);
    
    HashableLine[] thisLines = new HashableLine[] {AB, BC, CA};
    double[][] thisHeights = new double[][] {{this.hA, this.hB}, {this.hB, this.hC}, {this.hC, this.hA}};
    HashableLine[] otherLines = new HashableLine[] {ABo, BCo, CAo};
    double[][] otherHeights = new double[][] {{other.hA, other.hB}, {other.hB, other.hC}, {other.hC, other.hA}};

    for(int i = 0; i < thisLines.length; i++) {
      for(int j = 0; j < otherLines.length; j++) {
        HashableLine l1 = thisLines[i];
        HashableLine l2 = otherLines[j];
//        && !l1.getP1().equals(l2.getP1())
//        && !l1.getP1().equals(l2.getP2())
//        && !l1.getP2().equals(l2.getP1())
//        && !l1.getP2().equals(l2.getP2())
        if(!l1.equals(l2) && l1.intersectsLine(l2)) {
          Point2D ip = l1.intersectionPoint(l2);
          
          double segmentLength1 = l1.getP1().distance(l1.getP2());
          double heightDiff1 = thisHeights[i][1] - thisHeights[i][0];
         
          double dist_source_ip1 = l1.getP1().distance(ip);
          double dist_ratio1 = dist_source_ip1/segmentLength1;
         
          double ip_height1 =  thisHeights[i][0] + heightDiff1 * dist_ratio1;
          
          
          double segmentLength2 = l2.getP1().distance(l2.getP2());
          double heightDiff2 = otherHeights[j][1] - otherHeights[j][0];
         
          double dist_source_ip2 = l2.getP1().distance(ip);
          double dist_ratio2 = dist_source_ip2/segmentLength2;
         
          double ip_height2 =  otherHeights[j][0] + heightDiff2 * dist_ratio2;
          
          double ip_height = (ip_height1 + ip_height2) / 2;
          
          ElevationPoint ep = new ElevationPoint(ip, ip_height);
          intersections.add(ep);
        }
      }
    }
    
    return intersections.toArray(new ElevationPoint[intersections.size()]);
  }

  public double area() {
//    proj.project(this);
    double A = Math.sqrt(Math.pow(b.getX()-a.getX(), 2) + Math.pow(b.getY()-a.getY(), 2));
    double B = Math.sqrt(Math.pow(c.getX()-b.getX(), 2) + Math.pow(c.getY()-b.getY(), 2));
    double C = Math.sqrt(Math.pow(a.getX()-c.getX(), 2) + Math.pow(a.getY()-c.getY(), 2));

    double s = (A+B+C)/2.0;

    return Math.sqrt(s * (s - A) * (s - B) * (s - C));
  }
  
  public double areaKM2() {
    double ab = Calculations.distanceGeogrCoord(a, b);
    double bc = Calculations.distanceGeogrCoord(b, c);
    double ac = Calculations.distanceGeogrCoord(a, c);

    double s = (ab + bc + ac)/2.0;
    double area = Math.sqrt(s * (s - ab) * (s - bc) * (s - ac));

    return area;
  }

  public double[] planeParameters(double x, double y) {

    double x21 = this.b.getX() - this.a.getX();
    double x31 = this.c.getX() - this.a.getX();
    double y21 = this.b.getY() - this.a.getY();
    double y31 = this.c.getY() - this.a.getY();

    double frac = y21/x21;

    double s = (y - this.a.getY() - frac * x + frac * this.a.getX()) / (y31 - x31 * frac);

    double r = (x - this.a.getX() - s * x31) / x21;

    return new double[] {r, s};
  }
  /**
   * calculates the hessian normal form of the triangle plane
   * ax + by + cz + d = 0
   * @return [a, b, c, d]
   */

  public double[] planeParameters() {
    double alpha = (b.getY()-a.getY())*(hC-hA)-(c.getY()-a.getY())*(hB-hA);
    double beta = (hB-hA)*(c.getX()-a.getX())-(hC-hA)*(b.getX()-a.getX());
    double gamma = (b.getX()-a.getX())*(c.getY()-a.getY())-(c.getX()-a.getX())*(b.getY()-a.getY());
    double delta = -(alpha*a.getX()+beta*a.getY()+gamma*hA);
    return new double[]{alpha, beta, gamma, delta};
  }

  /**
   * computes the center of the circumcircle of the triangle by intersection of
   * the perpendicular bisectors
   * @return
   */
  public Point2D circumCenter() {
    // mean point of the edges ab and bc e.g. for ab: Mab = a + 0.5*ab
    Point2D Mab = new Point2D.Double(a.getX()+0.5*(b.getX()-a.getX()), 
        a.getY()+0.5*(b.getY()-a.getY()));
    Point2D Mbc = new Point2D.Double(b.getX()+0.5*(c.getX()-b.getX()), 
        b.getY()+0.5*(c.getY()-b.getY()));
    
    // normal vectors
    Point2D Nab = new Point2D.Double(-(b.getY()-a.getY()), (b.getX()-a.getX()));
    Point2D Nbc = new Point2D.Double(-(c.getY()-b.getY()), (c.getX()-b.getX()));

    // normalization factor
    double nab = Math.sqrt(Math.pow(Nab.getX(), 2) + Math.pow(Nab.getY(), 2));
    double nbc = Math.sqrt(Math.pow(Nbc.getX(), 2) + Math.pow(Nbc.getY(), 2));
    
    // normalized normal vector
    Nab = new Point2D.Double(Nab.getX()/nab, Nab.getY()/nab);
    Nbc = new Point2D.Double(Nbc.getX()/nbc, Nbc.getY()/nbc);
            
    // observation vector l_lambda
    double l1 = Mab.getX() - Mbc.getX();
    double l2 = Mab.getY() - Mbc.getY();
    
    // N = AT*A
    double N11 = Nbc.getX()*Nbc.getX() + Nbc.getY()*Nbc.getY();
    double N12 = -Nbc.getX()*Nab.getX() - Nbc.getY()*Nab.getY();
    double N21 = N12;
    double N22 = Nab.getX()*Nab.getX() + Nab.getY()*Nab.getY();
    
    // n = AT*l
    double n1 = Nbc.getX()*l1 + Nbc.getY()*l2;
    double n2 = -Nab.getX()*l1 - Nab.getY()*l2;
    
    // solve x = N⁻¹*n for x1
    double x1 = 1/(N11-N12/N22*N21)*n1 + (-(1/(N11-N12/N22*N21))*N12/N22)*n2;
    
    // center of the circle
    double x_center = Mbc.getX()+x1*Nbc.getX();
    double y_center = Mbc.getY()+x1*Nbc.getY();
    
    Point2D center = new Point2D.Double(x_center, y_center);
    return center;    
  }
  
  /**
   * computes the bounding box/envelope of a circle
   * @param center
   * @param radius
   * @return env: envelope of a circle
   */
  public Envelope envCircle(Point2D center, double radius){
    
    // left and right corner of bounding box over cirumcircle
    Point2D lc = new Point2D.Double(center.getX()-radius, center.getY()-radius);
    Point2D rc = new Point2D.Double(center.getX()+radius, center.getY()+radius);   
    
    Envelope envCircle = new Envelope();
    envCircle.expandToInclude(lc.getX(), lc.getY());
    envCircle.expandToInclude(rc.getX(), rc.getY());
    
    return envCircle;
  }

  @Override
  public String toString() {
    return "Triangle [a=" + a + ", b=" + b + ", c=" + c + ", env=" + env + "]\n";
  }
  @Override
  public boolean equals(Object obj) {
    if(obj == this) {
      return true;
    }
    if(!(obj instanceof Triangle)) {
      return false;
    }
    Triangle t = (Triangle) obj;

    if(!this.getA().equals(t.getA()) && !this.getA().equals(t.getB()) && !this.getA().equals(t.getC())) return false;
    if(!this.getB().equals(t.getA()) && !this.getB().equals(t.getB()) && !this.getB().equals(t.getC())) return false;
    if(!this.getC().equals(t.getA()) && !this.getC().equals(t.getB()) && !this.getC().equals(t.getC())) return false;
    return true;
  }

  @Override
  public int hashCode() {
    TreeSet<Point2D> hashTree = new TreeSet<Point2D>(new Point2DComparator());
    hashTree.add(a);
    hashTree.add(b);
    hashTree.add(c);
    return Objects.hash(hashTree.toArray());
  }
  
  public Object clone()throws CloneNotSupportedException{  
    return (Triangle)super.clone();  
 }

}
