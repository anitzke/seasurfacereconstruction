package geometry;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.Objects;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class HashableLine extends Line2D {

  private double x1, x2, y1, y2;


  public HashableLine(){
    super();
  }

  public HashableLine(double x1, double y1, double x2, double y2) {
    this.x1 = x1;
    this.x2 = x2;
    this.y1 = y1;
    this.y2 = y2;

  }

  public HashableLine(Point2D p1, Point2D p2) {
    this.x1 = p1.getX();
    this.x2 = p2.getX();
    this.y1 = p1.getY();
    this.y2 = p2.getY();
  }

  public HashableLine(Line2D l) {
    this.x1 = l.getX1();
    this.x2 = l.getX2();
    this.y1 = l.getY1();
    this.y2 = l.getY2();
  }

  @Override
  public boolean equals(Object obj) {
    if(obj == this) {
      return true;
    }
    if(!(obj instanceof HashableLine)) {
      return false;
    }
    HashableLine l = (HashableLine) obj;


    if(this.getP1().equals(l.getP1()) && this.getP2().equals(l.getP2())) return true;
    if(this.getP2().equals(l.getP1()) && this.getP1().equals(l.getP2())) return true;
    return false;
  }

  @Override
  public Rectangle2D getBounds2D() {
    double x = Math.min(x1, x2);
    double y = Math.min(y1, y2);
    double w = Math.max(x1, x2) - x;
    double h = Math.max(y1, y2) - y;
    return new Rectangle2D.Double(x, y, w, h);
  }

  @Override
  public Point2D getP1() {
    return new Point2D.Double(x1, y1);
  }

  @Override
  public Point2D getP2() {
    return new Point2D.Double(x2, y2);
  }

  @Override
  public double getX1() {
    return x1;
  }

  @Override
  public double getX2() {
    return x2;
  }

  @Override
  public double getY1() {
    return y1;
  }

  @Override
  public double getY2() {
    return y2;
  }

  @Override
  public String toString() {
    return "Line [x1=" + x1 + ", y1=" + y1 + ", x2=" + x2 + ", y2=" + y2 + "]";
  }

  @Override
  public void setLine(double x1, double y1, double x2, double y2) {
    this.x1 = x1;
    this.x2 = x2;
    this.y1 = y1;
    this.y2 = y2;
  }

  @Override
  public int hashCode() {
    if(x1 > x2 || (x1 == x2 && y1 > y2)) return Objects.hash(x2, y2, x1, y1);
    return Objects.hash(x1, y1, x2, y2);
  }

  public void turn() {
    double temp = x1;
    x1 = x2;
    x2 = temp;
    temp = y1;
    y1 = y2;
    y2 = temp;
  }
  
  public Point2D intersectionPoint(Line2D other) {
    
    if(!this.intersectsLine(other) || this.equals(other)) {
      return null;
    }
    
    double denominator = (this.x1 - this.x2) * (other.getY1() - other.getY2()) - 
        (this.y1 - this.y2) * (other.getX1() - other.getX2());
        
    double t = (this.x1 - other.getX1()) * (other.getY1() - other.getY2()) - 
        (this.y1 - other.getY1()) * (other.getX1() - other.getX2());
    t /= denominator;
        
    double x = this.x1 + t*(this.x2 - this.x1);
    double y = this.y1 + t*(this.y2 - this.y1);
    
    return new Point2D.Double(x,y);
  }

}
