package processing;

import java.awt.geom.Point2D;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class ElevationPoint {
	
	public Point2D p;
	public double h;
	
	public ElevationPoint(Point2D p, double h) {
		this.p = p;
		this.h = h;
	}
	
	public String toString() {
		return p.toString() + ", " + h;
	}

}
