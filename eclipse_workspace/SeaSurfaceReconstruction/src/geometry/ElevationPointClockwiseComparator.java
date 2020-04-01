package geometry;

import java.awt.geom.Point2D;
import java.util.Comparator;

import processing.ElevationPoint;
import processing.ValueComparator;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class ElevationPointClockwiseComparator implements Comparator<ElevationPoint>{

	private Point2D pivot; 
	
	public ElevationPointClockwiseComparator() {
		pivot = new Point2D.Double();
	}
	
	public ElevationPointClockwiseComparator(Point2D p) {
		this.pivot = p;
	}
	
	public ElevationPointClockwiseComparator(ElevationPoint p) {
		this.pivot = p.p;
	}
	
	@Override
	public int compare(ElevationPoint p1, ElevationPoint p2) {
		ValueComparator vc = new ValueComparator();
		double a1 = Math.atan2(p1.p.getY()-pivot.getY(), p1.p.getX()-pivot.getX());
		if(a1 > 0) {
			a1 = Math.PI*2 - a1;
		}else {
			a1 = Math.abs(a1);
		}
		double a2 = Math.atan2(p2.p.getY()-pivot.getY(), p2.p.getX()-pivot.getX());
		if(a2 > 0) {
			a2 = Math.PI*2 - a2;
		}else {
			a2 = Math.abs(a2);
		}
//		if(a1 == a2)return 0;
//		return(a1 > a2) ? 1 : -1;
		
		return vc.compare(a1, a2);
	}
	
}
