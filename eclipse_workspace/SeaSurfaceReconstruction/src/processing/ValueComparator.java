package processing;

import java.util.Comparator;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class ValueComparator implements Comparator<Double>{

  @Override
  public int compare(Double o1, Double o2) {
    double tolerance = Math.pow(10, -9);
    if(o1 - o2 < -tolerance) {
      return -1;
    }else if(o1 - o2 > tolerance) {
      return 1;
    }
    return 0;
  }

}
