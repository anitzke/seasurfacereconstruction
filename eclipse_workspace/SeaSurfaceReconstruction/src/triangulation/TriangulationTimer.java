package triangulation;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class TriangulationTimer {

  private double startFormulating;
  private double endFormulating;
  private double startSolving;
  private double endSolving;

  public TriangulationTimer() {

  }

  public void recordFormulationStart() {
    startFormulating = System.nanoTime()/1000;
  }

  public void recordFormulationEnd() {
    endFormulating = System.nanoTime()/1000;
  }

  public void recordSolveStart() {
    startSolving = System.nanoTime()/1000;
  }

  public void recordSolveEnd() {
    endSolving = System.nanoTime()/1000;
  }

  public double getFormulationTime() {
    return endFormulating - startFormulating;
  }

  public double getSolveTime() {
    return endSolving - startSolving;
  }

  public double[] getTimes() {
    return new double[] {endFormulating - startFormulating, endSolving - startSolving};
  }

}
