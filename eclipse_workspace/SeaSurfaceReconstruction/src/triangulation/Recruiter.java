package triangulation;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class Recruiter {

  private RawData mData;

  public Recruiter(RawData data) {
    this.mData = data;
  }

  /**
   * Creates n TriangulationWorkers
   * @param n
   * @return List of TriangulationWorkers
   */
  public List<TriangulationWorker> recruitWorkers(int n) {

    if(mData instanceof RawData.Evaluate) {
      switch(DataCollector.getParameter("-mode")) {
      case "reconstruct":
        return recruitReconstructors(n);
      case "evaluate":
        return recruitEvaluators(n);
      default:
        break;  
      }
    }else if(mData instanceof RawData.Triangulate) {
      return recruitTriangulators(n);
    }else if(mData instanceof RawData.Validate) {
      return recruitValidators(n);
    }else {
      System.out.println("what!? how!?");
    }	
    return null;
  }

  /**
   * creates n TriangulationWorker.Triangulate
   * @param n
   * @return list of TriangulationWorkers
   */
  private List<TriangulationWorker> recruitTriangulators(int n) {
    List<TriangulationWorker> workers = new ArrayList<TriangulationWorker>();
    RawData.Triangulate tempData = (RawData.Triangulate) mData;

    int tasks = (int) ((double)tempData.identifiers.size() / (double)n + 0.5);

    for(int i = 0; i < n; i++) {
      int startIdx = i * tasks;
      int endIdx = ((i+1) * tasks < tempData.identifiers.size()) ? (i+1) * tasks : tempData.identifiers.size();

      TriangulationWorker.Triangulate worker = new TriangulationWorker.Triangulate();
      worker.identifiers = tempData.identifiers.subList(startIdx, endIdx);
      worker.points = tempData.points.subList(startIdx, endIdx);
      if(tempData.rasters.size() > 1) {
        worker.rasters = tempData.rasters.subList(startIdx, endIdx);
      }else {
        worker.rasters = tempData.rasters;
      }
      worker.pointBounds = tempData.pointBounds;
      worker.rasterBounds = tempData.rasterBounds;
      worker.tc = tempData.tc;

      workers.add(worker);
    }

    return workers;
  }

  /**
   * creates n TriangulationWorker.Evaluate
   * @param n
   * @return list of TriangulationWorkers
   */
  private List<TriangulationWorker> recruitEvaluators(int n) {
    List<TriangulationWorker> workers = new ArrayList<TriangulationWorker>();
    RawData.Evaluate tempData = (RawData.Evaluate) mData;
    
    if(n != 1) {
      int num_res = tempData.epochI.size()*tempData.epochJ.size();
      int tasks = (int) ((double)num_res / (double)n + 0.5);

      for(int i = 0; i < n; i++) {
        int startIdx = i * tasks;
        int endIdx = ((i+1) * tasks < num_res) ? (i+1) * tasks : num_res;

        TriangulationWorker.Evaluate worker = new TriangulationWorker.Evaluate();
        worker.identifiers = tempData.identifiers;
        worker.epochI = tempData.epochI;
        worker.epochJ = tempData.epochJ;
        worker.points = tempData.points.subList(startIdx, endIdx);
        if(tempData.rasters.size() > 1) {
          worker.rasters = tempData.rasters.subList(startIdx, endIdx);
        }else {
          worker.rasters = tempData.rasters;
        }
        worker.pointBounds = tempData.pointBounds;
        worker.rasterBounds = tempData.rasterBounds;
        worker.tc = tempData.tc;
        worker.proj = tempData.proj;
        worker.setSingleEvaluationTriangulation(tempData.singleEvaluation);

        workers.add(worker);
      } 
    } else {
      TriangulationWorker.Evaluate worker = new TriangulationWorker.Evaluate();
      worker.identifiers = tempData.identifiers;
      worker.epochI = tempData.epochI;
      worker.epochJ = tempData.epochJ;
      worker.points = tempData.points;
      worker.rasters = tempData.rasters;
      worker.pointBounds = tempData.pointBounds;
      worker.rasterBounds = tempData.rasterBounds;
      worker.tc = tempData.tc;
      worker.proj = tempData.proj;
      worker.setSingleEvaluationTriangulation(tempData.singleEvaluation);

      workers.add(worker);
    }

    return workers;
  }

  /**
   * creates n TriangulationWorker.Reconstruct
   * @param n
   * @return list of TriangulationWorkers
   */
  private List<TriangulationWorker> recruitReconstructors(int n) {
    List<TriangulationWorker> workers = new ArrayList<TriangulationWorker>();
    RawData.Evaluate tempData = (RawData.Evaluate) mData;

    if(n != 1) {
      int num_res = tempData.epochI.size()*tempData.epochJ.size();
      int tasks = (int) ((double)num_res / (double)n + 0.5);

      for(int i = 0; i < n; i++) {
        int startIdx = i * tasks;
        int endIdx = ((i+1) * tasks < num_res) ? (i+1) * tasks : num_res;

        TriangulationWorker.Reconstruct worker = new TriangulationWorker.Reconstruct();
        worker.identifiers = tempData.identifiers;
        worker.epochI = tempData.epochI;
        worker.epochJ = tempData.epochJ;
        worker.points = tempData.points.subList(startIdx, endIdx);
        if(tempData.rasters.size() > 1) {
          worker.rasters = tempData.rasters.subList(startIdx, endIdx);
        }else {
          worker.rasters = tempData.rasters;
        }
        worker.pointBounds = tempData.pointBounds;
        worker.rasterBounds = tempData.rasterBounds;
        worker.tc = tempData.tc;
        worker.proj = tempData.proj;

        workers.add(worker);
    } 
    } else {
      TriangulationWorker.Reconstruct worker = new TriangulationWorker.Reconstruct();
      worker.identifiers = tempData.identifiers;
      worker.epochI = tempData.epochI;
      worker.epochJ = tempData.epochJ;
      worker.points = tempData.points;
      worker.rasters = tempData.rasters;
      worker.pointBounds = tempData.pointBounds;
      worker.rasterBounds = tempData.rasterBounds;
      worker.tc = tempData.tc;
      worker.proj = tempData.proj;

      workers.add(worker);
    }

    return workers;
  }
  
  /**
   * creates n TriangulationWorker.Validate
   * @param n
   * @return list of TriangulationWorkers
   */
  private List<TriangulationWorker> recruitValidators(int n) {
    List<TriangulationWorker> workers = new ArrayList<TriangulationWorker>();
    RawData.Validate tempData = (RawData.Validate) mData;

    int tasks = (int) ((double)tempData.identifiers.size() / (double)n + 0.5);

    for(int i = 0; i < n; i++) {
      int startIdx = i * tasks;
      int endIdx = ((i+1) * tasks < tempData.identifiers.size()) ? (i+1) * tasks : tempData.identifiers.size();

      TriangulationWorker.Validate worker = new TriangulationWorker.Validate();
      worker.identifiers = tempData.identifiers.subList(startIdx, endIdx);
      worker.triangulations = tempData.triangulations.subList(startIdx, endIdx);

      workers.add(worker);
    }

    return workers;
  }

}
