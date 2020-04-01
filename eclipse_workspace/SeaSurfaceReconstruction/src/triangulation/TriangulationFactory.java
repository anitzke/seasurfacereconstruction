package triangulation;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class TriangulationFactory {

  private String[] inputParameters;
  private RawData mData;
  private List<TriangulationWorker> workers = new ArrayList<TriangulationWorker>();
  private int nJobs = 1;

  public TriangulationFactory() {

  }

  public TriangulationFactory(String[] parameters) {
    setParameters(parameters);
  }

  public TriangulationFactory(RawData data, int nJobs) {
    this.mData = data;
    this.nJobs = nJobs;
  }

  public void setParameters(String[] parameters) {
    this.inputParameters = parameters;
    setRawData();
  }

  public void setRawData() {
    DataCollector collector = new DataCollector(inputParameters);
    mData = collector.getData();
    nJobs = Integer.parseInt(collector.getParameter("-threads")); 
  }

  /**
   * Starts all TriangulationWorkers
   * @return ResultData
   */
  public ResultData process() {
    ResultData finalOutput = new ResultData();
    Recruiter recruiter = new Recruiter(mData);
    workers = recruiter.recruitWorkers(nJobs);

    List<ResultData> partialOutputs = new ArrayList<ResultData>();
    Thread[] threads = new Thread[workers.size()];
    for(int i = 0; i < threads.length; i++) {
      final Integer fi = i;
      threads[i] = new Thread(new Runnable() {
        @Override
        public void run() {
          // run worker on data
          partialOutputs.add(workers.get(fi).work());
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

    // compile the results
    for(ResultData out : partialOutputs) {
      finalOutput.addNames(out.getNames());
      finalOutput.addAnomalies(out.getAnomalies());
      finalOutput.addSquaredAnomalies(out.getSquaredAnomalies());
      finalOutput.addStatistics(out.getStatistics());
      finalOutput.addTriangulations(out.getTriangulations());
      finalOutput.addValidations(out.getValidations());
      finalOutput.addTimes(out.getTimes());
    }

    return finalOutput;
  }


}
