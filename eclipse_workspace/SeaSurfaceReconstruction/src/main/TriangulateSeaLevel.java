package main;

import java.sql.Timestamp;

import triangulation.DataCollector;
import triangulation.ResultData;
import triangulation.TriangulationFactory;

/**
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class TriangulateSeaLevel {

  public static void main(String[] args) {
    System.getProperty("user.dir");
    
    // start time
    Timestamp timestamp = new Timestamp(System.currentTimeMillis());
    System.out.println(timestamp);

    // read altimeter data and tide gauges from files
    DataCollector collector = new DataCollector(args);

    // create workers
    TriangulationFactory factory = new TriangulationFactory(collector.getData(), Integer.parseInt(DataCollector.getParameter("-threads")));
    // run workers on data
    ResultData results = factory.process();

    // store results to files
//    DataDistributor distributor = new DataDistributor(results, DataCollector.getParameter("-store"));
//    System.out.println("write all results within main");
//    distributor.distribute();
    
    // end time
    timestamp = new Timestamp(System.currentTimeMillis());
    System.out.println(timestamp);
  }
}


