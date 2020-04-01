package triangulation;

import java.util.ArrayList;
import java.util.List;

/**
 * container-class for output data
 * 
 * @author Alina Förster
 *
 * @date Apr 1, 2020
 */
public class ResultData {

  private List<String> names = new ArrayList<String>();
  private List<double[][]> anomalies = new ArrayList<double[][]>();
  private List<double[][]> squaredAnomalies = new ArrayList<double[][]>();
  private List<String> statistics = new ArrayList<String>();
  private List<TriangulationInstance> triangulations = new ArrayList<TriangulationInstance>();
  private List<String> validations = new ArrayList<String>();
  private List<String> times = new ArrayList<String>();
  private List<double[][]> interpolatedValues = new ArrayList<double[][]>();
  private List<String> interpolationMeans = new ArrayList<String>();
  private List<double[][]> sqaDiffList = new ArrayList<double[][]>();


  public ResultData() {
  }

  public List<String> getNames() {
    return names;
  }

  public List<double[][]> getAnomalies() {
    return anomalies;
  }

  public List<double[][]> getSquaredAnomalies() {
    return squaredAnomalies;
  }
  
  public List<double[][]> getInterpolatedValues() {
    return interpolatedValues;
  }

  public List<String> getStatistics() {
    return statistics;
  }

  public List<TriangulationInstance> getTriangulations() {
    return triangulations;
  }

  public List<String> getValidations() {
    return validations;
  }

  public List<String> getTimes() {
    return times;
  }
  
  public List<String> getInterpolationMeans() {
    return interpolationMeans;
  }

  public List<double[][]> getSqaDiffs() {
    return sqaDiffList;
  }

  public void setSqaDiffs(List<double[][]> sqaDiffList) {
    this.sqaDiffList = sqaDiffList;
  }

  public void setNames(List<String> names) {
    this.names = names;
  }

  public void setAnomalies(List<double[][]> anomalies) {
    this.anomalies = anomalies;
  }

  public void setSquaredAnomalies(List<double[][]> squaredAnomalies) {
    this.squaredAnomalies = squaredAnomalies;
  }
  
  public void setInterpolatedValues(List<double[][]> interpolatedValues) {
    this.interpolatedValues = interpolatedValues;
  }

  public void setStatistics(List<String> statistics) {
    this.statistics = statistics;
  }

  public void setTriangulations(List<TriangulationInstance> triangulations) {
    this.triangulations = triangulations;
  }

  public void setValidations(List<String> validations) {
    this.validations = validations;
  }

  public void setTimes(List<String> times) {
    this.times = times;
  }
  
  public void setInterpolationMeans(List<String> interpolationMeans) {
    this.interpolationMeans = interpolationMeans;
  }

  public void addNames(List<String> names) {
    this.names.addAll(names);
  }

  public void addAnomalies(List<double[][]> anomalies) {
    this.anomalies.addAll(anomalies);
  }

  public void addSquaredAnomalies(List<double[][]> sqAnomalies) {
    this.squaredAnomalies.addAll(sqAnomalies);
  }
  
  public void addInterpolatedValues(List<double[][]> interpolations) {
    this.interpolatedValues.addAll(interpolations);
  }

  public void addStatistics(List<String> statistics) {
    this.statistics.addAll(statistics);
  }

  public void addTriangulations(List<TriangulationInstance> triangulations) {
    this.triangulations.addAll(triangulations);
  }

  public void addValidations(List<String> validations) {
    this.validations.addAll(validations);
  }

  public void addTimes(List<String> times) {
    this.times.addAll(times);
  }
  
  public void addInterpolationMeans(List<String> iMeans) {
    this.interpolationMeans.addAll(iMeans);
  }
  
  public void addSqaDiffs(List<double[][]> sqa) {
    this.sqaDiffList.addAll(sqa);
  }
}
