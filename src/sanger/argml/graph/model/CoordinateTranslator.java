package sanger.argml.graph.model;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSetException;

public class CoordinateTranslator extends Environmental{

	protected NaturalDomain source;
	protected NaturalDomain target;
	protected int[] embedding;
	
	public CoordinateTranslator(Environment env, NaturalDomain source, NaturalDomain target, int[] embedding){
		super(env);
		this.source = source;
		this.target = target;
		this.embedding = embedding;
	}
	
	protected double coverage(int i) throws NaturalSetException{
		return coordinate(i+1) - coordinate(i);
	}
	
	protected double overlappingCoverage(int i, double min, double max) throws NaturalSetException{
		return	((min < coordinate(i+1)) && (max > coordinate(i))) ? (Math.min(coordinate(i+1), max) - Math.max(coordinate(i), min)) : 0.0;
	}
	
	protected double overlappingCoverageFraction(int i, double min, double max) throws NaturalSetException{
		return overlappingCoverage(i, min, max) / coverage(i);
	}
	
	public double[] tiledAverage(double tile, double[] parameter) throws NaturalSetException {
		double[] average = new double[(int)Math.ceil(target.closureCardinality() / (double)tile)];
		int left = source.min();
		int count = 0;
		double mark = target.min();

		for(int i : source){
			while(coordinate(i) > mark + tile){
				for(int k=left; k<i; k++){
					average[count] +=  parameter[source.toRelativeCoordinate(k)] / (i-left);
				}
				count++;
				mark += tile;
				left = i-1;
			}
		}
		
		// Complete the last, incomplete, tile
		if(mark < target.max()){
			for(int k=left; k<source.max(); k++){
				average[count] +=  parameter[source.toRelativeCoordinate(k)] / (source.max()-left);
			}
		}

		return  average;
	}
		
	public double[] tiledDensity(double tile, double[] parameter) throws NaturalSetException {
		double[] density = new double[(int)Math.ceil(target.closureCardinality() / (double)tile)];
		int left = source.min();
		int count = 0;
		double mark = target.min();

		for(int i : source){
			while(coordinate(i) > mark + tile){
				for(int k=left; k<i; k++){
					
					density[count] += overlappingCoverageFraction(k, mark, mark + tile) * parameter[source.toRelativeCoordinate(k)];
				}
				count++;
				mark += tile;
				left = i-1;
			}
		}
		
		// Complete the last, incomplete, tile
		if(mark < target.max()){
			for(int k=left; k<source.max(); k++){
				density[count] += overlappingCoverageFraction(k, mark, target.max()) * parameter[source.toRelativeCoordinate(k)];
			}
		}

		return  density;
	}

	public double[] sourceInTargetDensity(double tile) throws NaturalSetException {		
		double[] density = new double[(int)Math.ceil(target.closureCardinality() / (double)tile)];
		int left = source.min();
		int count = 0;
		double mark = target.min();

		for(int i : source){
			while(coordinate(i) > mark + tile){
				for(int k=left; k<i; k++){
					density[count] += overlappingCoverageFraction(k, mark, mark + tile);
				}
				count++;
				mark += tile;
				left = i-1;
			}
		}
		
		// Complete the last, incomplete, tile
		if(mark < target.max()){
			for(int k=left; k<source.max(); k++){
				density[count] += overlappingCoverageFraction(k, mark, target.max());
			}
		}

		return  density;
	}

	protected int coordinate(int i) throws NaturalSetException {
		return embedding[source.toRelativeCoordinate(i)];
	}
}
