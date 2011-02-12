/* 
 * Foliage. An Ancestral Recombination Graph Manipulation Library.
 * 
 * Copyright (c) 2008 Genome Research Ltd.
 * 
 * Author: Lior Galanti <lior.galanti@gmail.com>
 * 
 * This file is part of Foliage.
 * Foliage is free software; you can redistribute it and/or modify it under the terms of 
 * the GNU General Public License as published by the Free Software Foundation; 
 * either version 2 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 * 
 */

package sanger.argml.graph.model;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import javax.xml.stream.XMLStreamException;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.argml.format.ArgmlDocument;
import sanger.argml.io.XmlOutput;
import sanger.argml.statistic.Calculator;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetCollection;
import sanger.math.set.NaturalSetException;

public class Statistics extends Environmental{
	protected String name;
	protected NaturalDomain snpDomain;
	protected NaturalDomain haplotypeDomain;
	protected NaturalDomain basePairDomain;
	
	protected int[] markerPositions;	
	protected double[] recombinationCount;
	protected double[][] treeCorrelationCount;
	
	protected int args;
	//private double treeCorrelationFactor;
	private double recombinationRateFactor;
	private double recombinationMax;
	private double minvalue;
	private double maxvalue;
	private double range;
	
	
	
	// Constractors
	
	protected Statistics(Environment env, String name){
		super(env);
		this.name = name;
		snpDomain = null;
		haplotypeDomain = null;
		basePairDomain = null;
		
		markerPositions = null;	
		recombinationCount = null;
		args = 0;
		recombinationMax = 0;
		treeCorrelationCount = null;
	}
		
	public Statistics(Environment env, String name, HaplotypeReader haplotypeInput) throws IOException{
		super(env);
		this.name = name;
		snpDomain = haplotypeInput.getSnpDomain();
		haplotypeDomain = haplotypeInput.getHaplotypeDomain();
		basePairDomain = haplotypeInput.getBasePairDomain();
		markerPositions = haplotypeInput.getMarkerPositions();
		haplotypeInput.close();
		
		recombinationCount = new double[snpDomain.closureCardinality()];
		args = 0;
		recombinationMax = 0;
		treeCorrelationCount = new double[snpDomain.closureCardinality()][];		
		for(int i=0; i<snpDomain.closureCardinality(); i++) { treeCorrelationCount[i] = new double[i+1]; }
	}
		

	public void addStatistics(Statistics other) throws NaturalSetException{
		if(snpDomain.equals(other.snpDomain) && haplotypeDomain.equals(other.haplotypeDomain)){
			for(int i=0; i<recombinationCount.length; i++){
				recombinationCount[i] += other.recombinationCount[i];
			}
			
			for(int i=0; i<treeCorrelationCount.length; i++){
				for(int j=0; j<treeCorrelationCount[i].length; j++){
					treeCorrelationCount[i][j] += other.treeCorrelationCount[i][j];
				}
			}
			
			args += other.args;			
			updateFactors();
			
		} else {
			throw new NaturalSetException("Other Statistics domain incompatible");
		}
	}
	
	public void addGenealogy(Genealogy genealogy) throws NaturalSetException {
		if ( snpDomain.equals(genealogy.snpDomain()) && haplotypeDomain.equals(genealogy.haplotypeDomain())){
			updateRecombinationCount(genealogy);
			updateLocalTreesCorrelation(genealogy);
			args++;
			updateFactors();
		} else {
			throw new NaturalSetException("Statistics and Genealogy must be of the same base-pair and haplotype domains");
		}
	}
	
		
	public static Statistics clip(Statistics instance, NaturalSet region) throws NaturalSetException{
		Statistics fragment = new Statistics(instance.env(), instance.name);
		if(instance.snpDomain.equals(region.domain())){
			fragment.args = instance.args;
			fragment.snpDomain = instance.snpDomain.createSubDomain(instance.snpDomain.project(region));
			NaturalSet bpRegion = new NaturalSet(instance.basePairDomain, instance.markerPositions[instance.snpDomain.toRelativeCoordinate(region.min())], instance.markerPositions[instance.snpDomain.toRelativeCoordinate(region.max())]);
			fragment.basePairDomain = instance.basePairDomain.createSubDomain(instance.basePairDomain.project(bpRegion));
			fragment.haplotypeDomain = instance.haplotypeDomain;
	
			fragment.markerPositions = new int[fragment.snpDomain.closureCardinality()];		
			fragment.recombinationCount = new double[fragment.snpDomain.closureCardinality()];
			fragment.treeCorrelationCount = new double[fragment.snpDomain.closureCardinality()][];
	
			if(fragment.snpDomain.isContinuous()){
				for(int i=fragment.snpDomain.min(); i<=fragment.snpDomain.max(); i++){
					fragment.markerPositions[fragment.snpDomain.toRelativeCoordinate(i)] = instance.coordinate(i);
					fragment.recombinationCount[fragment.snpDomain.toRelativeCoordinate(i)] = instance.recombinationCount(i);
				}
				
			} else {
				for(int i=fragment.snpDomain.min(); i<=fragment.snpDomain.max(); i++){
					if(fragment.snpDomain.contains(i)) {
						fragment.markerPositions[fragment.snpDomain.toRelativeCoordinate(i)] = instance.coordinate(i);
						fragment.recombinationCount[fragment.snpDomain.toRelativeCoordinate(i)] = instance.recombinationCount(i);
	
					} else {
						fragment.markerPositions[fragment.snpDomain.toRelativeCoordinate(i)] = 0;
						fragment.recombinationCount[fragment.snpDomain.toRelativeCoordinate(i)] = 0;
					}
				}
			}
			
	
			for(int k=0; k<fragment.snpDomain.closureCardinality(); k++){
				fragment.treeCorrelationCount[k] = new double[k+1];			
				int end = fragment.snpDomain.toAbsoluteCoordinate(k);
				if(fragment.snpDomain.isContinuous()){
					for(int i=fragment.snpDomain.min(); i<=end; i++){
						fragment.treeCorrelationCount[k][fragment.snpDomain.toRelativeCoordinate(i)] = instance.localTreeCorrelationCount(end,i);
					}
					
				} else {
					for(int i=fragment.snpDomain.min(); i<=end; i++){
						if(fragment.snpDomain.contains(i)) {
							fragment.treeCorrelationCount[k][fragment.snpDomain.toRelativeCoordinate(i)] = instance.localTreeCorrelationCount(end,i);
						} else {
							fragment.treeCorrelationCount[k][fragment.snpDomain.toRelativeCoordinate(i)] = 0;
						}
					}
				}			
			}

			fragment.updateFactors();	

		} else {
			throw new NaturalSetException("Statistics should be clipped on its on snp domain.");
		}
		return fragment;		
	}

	public double recombinationMax() throws NaturalSetException{
		return (double)recombinationMax * recombinationRateFactor;
	}

	
	// Accessors

	public NaturalDomain haplotypeDomain() {
		return haplotypeDomain;
	}
	
	public NaturalDomain snpDomain() {
		return snpDomain;
	}
	
	public NaturalDomain basePairDomain() {
		return basePairDomain;
	}
	
	public int numberOfGenealogies() {
		return args;
	}

	public int coordinate(int i) throws NaturalSetException{
		return markerPositions[snpDomain.toRelativeCoordinate(i)];
	}

	private double localTreeCorrelationCount(int i, int j) throws NaturalSetException{
		int si = Math.max(i, j);
		int sj = Math.min(i, j);
		return (double)treeCorrelationCount[snpDomain.toRelativeCoordinate(si)][snpDomain.toRelativeCoordinate(sj)];
	}

	public double localTreeCorrelation(int i, int j) throws NaturalSetException{
		int si = Math.max(i, j);
		int sj = Math.min(i, j);
		return normalizeDistance((double)treeCorrelationCount[snpDomain.toRelativeCoordinate(si)][snpDomain.toRelativeCoordinate(sj)]);
	}
	
	private double recombinationCount(int i) throws NaturalSetException{
		return (double)recombinationCount[snpDomain.toRelativeCoordinate(i)];
	}

	public double recombinationRate(int i) throws NaturalSetException{
		return recombinationCount(i) * recombinationRateFactor;
	}
	
	public double[] recombinationRates() {
		double[] result = new double[recombinationCount.length];
		for(int i=0; i<result.length; i++){
			result[i] = recombinationCount[i] * recombinationRateFactor;
		}
		return result;
	}

	// XML Writer Procedures
	public void writeDocument(XmlOutput out) throws XMLStreamException, FileNotFoundException, IOException{
		ArgmlDocument doc = new ArgmlDocument(out);
		doc.writeStartDocument();
		writeElement(out);
		doc.writeEndDocument();
	}
	
	public void writeElement(XmlOutput out) throws XMLStreamException, FileNotFoundException, IOException{
		out.writer().writeStartElement("", "statistics");
		out.writer().writeAttribute("genealogies", String.valueOf(args));
		out.writer().writeAttribute("name", name);
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);

		haplotypeDomain.writeXml(out, "haplotype");
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		snpDomain.writeXml(out, "snp");
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		basePairDomain.writeXml(out, "basepair");
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		
		StringBuilder sb = new StringBuilder();

		out.writer().writeStartElement("", "marker");
		for(int m : markerPositions){
			sb.append(' ');
			sb.append(m);
		} sb.delete(0, 1);
		out.writer().writeCharacters(sb.toString());
		out.writer().writeEndElement();
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		
		out.writer().writeStartElement("", "recombination");
		sb = new StringBuilder();
		for(double r : recombinationCount){
			sb.append(' ');
			sb.append(r);
		} sb.delete(0, 1);
		out.writer().writeCharacters(sb.toString());
		out.writer().writeEndElement();
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		
		out.writer().writeStartElement("", "localtreecorrelation");
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		for(double[] row : treeCorrelationCount){
			sb = new StringBuilder();
			for(double cell : row){
				sb.append(' ');
				sb.append(cell);
			} sb.delete(0, 1);
			out.writer().writeStartElement("row");
			out.writer().writeCharacters(sb.toString());
			out.writer().writeEndElement();
			out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		}
		out.writer().writeEndElement();
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		out.writer().writeEndElement();
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
	}
		
	// Service function	
	private void updateRecombinationCount(Genealogy genealogy) throws NaturalSetException {
		int[] distribution = genealogy.recombinationRates();
		for(int i=0; i<recombinationCount.length; i++){
			recombinationCount[i] += distribution[i];
		}
	}
	
	private void updateLocalTreesCorrelation(Genealogy genealogy) throws NaturalSetException {
		HashMap<NaturalSet, NaturalSetCollection>  localTreesBiPartitions = genealogy.localTreesBiPartitions();
		ArrayList<NaturalSet> frames = new ArrayList<NaturalSet>(localTreesBiPartitions.keySet());
		Collections.sort(frames);
		
		for(NaturalSet x : frames){
			for(NaturalSet y : frames){
				if(x.compareTo(y) >= 0){
					NaturalSetCollection xbp = localTreesBiPartitions.get(x);
					NaturalSetCollection ybp = localTreesBiPartitions.get(y);
					int value = xbp.intersectCount(ybp);
					
					if(env().stringProperty("DistanceMetric").equals("bs")){
						value = xbp.size() + ybp.size() - haplotypeDomain().cardinality() - value;
					}
					
					int xmin = x.min(), xmax = x.max(), ymin = y.min(), ymax = y.max();
					for(int i=xmin; i<=xmax; i++){
						for(int j=ymin; j<=ymax && j<=i; j++){
							if(snpDomain.contains(i) && snpDomain.contains(j)){
								treeCorrelationCount[snpDomain.toRelativeCoordinate(i)][snpDomain.toRelativeCoordinate(j)] += value;
							}
						}
					}					
				}
			}
		}
	}

	protected void updateFactors(){
		if(env().stringProperty("DistanceMetric").equals("bs")){			
			maxvalue = Integer.MIN_VALUE;
			minvalue = Integer.MAX_VALUE;
			for(int i=0; i<treeCorrelationCount.length; i++){
				for(int j=0; j<treeCorrelationCount[i].length; j++){
					maxvalue = Math.max(treeCorrelationCount[i][j], maxvalue);
					minvalue = Math.min(treeCorrelationCount[i][j], minvalue);
				}
			}
		}
		
		if(env().stringProperty("DistanceMetric").equals("sbp")){			
			maxvalue = (double)args * (haplotypeDomain.cardinality() - 3);
			minvalue = 0;
		}
		
		range = maxvalue - minvalue;
		recombinationRateFactor = 1.0 / (double)args;		
		recombinationMax = Calculator.max(recombinationCount);
	}

	protected double normalizeDistance(double value){
		double v = (value - minvalue) / range;
		if(env().stringProperty("DistanceMetric").equals("bs")){ v = 1 - v; }
		return v;
	}
	
	public String toString() {
		StringBuilder display = new StringBuilder();
		display.append(name);
		display.append(" {");
		display.append(" SNPS: ");
		display.append(snpDomain );
		display.append(", SEQS: ");
		display.append(haplotypeDomain );					
		display.append(", ARGS: ");
		display.append(args);					
		display.append(", MAXRECS: ");
		display.append(recombinationMax * recombinationRateFactor);
		display.append(" }");					
		return display.toString();
	}

	
	/**
	 * @return the name
	 */
	public String name() {
		return name;
	}
	
	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	
	public int[] markerPositions() {
		return markerPositions;
	}

}
