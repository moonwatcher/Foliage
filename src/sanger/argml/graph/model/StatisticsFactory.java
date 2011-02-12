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

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import sanger.argml.environment.Environmental;
import sanger.argml.io.XmlInput;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetException;

public class StatisticsFactory extends Environmental implements Iterable<Statistics> {
	private XmlInput input;
	private XMLStreamReader parser;
	private boolean more;
	private String nextName;
	
	private Pattern nameFilter = null;
	private Integer lowerSnpFilter = null;
	private Integer upperSnpFilter = null;
	
	
	public StatisticsFactory(XmlInput input, String filter) throws XMLStreamException, NaturalSetException{
		super(input.env());
		if(filter!=null) nameFilter = Pattern.compile(filter);
		this.input = input;
		this.parser = input.parser();
		this.more = true;
		boolean stop = false;
		while(!stop){
			switch (parser.next()) {
				case XMLStreamConstants.START_ELEMENT:
					if (parser.getLocalName().equals("statistics")) {
						nextName = parser.getAttributeValue("", "name");
						if(nameFilter == null || nameFilter.matcher(nextName).matches()){
							stop = true;
						}
					}
				break;
				
				case XMLStreamConstants.END_DOCUMENT:
					stop = true;
					more = false;
				break;
			}
		}			
	}
	
	private Statistics readNext() throws XMLStreamException, NaturalSetException{
		Statistics s = null;
		
		s = new Statistics(input.env(), parser.getAttributeValue("", "name"));
		s.args = Integer.parseInt(parser.getAttributeValue("", "genealogies"));
		int index = 0;
		boolean stop = false;
		boolean skip = false;
						
		while(!stop) {
			switch (parser.next()) {
				case XMLStreamConstants.START_ELEMENT:
					if(!skip && parser.getLocalName().equals("marker")) {
						s.markerPositions = new int[s.snpDomain().closureCardinality()];
						String[] values = parser.getElementText().split("[\\s]+");
						for(int i=0; i <values.length; i++){
							s.markerPositions[i] = Integer.parseInt(values[i]);
						}
						
					} else if(!skip && parser.getLocalName().equals("domain")) {
							String d = parser.getAttributeValue("", "name");
							if(d.equals("haplotype")) s.haplotypeDomain = NaturalDomain.parser(parser);
							if(d.equals("snp")) s.snpDomain = NaturalDomain.parser(parser);
							if(d.equals("basepair")) s.basePairDomain = NaturalDomain.parser(parser);
							
					} else if (!skip && parser.getLocalName().equals("recombination")) {
						s.recombinationCount = new double[s.snpDomain().closureCardinality()];
						String[] values = parser.getElementText().split("[\\s]+");
						for(int i=0; i <values.length; i++){
							s.recombinationCount[i] = Double.parseDouble(values[i]);
						}
						
					} else if (!skip && parser.getLocalName().equals("localtreecorrelation")) {
						s.treeCorrelationCount = new double[s.snpDomain.closureCardinality()][];
						for(int i=0; i<s.snpDomain.closureCardinality(); i++) { s.treeCorrelationCount[i] = new double[i+1]; }
						
					} else if (!skip && parser.getLocalName().equals("row")) {
						String[] values = parser.getElementText().split("[\\s]+");
						for(int i=0; i <values.length; i++){
							s.treeCorrelationCount[index][i] = Double.parseDouble(values[i]);
						}
						index++;
					} else if(parser.getLocalName().equals("statistics")) {
						nextName = parser.getAttributeValue("", "name");
						if(nameFilter == null || nameFilter.matcher(nextName).matches()){
							stop = true;
						} else {
							skip = true;
						}
					}
				break;
				
				case XMLStreamConstants.END_DOCUMENT:
					stop = true;
					more = false;
				break;				
			}
		}
		s.updateFactors();
		return s;
	}
	
	public Iterator<Statistics> iterator(){
		return new StatisticsIterator(this);
	}
	
	private class StatisticsIterator implements Iterator<Statistics> {
		StatisticsFactory factory;

		public StatisticsIterator(StatisticsFactory factory) {
			this.factory = factory;
		}
		
		public boolean hasNext() {
			return (factory.more);
		}

		public Statistics next() {
			Statistics s = null;
			if(hasNext()){
				try { 
					s = factory.readNext();
					s = clip(s);
					
				} catch (Exception e) {
					factory.env().log().printError(e);
					factory.more = false;
				}
			} else {
				throw new NoSuchElementException("No more statistics elements.");
			}
			return s;
		}

		public void remove() {}
	}

	public String toString() {
		StringBuilder display = new StringBuilder();
		display.append("Statistics Factory{");				
		display.append(" }");					
		return display.toString();
	}	
	
	private Statistics clip(Statistics s) throws NaturalSetException{
		if(lowerSnpFilter!=null || upperSnpFilter!=null){
			NaturalSet set = new NaturalSet(
				s.snpDomain(), 
				lowerSnpFilter!=null ? lowerSnpFilter : s.snpDomain().min(), 
				upperSnpFilter!=null ? upperSnpFilter : s.snpDomain().max());
			s = Statistics.clip(s, set);			
		}
		return s;
	}

	/**
	 * @return the lowerSnpFilter
	 */
	public Integer getLowerSnpFilter() {
		return lowerSnpFilter;
	}

	/**
	 * @param lowerSnpFilter the lowerSnpFilter to set
	 */
	public void setLowerSnpFilter(Integer lowerSnpFilter) {
		this.lowerSnpFilter = lowerSnpFilter;
	}

	/**
	 * @return the upperSnpFilter
	 */
	public Integer getUpperSnpFilter() {
		return upperSnpFilter;
	}

	/**
	 * @param upperSnpFilter the upperSnpFilter to set
	 */
	public void setUpperSnpFilter(Integer upperSnpFilter) {
		this.upperSnpFilter = upperSnpFilter;
	}
}
