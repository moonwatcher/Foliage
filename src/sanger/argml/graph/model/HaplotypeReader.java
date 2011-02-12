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
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sanger.argml.environment.Environment;
import sanger.argml.io.TextInput;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSetException;

public class HaplotypeReader implements Iterable<Haplotype> {
	private TextInput input;
	private enum HaplotypeFileLineType {HEADER, MARKER, HAPLOTYPE, EOF, CLOSED}	
	private HaplotypeFileLineType state;

	private static final Pattern headerPattern = Pattern.compile("^([0-9]+) ([0-9]+) ([0-9]+)$");
	private static final Pattern markerPattern = Pattern.compile("^([0-9]+)$");
	private Pattern haplotypePattern = null;
	
	private NaturalDomain haplotypeDomain = null;	
	private NaturalDomain snpDomain = null;
	private NaturalDomain basePairDomain = null;
	
	private int[] markerPositions = null;
	private int mindex = 0;
	private Haplotype currentHaplotype = null;
	
	public HaplotypeReader(Environment env, TextInput input) throws FileNotFoundException{
		this.input = input;
		this.state = HaplotypeFileLineType.HEADER;
		while (state != HaplotypeFileLineType.HAPLOTYPE && state != HaplotypeFileLineType.EOF){
			readNext();
		}
	}
	
	public void close() throws IOException {
		input.close();
		state = HaplotypeFileLineType.CLOSED;
		haplotypePattern = null;
		haplotypeDomain = null;
		snpDomain = null;
		basePairDomain = null;
		markerPositions = null;
		currentHaplotype = null;
		mindex = 0;
	}
		
	public HaplotypeSet readHaplotypeSet() throws NaturalSetException{
		HaplotypeSet hs = new HaplotypeSet(input.env(), haplotypeDomain, snpDomain, markerPositions);
		for(Haplotype haplotype : this){
			hs.add(haplotype);
		}
		return hs;
	}
	
	private void processHeader(Matcher matcher) throws NaturalSetException{
		int markers = Integer.decode(matcher.group(3));		
		snpDomain = new NaturalDomain(markers);
		
		int haplotypes = Integer.decode(matcher.group(1)) + Integer.decode(matcher.group(2));		
		haplotypeDomain = new NaturalDomain(haplotypes);
		
		haplotypePattern = Pattern.compile("^[01M]{" + markers + "}$");		
		state = HaplotypeFileLineType.MARKER;		
		markerPositions = new int[markers];
	}
	
	// Haplotype iterator
	
	private void readNext() {
		try {	
			String line = input.reader().readLine();
			if(line != null){
				Matcher matcher = headerPattern.matcher(line);
				if(state == HaplotypeFileLineType.HEADER && matcher.matches()){
					processHeader(matcher);
					
				} else if (state != HaplotypeFileLineType.HEADER) {
					if(state != HaplotypeFileLineType.HAPLOTYPE) {
						matcher = markerPattern.matcher(line);
						if(matcher.matches()){
							markerPositions[mindex]= Integer.decode(matcher.group(1));
							mindex++;
							if(mindex == snpDomain.cardinality()){
								state = HaplotypeFileLineType.HAPLOTYPE;
								basePairDomain = new NaturalDomain(markerPositions[0], markerPositions[markerPositions.length -1]);
								readNext();
							}
						} 
					} else {
						matcher = haplotypePattern.matcher(line);
						if(matcher.matches()){
							currentHaplotype = Haplotype.decode(line, snpDomain);
						}
					}
				}
			} else {
				state = HaplotypeFileLineType.EOF;
				input.reader().close();
			}
		} catch (Exception e) {
			input.env().log().printError(e);
		}
	}		

	public Iterator<Haplotype> iterator(){
		return new HaplotypeIterator(this);
	}
	
	public class HaplotypeIterator implements Iterator<Haplotype> {
		HaplotypeReader factory;

		public HaplotypeIterator(HaplotypeReader factory) {
			this.factory = factory;
		}
		
		public boolean hasNext() {
			return factory.state == HaplotypeFileLineType.HAPLOTYPE;
		}

		public Haplotype next() {
			Haplotype result = null;
			if(hasNext()){
				result = factory.currentHaplotype;
				factory.readNext();
			} else {
				throw new NoSuchElementException("No more haplotypes in file.");
			}
			return result;
		}

		public void remove() {}
	}

	// Accessor functions
	
	/**
	 * @return the markerPositions
	 */
	public int[] getMarkerPositions() {
		return markerPositions;
	}

	/**
	 * @return the basePairDomain
	 */
	public NaturalDomain getBasePairDomain() {
		return basePairDomain;
	}
	
	/**
	 * @return the haplotypeDomain
	 */
	public NaturalDomain getHaplotypeDomain() {
		return haplotypeDomain;
	}
	
	/**
	 * @return the snpDomain
	 */
	public NaturalDomain getSnpDomain() {
		return snpDomain;
	}
	
}
