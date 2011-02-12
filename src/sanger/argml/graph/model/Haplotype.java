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

import org.apache.solr.util.OpenBitSet;

import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetException;

/**
 * A biallelic haplotype.
 * The {@link #getMissing() missing} set indicates which values are in fact unknown
 * and should be ignored when reading the {@link #getMarkers() marker} values set.  
 * @author lg8
 *
 */
public class Haplotype {
	private NaturalSet markers;
	private NaturalSet missing;

	/**
	 * Construct a new haplotype with the supplied markers and missing values.
	 * @param markers Biallelic marker values.
	 * @param missing Positions where the value is unknown.
	 */
	protected Haplotype(NaturalSet markers, NaturalSet missing){
		this.markers = markers;
		this.missing = missing;
	}
	
	/**
	 * Decode a binary string into a haplotype.
	 * line should match the <code>[01M]*<code> expression.
	 * @param line String to decode
	 * @param domain Domain determining the dimensions of the haplotype.
	 * @return A new haplotypes.
	 */
	public static Haplotype decode(String line, NaturalDomain domain){
		return new Haplotype(decodeMarkers(line, domain), decodeMissing(line, domain));
	}

	/**
	 * Decode marker values from the string.
	 * @param line
	 * @param domain
	 * @return The markers
	 */
	private static NaturalSet decodeMarkers(String line, NaturalDomain domain){
		OpenBitSet bitSet = new OpenBitSet(domain.closureCardinality());
		for(int i=0; i<line.length(); i++){
			if('1' == line.charAt(i)) bitSet.fastSet(i);
		}
		return new NaturalSet(domain, bitSet); 
	}

	/**
	 * Decode the missing values from the string.
	 * @param line
	 * @param domain
	 * @return The missing positions
	 */
	private static NaturalSet decodeMissing(String line, NaturalDomain domain){
		OpenBitSet bitSet = new OpenBitSet(domain.closureCardinality());
		for(int i=0; i<line.length(); i++){
			if('M' == line.charAt(i)) bitSet.fastSet(i);
		} 
		return new NaturalSet(domain, bitSet); 
	}
	
	/**
	 * Project the haplotype onto a different domain.
	 * For more details about projecting a set see the domain 
	 * {@link NaturalDomain#project(NaturalSet) project} method.
	 * @param domain Domain to prject the haplotype to.
	 * @param haplotype Haplotype to project.
	 * @return A new haplotype in <code>domain</code> coordinates.
	 * @throws NaturalSetException if the haplotype is not defined over the same space as <code>domain</code>.
	 */
	public static Haplotype project(NaturalDomain domain, Haplotype haplotype) throws NaturalSetException{
		return new Haplotype(domain.project(haplotype.markers), domain.project(haplotype.missing));
	}

	/**
	 * @return the markers
	 */
	public NaturalSet getMarkers() {
		return markers;
	}

	/**
	 * @return the missing
	 */
	public NaturalSet getMissing() {
		return missing;
	}

	/**
	 * @return The haplotype's domain.
	 */
	public NaturalDomain domain(){
		return markers.domain();
	}
}
