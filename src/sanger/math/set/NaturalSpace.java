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

package sanger.math.set;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import org.apache.solr.util.OpenBitSet;

import sanger.argml.io.XmlOutput;

/**
 * A continuous and finit space of natural numbers.
 * A <code>NaturalSpace</code> is always continuous and provides the facility to define {@link NaturalDomain} objects contained in it.
 * It provides {@link NaturalDomain} objects with facilities for move between the relative and absolute coordinate system 
 * so that memory complexity of <code>NaturalSpace</code> does not depand on its <code>min</code> and <code>max</code> values
 * but only on its <code>cardinality</code>. 
 * @author Lior Galanti
 */
public class NaturalSpace {
	private static final Pattern codec = Pattern.compile("^\\[([0-9]+),([0-9]+)\\]$");

	private int min;
	private int max;
	private int cardinality;
	private OpenBitSet one;
	private int last;
	
	/**
	 * Constructs a <code>NaturalSpace</code> spanning from <code>min</code> to <code>max</code>.
	 * A <code>NaturalSpace</code> cannot be empty.
	 * @param min The smallest element in the space
	 * @param max The largest element in the space
	 * @throws NaturalSetException If the space is empty
	 */
	public NaturalSpace(int min, int max) throws NaturalSetException{
		if(min > Integer.MIN_VALUE){
			this.min = min;
			this.max = max;
			this.last = max - min;
			this.cardinality = last + 1;
			
			this.one = new OpenBitSet(cardinality);
			this.one.flip(0, cardinality);
			if(min > max){
				throw new NaturalSetException("Natural Space can not be empty");
			}
		} else {
			throw new NaturalSetException("Natural Space min too small.");
		}
	}
	
	/**
	 * Translates <code>position</code> to the relative coordinate system.
	 * In the relative coordinate system the smallest element is always <code>0</code>.
	 * @param position A natural number in the absolute coordinate system.
	 * @return The relative position.
	 * @throws NaturalSetException If the <code>position</code> is not contained in the space.
	 */
	protected int toRelativeCoordinate(int position) throws NaturalSetException{
		int value = position - min;
		if(!contains(position)){ throw new NaturalSetException("Coordinate not in Space"); }
		return value;
	}
			
	/**
	 * Translates <code>position</code> to the absolute cooridinate system.
	 * In the absolute coordinate system the smallest element is the <code>min</code>.
	 * @param position A natural number in the relative coordinate system.
	 * @return The absolute position.
	 * @throws NaturalSetException If the <code>position</code> is not contained in the space.
	 */
	protected int toAbsoluteCoordinate(int position) throws NaturalSetException{
		int value = position + min;
		if(!contains(value)){ throw new NaturalSetException("Coordinate not in Space"); }
		return value;
	}
		
	/**
	 * @return The largest element in the space.
	 */
	public int max() {
		return max;
	}

	/**
	 * @return The smallest element in the space.
	 */
	public int min() {
		return min;
	}

	/** 
	 * @return True if the two spaces contain the same set of natural numbers.
	 */
	public boolean equals(Object o) {
		boolean result = false;
		if (o instanceof NaturalSpace) {
			NaturalSpace other = (NaturalSpace) o;
			if(min == other.min && max == other.max) { 
				result = true; 
			}
		}
		return result;
	}
	
	/**
	 * @return The largest element in the relative coordinates.
	 */
	public int last(){
		return last;
	}
	
	/**
	 * @return An {@link OpenBitSet} with elements in the space.
	 * The <code>one</code> vector is in the relative coordinate system, 
	 * so the first bit is the <code>min<code> element. 
	 * This is useful when performing binary operations on sets in the space. 
	 */
	public OpenBitSet one(){
		return one;
	}
	
	/**
	 * @return The cardinality of the space.
	 */
	public int cardinality(){
		return cardinality;
	}
	
	/**
	 * Tests if a number is contained in the space.
	 * @param position
	 * @return True if the position is contained in the space.
	 */
	public boolean contains(int position){
		return (min <= position && position <= max);
	}

	/**
	 * @return a {@link NaturalDomain} containing all the elements in the space.
	 */
	public NaturalDomain createCompleteDomain(){
		return new NaturalDomain(this, min, max);
	}
	
	public void writeXml(XmlOutput out) throws XMLStreamException{
		out.writer().writeStartElement("space");
		out.writer().writeAttribute("min", String.valueOf(min()));
		out.writer().writeAttribute("max", String.valueOf(max()));
		out.writer().writeEndElement();
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		sb.append(min);
		sb.append(',');
		sb.append(max);			
		sb.append(']');
		return sb.toString();
	}
	
	/**
	 * @param value a string encoding the space
	 * @return a <code>NaturalSpace</code> object decoded from the string
	 * @throws NaturalSetException if the string cannot be decoded to a <code>NaturalSpace</code>
	 */
	public static NaturalSpace decode(String value) throws NaturalSetException{
		NaturalSpace result = null;
		Matcher matcher = codec.matcher(value);
		if(matcher.matches()){
			result = new NaturalSpace(Integer.decode(matcher.group(1)), Integer.decode(matcher.group(2)));
		} else {
			throw new NaturalSetException("String does not decode to a Natural Space");
		}
		return result;
	}
	
	/**
	 * @param value string to decode
	 * @return an {@link OpenBitSet} in the relative coordinate system
	 */
	public static OpenBitSet decodeBinaryString(String value){
		OpenBitSet bitSet = new OpenBitSet(value.length());
		for(int i=0; i<value.length(); i++){
			if('1' == value.charAt(i)){bitSet.fastSet(i);}
		} 
		return bitSet;
	}
}
