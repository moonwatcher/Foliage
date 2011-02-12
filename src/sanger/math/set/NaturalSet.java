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

import org.apache.solr.util.BitSetIterator;
import org.apache.solr.util.OpenBitSet;

/**
 * A finite set of natural numbers.
 * @author Lior Galanti
 *
 */
public class NaturalSet extends AbstractSet implements Comparable<NaturalSet> {
	private static final Pattern codec = Pattern.compile("^\\[((([0-9]+),([0-9]+))|([01]+))\\]:(\\[(([0-9]+,[0-9]+)|([01]+))\\]:\\[[0-9]+,[0-9]+\\])$");
	protected NaturalDomain domain;

	/**
	 * Constructs an empty <code>NaturalSet</code>.
	 * @param domain
	 */
	protected NaturalSet(NaturalDomain domain){
		super();
		this.domain = domain;
	}
	
	/**
	 * Constructs the <code>NaturalSet</code> <code>[min, max]</code> over <code>domain</code>.
	 * The resulting set is intersected with the domain to make sure it is covered by the space.
	 * @param domain
	 * @param min
	 * @param max
	 */
	public NaturalSet(NaturalDomain domain, int min, int max){
		super(Math.max(min, domain.min()), Math.min(max, domain.max()));
		this.domain = domain;
		if(!this.domain.isContinuous()){
			vectorize();
			this.map.intersect(domain.one());	
		}
		normalize();
	}
			
	/**
	 * Constructs a <code>NaturalSet</code> from the {@link OpenBitSet} over <code>domain</code>.
	 * The provided <code>map</code> is decoded in the <code>domain</code> relative coordinate system.
	 * @param domain
	 * @param map
	 */
	public NaturalSet(NaturalDomain domain, OpenBitSet map){
		super(map);
		this.domain = domain;
		this.map.intersect(this.domain.one());
		normalize();
	}
	
	/**
	 * Creates a copy of the <code>NaturalSet</code> using the same domain.
	 * @return a new <code>NaturalSet</code> object.
	 */
	public NaturalSet clone(){
		NaturalSet clone;
		if(isContinuous()){
			clone = new NaturalSet(domain, min, max);
		} else {
			clone = new NaturalSet(domain, (OpenBitSet)map.clone());
		}
		return clone;
	}

	public int min(){
		int result = min;
		if(!isContinuous()){
			try {
				result = domain.toAbsoluteCoordinate(map.nextSetBit(0));	
			} catch (NaturalSetException e) {
				e.printStackTrace();
			}					
		}
		return result;
	}
	
	public int max(){
		int result = max;
		if(!isContinuous()){
			BitSetIterator iterate = new BitSetIterator(map);			
			int next = iterate.next();
			int last = next;
			while(next != -1){
				last = next;
				next = iterate.next();
			}
			try {			
				result = domain.toAbsoluteCoordinate(last);
			} catch (NaturalSetException e) {
				e.printStackTrace();
			}				
		}
		return result;
	}
	
	public OpenBitSet map(){
		OpenBitSet vector = null;
		if(isContinuous()){
			if(!isEmpty()){
				vector = new OpenBitSet(domain.closureCardinality());
				try {
					int lmin = domain.toRelativeCoordinate(min);
					int lmax = domain.toRelativeCoordinate(max);
					for(int i = lmin; i <= lmax; i++){
						vector.fastSet(i);
					}
				} catch (NaturalSetException e) {
					e.printStackTrace();
				}	
			} else {
				vector = new OpenBitSet(domain.closureCardinality());
			}
		} else {
			vector = map;
		}
		return vector;			
	}
	
	/**
	 * The order relation is defined as the natural order on the binary number space.
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(NaturalSet other){
		int result = 0;
		try{
			NaturalSet d = NaturalSet.xor(this, other);
			OpenBitSet test = d.map();				
			int msb = test.nextSetBit(0);
			if(msb >= 0){ result = this.contains(domain.toAbsoluteCoordinate(msb))? -1 : 1; }				
		} catch (NaturalSetException e) {
			throw new ClassCastException("Cannot compare sets of diffrent domain");
		}
		return result;
	}

	public boolean contains(int position) {
		boolean result = false;
		if(!isEmpty()){
			if(isContinuous()){
				result = (min <= position && position <= max);
			} else {
				try {
					result = (domain.contains(position) && map.fastGet(domain.toRelativeCoordinate(position)));
				} catch (NaturalSetException e) {
					e.printStackTrace();
				}	
			}
		}
		return result;
	}
				
	public int cardinality() {
		int result = 0;
		if(!isEmpty()) {
			if(isContinuous()) { 
				result = max - min + 1;
			} else {
				result = (int)map.cardinality();
			}
		}
		return result;
	}


	protected OpenBitSet vector() {
		return (OpenBitSet)map().clone();
	}
	
	protected void vectorize(){
		if(isContinuous()) {
			if(!isEmpty()){
				map = new OpenBitSet(domain.closureCardinality());
				try {
					int lmin = domain.toRelativeCoordinate(min);
					int lmax = domain.toRelativeCoordinate(max);
					for(int i = lmin; i <= lmax; i++){
						map.fastSet(i);
					}
					min = Integer.MIN_VALUE;;
					max = Integer.MIN_VALUE;;
				} catch (NaturalSetException e) {
					e.printStackTrace();
				}	
			} else {
				map = new OpenBitSet(domain.closureCardinality());
			}
		}
	}
		
	protected void normalize(){
		if(!isContinuous()){
			int first = map.nextSetBit(0);
			if(first == -1) {
				clear();
			
			} else {
				OpenBitSet inversed = ((OpenBitSet)map.clone());
				inversed.xor(domain.one());
				int end = inversed.nextSetBit(first);					
				if(end == -1 || map.nextSetBit(end) == -1){
					try {
						min = domain.toAbsoluteCoordinate(first);
						max = domain.toAbsoluteCoordinate((end == -1) ? domain.last() : end - 1);
						map = null;
					} catch (NaturalSetException e) {
						e.printStackTrace();
					}	
				} 
			}
		}
	}
		
	public static NaturalSet inverse(NaturalSet o){
		OpenBitSet inv = o.vector();
		inv.xor(o.domain().one());
		return new NaturalSet(o.domain(), inv);
	}
	
	public NaturalSet inverse(){
		if(!isEmpty()){
			this.vectorize();
			map.xor(domain().one());
			this.normalize();
		} else {
			map = (OpenBitSet)domain.one().clone();
		}
		return this;
	}
	
	
	
	public OpenBitSet toOpenBitSet(){
		OpenBitSet vector = null;
		if(isContinuous()){
			if(!isEmpty()){
				vector = new OpenBitSet(domain.closureCardinality());
				try {
					int lmin = domain.toRelativeCoordinate(min);
					int lmax = domain.toRelativeCoordinate(max);
					for(int i = lmin; i <= lmax; i++){
						vector.fastSet(i);
					}
				} catch (NaturalSetException e) {
					e.printStackTrace();
				}	
			} else {
				vector = new OpenBitSet(domain.closureCardinality());
			}
		} else {
			vector = (OpenBitSet)map.clone();
		}
		return vector;			
	}
			
	/**
	 * Two <code>NaturalSet</code> objects are equal if they contain the same elements.
	 * Only sets from the same domain can be equal. 
	 * To compare sets of different domain use {@link NaturalDomain#project(NaturalSet)} first.
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	public boolean equals(Object o) {
		boolean result = false;
		if (o instanceof NaturalSet) {
			NaturalSet other = (NaturalSet) o;
			if(this.domain.equals(other.domain)){
				if(isContinuous() && other.isContinuous()){
					result = (min == other.min && max == other.max);
				} else {
					result = map().equals(other.map());
				}
			}
		}
		return result;
	}

	
	/**
	 * Union with <code>other</code> and return the number of add elements.
	 * @param other The <code>NaturalSet</code> to union with.
	 * @return The number of elements that were added to the set.
	 * Strictly: <code>other.cardinality() - NaturalSet.intersect(this, other)</code>
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public int unionAndCompare(final NaturalSet other) throws NaturalSetException{
		int compare = cardinality();
		union(other);
		compare = cardinality() - compare;
		return compare;
	}
	
	/**
	 * Count the number of elements in the intersect with <code>other</code>.
	 * This does not modify the set.
	 * @param other The other <code>NaturalSet</code>.
	 * @return The number of elements in the intersect with <code>other</code>.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public int intersectionCount(final NaturalSet other) throws NaturalSetException{
		return NaturalSet.intersect(this, other).cardinality();
	}
			
	/**
	 * Count the number of elements in the union with <code>other</code>.
	 * This does not modify the set.
	 * @param other The other <code>NaturalSet</code>.
	 * @return The number of elements in the union with <code>other</code>.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public int unionCount(final NaturalSet other) throws NaturalSetException{
		return NaturalSet.union(this, other).cardinality();
	}
			
	/**
	 * Intersect this set with <code>other</code>.
	 * @param other the set to intersect with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public NaturalSet intersect(final NaturalSet other) throws NaturalSetException{
		if(domain.equals(other.domain)){
			if(isEmpty() || other.isEmpty()) {
				clear();
			} else {
				if(isContinuous() && other.isContinuous()) {
					min =  Math.max(min, other.min);
					max = Math.min(max, other.max);
					if(min > max) { clear(); }
					
				} else {
					this.vectorize();
					map.intersect(other.map());
					this.normalize();
				}
			}
		} else {
			throw new NaturalSetException("Cannot intersect sets of different domains");
		}
		return this;		
	}
	
	/**
	 * Union this set with <code>other</code>.
	 * @param other the set to union with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public NaturalSet union(final NaturalSet other) throws NaturalSetException {
		if(domain.equals(other.domain)){
			if(!other.isEmpty()){
				if(!isEmpty()){
					this.vectorize();
					map.union(other.map());
					this.normalize();
				} else {
					copy(other);
				}
			}
		} else {
			throw new NaturalSetException("Cannot union sets of different domains");
		}
		return this;
	}

	/**
	 * XOR this set with  <code>other</code>.
	 * @param other the set to XOR with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public NaturalSet xor(final NaturalSet other) throws NaturalSetException {
		if(domain.equals(other.domain)){
			this.vectorize();
			map.xor(other.map());
			this.normalize();
		} else {
			throw new NaturalSetException("Cannot XOR sets of different domains");
		}
		return this;
		
	}		
	
	/**
	 * AND this set with  <code>other</code>.
	 * @param other the set to AND with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public NaturalSet and(final NaturalSet other) throws NaturalSetException {
		if(domain.equals(other.domain)){
			this.vectorize();
			map.and(other.map());
			this.normalize();
		} else {
			throw new NaturalSetException("Cannot AND sets of different domains");
		}
		return this;
		
	}		
	
	/**
	 * OR this set with  <code>other</code>.
	 * @param other the set to OR with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public NaturalSet or(final NaturalSet other) throws NaturalSetException {
		if(domain.equals(other.domain)){
			this.vectorize();
			map.or(other.map());
			this.normalize();
		} else {
			throw new NaturalSetException("Cannot OR sets of different domains");
		}
		return this;
		
	}		
	
	
	/**
	 * Construct a new <code>NaturalSet</code> composed of the intersect of the two sets.
	 * @param a
	 * @param b
	 * @return A new <code>NaturalSet</code> composed of the intersect of the two sets.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public static NaturalSet intersect(final NaturalSet a, final NaturalSet b) throws NaturalSetException{
		NaturalSet naturalSet = a.clone();
		naturalSet.intersect(b);
		return naturalSet;
	}
	
	/**
	 * Construct a new <code>NaturalSet</code> composed of the union of the two sets.
	 * @param a
	 * @param b
	 * @return A new <code>NaturalSet</code> composed of the union of the two sets.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public static NaturalSet union(final NaturalSet a, final NaturalSet b) throws NaturalSetException{
		NaturalSet naturalSet = a.clone();
		naturalSet.union(b);
		return naturalSet;
	}
	
	/**
	 * Construct a new <code>NaturalSet</code> composed of the xor of the two sets.
	 * @param a
	 * @param b
	 * @return A new <code>NaturalSet</code> composed of the xor of the two sets.
	 * @throws NaturalSetException if the two sets are not defined over the same domain.
	 */
	public static NaturalSet xor(final NaturalSet a, final NaturalSet b) throws NaturalSetException{
		NaturalSet naturalSet = a.clone();
		naturalSet.xor(b);
		return naturalSet;
	}
	
	
	/**
	 * Remove all the element strictly smaller than <code>position</code>.
	 * @param position
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 */
	public NaturalSet removeBelow(int position){
		if(!isEmpty() && domain.contains(position)){
			if(position <= max()){
				if(isContinuous()){
					min = Math.max(min, position);
				} else {
					try {
						int lposition = domain.toRelativeCoordinate(position);
						BitSetIterator iterate = new BitSetIterator(map);					
						int i = iterate.next();
						while(i > -1 && i < lposition){
							map.fastClear(i);
							i = iterate.next();
						}
						normalize();
					} catch (NaturalSetException e) {
						e.printStackTrace();
					}						
				}					
			} else {
				clear();
			}
		}
		return this;
	}
	
	/**
	 * Remove all the element strictly bigger than <code>position</code>.
	 * @param position
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 */
	public NaturalSet removeAbove(int position){
		if(!isEmpty() && domain.contains(position)){
			if(position >= min()){
				if(isContinuous()){
					max = Math.min(max, position);
				} else {
					if(position < max()){
						try {
							int lposition = domain.toRelativeCoordinate(position);
							BitSetIterator iterate = new BitSetIterator(map);
							int i = iterate.next(lposition + 1);
							while(i > -1 && i <= domain.max()){
								map.fastClear(i);
								i = iterate.next();
							}
							normalize();
						} catch (NaturalSetException e) {
							e.printStackTrace();
						}	
					}
				}					
			} else {
				clear();
			}
		}
		
		return this;
	}
							
	/**
	 * Construct a new <code>NaturalSet</code> and copy all the elements from this set who are bigger than or equal to <code>position</code>.
	 * @param position
	 * @return A new <code>NaturalSet</code> object.
	 */
	public NaturalSet copyFrom(int position) {
		NaturalSet naturalSet;
		if(!isEmpty()  && domain.contains(position)){
			naturalSet = this.clone();
			naturalSet.removeBelow(position);
		} else {
			naturalSet = new NaturalSet(domain);
		}
		return naturalSet;
	}
	
	
	
	/**
	 * Construct a new <code>NaturalSet</code> and copy all the elements from this set who are smaller than or equal to <code>position</code>.
	 * @param position
	 * @return A new <code>NaturalSet</code> object.
	 */
	public NaturalSet copyUpTo(int position){
		NaturalSet naturalSet;
		if(!isEmpty()  && domain.contains(position)){
			naturalSet = this.clone();
			naturalSet.removeAbove(position);
		} else {
			naturalSet = new NaturalSet(domain);
		}
		return naturalSet;
	}
			
	public String toBinaryString() {
		StringBuilder sb = new StringBuilder();
		OpenBitSet vector = map();
		for(int i=0; i<domain.closureCardinality(); i++) {
			sb.append(vector.fastGet(i)==false ? '0' : '1');
		}
		return sb.toString();
	}
	
					
	public String toString(){
		StringBuilder sb = new StringBuilder();
		
		try {
			sb.append('[');

			if(isContinuous() && !isEmpty()){
				sb.append(min);
				sb.append(',');
				sb.append(max);
			} else if (!isContinuous()) {
				sb.append(toBinaryString());
			} else {
				sb.append(topology());
			}
			
			sb.append("]:");
			sb.append(domain.toString());
		} catch (Exception e) {
			e.printStackTrace();
		}

		return sb.toString();
	}
	
	public String toCompactString(){
		StringBuilder sb = new StringBuilder();
		
		try {
			if(isContinuous() && !isEmpty()){
				sb.append(min);
				sb.append(",");
				sb.append(max);
			} else if (!isContinuous()) {
				sb.append(toBinaryString());
			} else {
				sb.append(topology());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		return sb.toString();
	}
	
	public NaturalDomain domain() {
		return domain;
	}

	/**
	 * @param value a string encoding the set.
	 * @return A <code>NaturalSet</code> object decoded from <code>value</code>.
	 * @throws NaturalSetException if the string cannot be decoded to a <code>NaturalSet</code>
	 */
	public static NaturalSet decode(String value) throws NaturalSetException{
		NaturalSet result = null;
		Matcher matcher = codec.matcher(value);
		if(matcher.matches()){
			NaturalDomain d = NaturalDomain.decode(matcher.group(6));
			if(matcher.group(5) != null){
				 result = new NaturalSet(d, NaturalSpace.decodeBinaryString(matcher.group(5)));
			} else {
				result = new NaturalSet(d, Integer.decode(matcher.group(3)), Integer.decode(matcher.group(4)));				
			}
		} else {
			throw new NaturalSetException("String does not decode to a Natural Domain");
		}
		return result;
	}	
}
