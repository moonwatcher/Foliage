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

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import org.apache.solr.util.BitSetIterator;
import org.apache.solr.util.OpenBitSet;

import sanger.argml.io.XmlOutput;


/**
 * A finite domain of natural numbers.
 * A <code>NaturalDomain</code> is defined over a {@link NaturalSpace} and does not have to be continuous.
 * Access to <code>min</code>, <code>max</code>, <code>map</code>, <code>one</code>, <code>cardinality</code> and <code>last</code> values 
 * can be concidered to be of <code>O(1)</code> complexity as they are internaly cached.  
 * @author Lior Galanti
 *
 */
public class NaturalDomain extends AbstractSet implements Iterable<Integer>{
	private static final Pattern codec = Pattern.compile("^\\[((([0-9]+),([0-9]+))|([01]+))\\]:(\\[[0-9]+,[0-9]+\\])$");

	private NaturalSpace space;
	
	private Integer _min;
	private Integer _max;
	private OpenBitSet _map;
	private OpenBitSet _one;
	private Integer _cardinality;
	private Integer _last;
	
	private void clearCache(){
		_min = null;
		_max = null;
		_map = null;
		_one = null;
		_cardinality = null;
		_last = null;
	}
	
	/**
	 * Constructs an empty <code>NaturalDomain</code> over <code>space</code>.
	 * @param space {@link NaturalSpace} containing the domain
	 */
	public NaturalDomain(NaturalSpace space){
		super();
		this.space = space;
	}
	
	/**
	 * Constructs the <code>NaturalDomain</code> <code>[0, cardinality -1]</code> over a new space of similar coordinates.
	 * @param cardinality 
	 * @throws NaturalSetException if the domain is empty
	 */
	public NaturalDomain(int cardinality) throws NaturalSetException{
		super(0, cardinality - 1);
		this.space = new NaturalSpace(0, cardinality - 1);
	}

	/**
	 * Constructs the <code>NaturalDomain</code> <code>[min, max]</code> over a new space of similar coordinates.
	 * @param min 
	 * @param max 
	 * @throws NaturalSetException if the domain is empty
	 */
	public NaturalDomain(int min, int max) throws NaturalSetException{
		super(min, max);
		this.space = new NaturalSpace(min, max);
	}

	/**
	 * Constructs the <code>NaturalDomain</code> <code>[min, max]</code> over <code>space</code>.
	 * The resulting domain is intersected with the space to make sure it is covered by the space.
	 * @param space
	 * @param min
	 * @param max
	 */
	public NaturalDomain(NaturalSpace space, int min, int max) {
		super(Math.max(min, space.min()) , Math.min(max, space.max()));
		this.space = space;
	}

	/**
	 * Constructs a <code>NaturalDomain</code> from the {@link OpenBitSet} over <code>space</code>.
	 * The provided <code>map</code> is decoded in the <code>space</code> relative coordinate system.
	 * @param space
	 * @param map
	 */
	public NaturalDomain(NaturalSpace space, OpenBitSet map){
		super(map);
		this.space = space;
		normalize();
	}

	/**
	 * Constructs a sub domain from the provided {@link NaturalSet}.  
	 * @param naturalSet
	 * @return a new <code>NaturalDomain</code> object.
	 * @throws NaturalSetException if the <code>NaturalSet</code> does not belong to the domain.
	 */
	public NaturalDomain createSubDomain(NaturalSet naturalSet) throws NaturalSetException{
		NaturalDomain result;
		if(equals(naturalSet.domain)){
			OpenBitSet sub = new OpenBitSet(space.cardinality());
			BitSetIterator iterate = new BitSetIterator(naturalSet.map());
			
			int dp = iterate.next();
			
			while(dp > -1){
				sub.fastSet(space.toRelativeCoordinate(toAbsoluteCoordinate(dp)));
				dp = iterate.next();
			}
			result = new NaturalDomain(space, sub);
		} else {
			throw new NaturalSetException("Cannot derive a sub-domain from a set of a diffrent domain");
		}
		
		return result;
	}
	
	/**
	 * Constructs sub domain with random selected elements.
	 * The method is non deterministic and sub domains are uniformaly distributed.
	 * @param cardinality Number of elements to select for the sub domain.
	 * @return a new <code>NaturalDomain</code> object.
	 * @throws NaturalSetException if requested cardinality is bigger than the domain's cardinality.
	 */
	public NaturalDomain createRandomSubDomain(int cardinality) throws NaturalSetException{
		OpenBitSet mask = new OpenBitSet(space.cardinality());
		if(cardinality <= cardinality()){
			while(mask.cardinality() < cardinality){
				int unit = (int)Math.round(Math.random() * last());
				if(contains(toRelativeCoordinate(unit))){
					mask.set(unit);
				}
			}
		} else {
			throw new NaturalSetException("Not enough elements in space");
		}
		return new NaturalDomain(space, mask);
	}
	
	/**
	 * Constructs a sub domain with random selected elements, non of which contained in <code>foreign</code> domain.
	 * The method is non deterministic and sub domains are uniformaly distributed.
	 * @param foreign <code>NaturalDomain</code> contaning elements that should not appear in the returned domain.
	 * @param cardinality Number of elements to select for the sub domain.
	 * @return a new <code>NaturalDomain</code> object.
	 * @throws NaturalSetException if requested cardinality is bigger than the <code>cardinality - foreign.cardinality</code>.
	 */
	public NaturalDomain createRandomSubDomain(NaturalDomain foreign, int cardinality) throws NaturalSetException{
		OpenBitSet mask = new OpenBitSet(space.cardinality());
		if(cardinality <= cardinality() - foreign.cardinality()){
			while(mask.cardinality() < cardinality){
				int unit = (int)Math.round(Math.random() * last());
				int relative = toRelativeCoordinate(unit);
				if(contains(relative) && !foreign.contains(relative)){
					mask.set(unit);
				}
			}
		} else {
			throw new NaturalSetException("Not enough elements in space");
		}
		return new NaturalDomain(space, mask);
	}
	
	/**
	 * Creates a copy of the <code>NaturalDomain</code> using the same space.
	 * @return a new <code>NaturalDomain</code> object.
	 */
	public NaturalDomain clone(){
		NaturalDomain clone;
		if(isContinuous()){
			clone = new NaturalDomain(space, min, max);
		} else {
			clone = new NaturalDomain(space, (OpenBitSet)map.clone());
		}
		return clone;
	}

	
	public int min(){
		if(_min == null){
			try {
				if(isContinuous()){
					_min = min;
				} else {
					_min = space.toAbsoluteCoordinate(map.nextSetBit(0));			
				}
			} catch (NaturalSetException e) {
				e.printStackTrace();
			}
		}
		return _min;
	}
	
	public int max(){
		if(_max == null){
			if(isContinuous()){
				_max = max;
			} else {
				try {
					BitSetIterator iterate = new BitSetIterator(map);			
					int next = iterate.next();
					int last = next;
					while(next != -1){
						last = next;
						next = iterate.next();
					}
								
					_max = space.toAbsoluteCoordinate(last);
				} catch (NaturalSetException e) {
					e.printStackTrace();
				}
			}
		}
		return _max;
	}

	public OpenBitSet map(){
		if(_map == null){
			if(isContinuous()){
				if(!isEmpty()){
					try {
						_map = new OpenBitSet(space.cardinality());
						_map.flip(space.toRelativeCoordinate(min), space.toRelativeCoordinate(max)+1);
					} catch (Exception e) {
						e.printStackTrace();
					}
					
				} else {
					_map = new OpenBitSet(space.cardinality());
				}
			} else {
				_map = (OpenBitSet)map.clone();
			}			
		}
		
		return _map;
	}
	

	public boolean contains(int position) {
		boolean result = false;
		try {
			if(!isEmpty()){
				if(isContinuous()){
					result = (min <= position && position <= max);
				} else {
					result = (space.contains(position) && map.fastGet(space.toRelativeCoordinate(position)));
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return result;
	}
	
	public int cardinality() {
		if(_cardinality == null){
			_cardinality = 0;
			if(!isEmpty()) {
				if(isContinuous()) { 
					_cardinality = max - min + 1;
				} else {
					_cardinality = (int)map.cardinality();
				}
			}
		}
		return _cardinality;
	}	
	
	protected void vectorize(){
		if(isContinuous()) {
			if(!isEmpty()){
				map = new OpenBitSet(space.cardinality());
				try {
					map.flip(space.toRelativeCoordinate(min), space.toRelativeCoordinate(max) + 1);
				} catch (NaturalSetException e) {
					e.printStackTrace();
				}
				min = Integer.MIN_VALUE;
				max = Integer.MIN_VALUE;
				
			} else {
				map = new OpenBitSet(space.cardinality());
			}
		}
	}
		
	protected void normalize(){
		if(!isContinuous()){
			int first = map.nextSetBit(0);
			if(first == -1) {
				clear();
			
			} else {
				try {
					OpenBitSet inversed = ((OpenBitSet)map.clone());
					inversed.xor(space.one());
					int end = inversed.nextSetBit(first);				
					if(end == -1 || map.nextSetBit(end) == -1){
						min = space.toAbsoluteCoordinate(first);
						max = space.toAbsoluteCoordinate((end == -1) ? space.last() : end - 1);
						map = null;
					}
				} catch (NaturalSetException e) {
					e.printStackTrace();
				}
			}
		}
	}

	protected OpenBitSet vector() {
		return (OpenBitSet)map().clone();
	}

	/**
	 * Tests for the two <code>NaturalDomain</code> objects are equal.
	 * Two <code>NaturalDomain</code> objects are concidered equal if they are defined over the same space and contain the same elements.
	 * @return True ifd both <code>NaturalDomain</code> objects are equal.
	 */
	public boolean equals(Object o) {
		boolean result = false;
		if (o instanceof NaturalDomain) {
			NaturalDomain other = (NaturalDomain) o;
			if(space.equals(other.space)){
				if(isContinuous() && other.isContinuous()){
					result = (min == other.min && max == other.max);
					
				} else {
					if(map().equals(other.map())) { 
						result = true; 
					}					
				}
			}
		} 
		return result;
	}
		
	/**
	 * @return the <code>NaturalSpace</code> over which this domain is defined.
	 */
	public NaturalSpace space(){
		return space;
	}
			
	public int last(){
		if(_last == null){
			_last = max() - min();
		}
		return _last;
	}
		
	/**
	 * Intersects the domain with the <code>other</code> domain.
	 * @param other The other <code>NaturalDomain</code> to intersect with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 */
	public NaturalDomain intersect(final NaturalDomain other){
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
		clearCache();
		return this;		
	}
	
	/**
	 * Union the domain with the <code>other</code> domain.
	 * @param other The other <code>NaturalDomain</code> to union with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 */
	public NaturalDomain union(final NaturalDomain other){
		if(!other.isEmpty()){
			if(!isEmpty()){
				this.vectorize();
				map.union(other.map());
				this.normalize();
			} else {
				copy(other);
			}
		}
		clearCache();
		return this;
	}

	/**
	 * Xor the domain with the <code>other</code> domain.
	 * @param other The other <code>NaturalDomain</code> to xor with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 */
	public NaturalDomain xor(final NaturalDomain other){
		this.vectorize();
		map.xor(other.map());
		this.normalize();
		clearCache();
		return this;
	}

	/**
	 * Construct a new <code>NaturalDomain</code> composed of the union of the two domain.
	 * @param a
	 * @param b
	 * @return A new <code>NaturalDomain</code> composed of the union of the two domains.
	 * @throws NaturalSetException if the two sets are not defined over the same space.
	 */
	public static NaturalDomain union(final NaturalDomain a, final NaturalDomain b) throws NaturalSetException{
		NaturalDomain domain = a.clone();
		domain.union(b);
		return domain;
	}	
	
	/**
	 * Construct a new <code>NaturalDomain</code> composed of the intersect of the two domain.
	 * @param a
	 * @param b
	 * @return A new <code>NaturalDomain</code> composed of the intersect of the two domains.
	 * @throws NaturalSetException if the two sets are not defined over the same space.
	 */
	public static NaturalDomain intersect(final NaturalDomain a, final NaturalDomain b) throws NaturalSetException{
		NaturalDomain domain = a.clone();
		domain.intersect(b);
		return domain;
	}	
	
	/**
	 * Construct and empty <code>NaturalSet</code> from this domain. 
	 * @return a new empty <code>NaturalSet</code> object.
	 */
	public NaturalSet createEmptyNaturalSet(){
		return new NaturalSet(this);
	}

	/**
	 * Constructs a <code>NaturalSet</code> containing all the elements in the domain.
	 * @return a new <code>NaturalSet</code> object.
	 */
	public NaturalSet createCompleteNaturalSet(){
		NaturalSet naturalSet = new NaturalSet(this, (OpenBitSet)one().clone());
		naturalSet.normalize();
		return naturalSet;
	}

	/**
	 * Constructs a <code>NaturalSet</code> with random selected elements.
	 * The method is non deterministic and <code>NaturalSet</code> are uniformaly distributed.
	 * @param cardinality Number of elements to select for the <code>NaturalSet</code>.
	 * @return a new <code>NaturalSet</code> object.
	 * @throws NaturalSetException if requested cardinality is bigger than the domain's cardinality.
	 */
	public NaturalSet createRandomNaturalSet(int cardinality) throws NaturalSetException{
		OpenBitSet mask = new OpenBitSet(closureCardinality());
		if(cardinality <= cardinality()){
			while(mask.cardinality() < cardinality){
				int unit = (int)Math.round(Math.random() * closureCardinality());
				if(contains(toRelativeCoordinate(unit))){
					mask.set(unit);
				}
			}
		}
		return new NaturalSet(this, mask);
	}

	/**
	 * Translates the position to the relative coordinate system.
	 * In the relative coordinate system the smallest element is always 0.
	 * @param position a position in the absolute coordinate system.
	 * @return the relative position.
	 * @throws NaturalSetException if the position is not contained in the domain.
	 */
	public int toRelativeCoordinate(int position) throws NaturalSetException{
		int value = position - min();
		if(!contains(position)){ throw new NaturalSetException("Coordinate not in domain"); }
		return value;
	}
			
	/**
	 * Translates the position to the absolute cooridinate system.
	 * In the absolute coordinate system the smallest element is the <code>min</code>.
	 * @param position a position in the relative coordinate system.
	 * @return the absolute position.
	 * @throws NaturalSetException if the position is not contained in the domain.
	 */
	public int toAbsoluteCoordinate(int position) throws NaturalSetException{
		int value = position + min();
		if(!contains(value)){ throw new NaturalSetException("Coordinate not in domain"); }
		return value;
	}
	
	public OpenBitSet one(){
		if(_one == null){
			try {
				_one = new OpenBitSet(closureCardinality());
				BitSetIterator iterate = new BitSetIterator(map());			
				int p = iterate.next();			
				while(p > -1){
					_one.fastSet(toRelativeCoordinate(space.toAbsoluteCoordinate(p)));
					p = iterate.next();
				}
			} catch (NaturalSetException e) {
				e.printStackTrace();
			}
		}		
		return _one;
	}

	/**
	 * Translate <code>source</code> to this domain's coordinates.
	 * @param source The <code>NaturalSet</code> to translate.
	 * @return A new <code>NaturalSet</code> object in this domain's coordinate system.
	 * @throws NaturalSetException if the <code>source</code> set is of a different <code>NaturalSpace</code>. 
	 */
	public NaturalSet project(final NaturalSet source) throws NaturalSetException{
		NaturalSet project = null;
		if(space.equals(source.domain.space)){
			if(!source.isEmpty()){				
				OpenBitSet projection = new OpenBitSet(closureCardinality());
				BitSetIterator iterate = new BitSetIterator(source.map());
				
				int sdp = iterate.next(Math.max(source.domain.toRelativeCoordinate(min()), 0));
				int sp = 0;
				
				while(sp <= max() && sdp > -1){
					sp = source.domain.toAbsoluteCoordinate(sdp);
					if(contains(sp)){
						projection.fastSet(toRelativeCoordinate(sp));
					}
					sdp = iterate.next();
				}
				project = new NaturalSet(this, projection);
				
			} else {
				project = createEmptyNaturalSet();
			}
		} else {
			throw new NaturalSetException("Cannot project between different spaces");
		}
		
		return project;
	}
	
	public String toBinaryString() {
		StringBuilder sb = new StringBuilder();
		for(int i= 0; i< space.cardinality(); i++) {
			sb.append(map().fastGet(i)==false ? '0' : '1');
		}
		return sb.toString();
	}
	
	public void writeXml(XmlOutput out, String name) throws XMLStreamException{
		out.writer().writeStartElement("", "domain");
		if(name!=null) out.writer().writeAttribute("name", name);
		space.writeXml(out);
		if(!isEmpty()){
			if(isContinuous()){
				out.writer().writeStartElement("", "continuous");
				out.writer().writeAttribute("min", String.valueOf(min()));
				out.writer().writeAttribute("max", String.valueOf(max()));
				out.writer().writeEndElement();
			} else {
				out.writer().writeStartElement("", "fragments");
				out.writer().writeAttribute("map", toBinaryString());
				out.writer().writeEndElement();
			}
		}
		out.writer().writeEndElement();		
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("[");
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
		sb.append(space);

		return sb.toString();
	}
	
	/**
	 * @param value a string encoding the domain.
	 * @return a <code>NaturalDomain</code> object decoded from <code>value</code>.
	 * @throws NaturalSetException if the string cannot be decoded to a <code>NaturalDomain</code>.
	 */
	public static NaturalDomain decode(String value) throws NaturalSetException{
		NaturalDomain result = null;
		Matcher matcher = codec.matcher(value);
		if(matcher.matches()){
			NaturalSpace s = NaturalSpace.decode(matcher.group(6));
			if(matcher.group(5) != null){
				 result = new NaturalDomain(s, NaturalSpace.decodeBinaryString(matcher.group(5)));
			} else {
				result = new NaturalDomain(s, Integer.decode(matcher.group(3)), Integer.decode(matcher.group(4)));				
			}
		} else {
			throw new NaturalSetException("String does not decode to a Natural Domain");
		}
		return result;
	}	

	public static NaturalDomain parser(XMLStreamReader parser) throws XMLStreamException, NaturalSetException{
		NaturalDomain domain = null;
		int smin=0, smax=0, dmin=0, dmax=0;
		OpenBitSet map = null;
		for(int event=parser.next(); domain==null; event=parser.next()) {
			switch (event) {
				case XMLStreamConstants.START_ELEMENT:
					if(parser.getLocalName().equals("space")) {
						smin = Integer.parseInt(parser.getAttributeValue("", "min"));
						smax = Integer.parseInt(parser.getAttributeValue("", "max"));
						
					} else if (parser.getLocalName().equals("continuous")) {
						dmin = Integer.parseInt(parser.getAttributeValue("", "min"));
						dmax = Integer.parseInt(parser.getAttributeValue("", "max"));
						
					} else if (parser.getLocalName().equals("fragments")) {
						map = NaturalSpace.decodeBinaryString(parser.getAttributeValue("", "map"));
					}
					break;

				case XMLStreamConstants.END_ELEMENT:
					if (parser.getLocalName().equals("domain")) {
						if(map == null){
							domain = new NaturalDomain(new NaturalSpace(smin, smax), dmin, dmax); 
						} else {
							domain = new NaturalDomain(new NaturalSpace(smin, smax),map); 							
						}
					}
				break;
			}
		}
		return domain;
	}
	
	public Iterator<Integer> iterator(){
		return new DomainIterator(this);
	}
	
	public class DomainIterator implements Iterator<Integer> {
		private NaturalDomain domain;
		private BitSetIterator iterate;
		private int next;

		public DomainIterator(NaturalDomain domain) {
			this.domain = domain;			
			this.iterate = new BitSetIterator(domain.map());
			this.next = iterate.next();			
		}
		
		public boolean hasNext() {
			return next!=-1;
		}

		public Integer next() {
			int value;
			try { value = domain.space.toAbsoluteCoordinate(next); } 
			catch (NaturalSetException e) { 
				throw new NoSuchElementException("No more elements in domain."); }
			next = iterate.next();
			return value;
		}

		public void remove() {}
	}	
}
