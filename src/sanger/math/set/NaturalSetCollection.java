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

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import org.apache.solr.util.OpenBitSet;



/**
 * A collection of {@link NaturalSet} objects.
 * @author Lior Galanti
 *
 */
public class NaturalSetCollection implements Iterable<NaturalSet> {
	private NaturalDomain domain;
	private ArrayList<NaturalSet> elements;
	
	/**
	 * Constructs a collection with an initial capacity for <code>size</code> elements.
	 * If <code>populate</code> is <code>true</code> the resulting object is populated with empty {@link NaturalSet} objects.
	 * @param domain The <code>NaturalDomain</code> all sets in the collection must be defined on.
	 * @param size Initial capacity of the collection. In case <code>populate</code> is <code>true</code> this is the number of empty sets the collection is populated with.
	 * @param populate Singals if to populate the collection with empty sets.
	 */
	public NaturalSetCollection(NaturalDomain domain, int size, boolean populate){
		this.domain = domain;
		this.elements = new ArrayList<NaturalSet>(size);
		if(populate){
			for(int i=0; i<size; i++){
				elements.add(new NaturalSet(domain, new OpenBitSet(domain.closureCardinality())));
			}
		}
	}
	
	/**
	 * Constructs a collection with an initial capacity for <code>size</code> elements.
	 * @param domain The <code>NaturalDomain</code> all sets in the collection must be defined on.
	 * @param size Initial capacity of the collection.
	 */
	public NaturalSetCollection(NaturalDomain domain, int size){
		this.domain = domain;
		this.elements = new ArrayList<NaturalSet>(size);
	}
	
	/**
	 * Constructs an empty collection.
	 * @param domain The <code>NaturalDomain</code> all sets in the collection must be defined on.
	 */
	public NaturalSetCollection(NaturalDomain domain){
		this.domain = domain;
		this.elements = new ArrayList<NaturalSet>();
	}
	
	/**
	 * Sort the collection according to the {@link NaturalSet#compareTo(NaturalSet)} order relation.
	 */
	public void sort(){
		Collections.sort(elements);		
	}
	
	/**
	 * Trims the capacity of this collection instance to be the current size. 
	 */
	public void trimToSize(){
		elements.trimToSize();
	}
	
	/**
	 * Add a {@link NaturalSet} to the collection.
	 * @param e
	 * @throws NaturalSetException if <code>this.domain().equals(e.domain())</code> does not evaluate to <code>true</code>.
	 */
	public void add(NaturalSet e) throws NaturalSetException{
		if(domain.equals(e.domain)){
			elements.add(e);
		} else {
			throw new NaturalSetException("Entry domain is incompatible with natural set set's domain.");
		}
	}
	
	/**
	 * @return The number of elements in the collection
	 */
	public int size(){
		return elements.size();
	}

	/**
	 * Computes the size of the intersect of the two collections.
	 * This does not modify either of the collections.
	 * @param other <code>NaturalSetCollection</code> to intersect with
	 * @return The cardinality of the intersect of the two collections.
	 * @throws NaturalSetException if <code>this.domain().equals(other.domain())</code> does not evaluate to <code>true</code>.
	 */
	public int intersectCount(final NaturalSetCollection other) throws NaturalSetException{
		int intersect = 0;
		if(this.domain.equals(other.domain)){
			ArrayList<NaturalSet> first = elements, second = other.elements;
			int firstPosition = 0, secondPosition = 0;
			int firstEnd = 	first.size() - 1, secondEnd = second.size() - 1;
	
			while(firstPosition < firstEnd && secondPosition < secondEnd) {
				int order = first.get(firstPosition).compareTo(second.get(secondPosition));			
				if(order == 0){
					intersect++;
					firstPosition++;
					secondPosition++;
				} else if (order > 0) {
					secondPosition++;
				} else {
					firstPosition++;
				}
			}
			
			if(firstPosition != firstEnd){
				first = other.elements;
				second = elements;
	
				secondPosition = firstPosition;
				firstPosition = secondEnd;
				secondEnd = firstEnd;
				firstEnd = firstPosition;
			} 
			
			NaturalSet last = first.get(firstEnd);
			NaturalSet current = second.get(secondPosition);
			
			int order = current.compareTo(last);	
			
			while(secondPosition <= secondEnd && order <= 0){
				if(order == 0){ intersect++; }
				secondPosition++;
			}
		} else {
			throw new NaturalSetException("NaturalSetCollections have to be defined over the same Domain to intersect");
		}
		return intersect;	
	}
	
	/**
	 * Computes the cardinality of the intersect of the two collections.
	 * This does not modify either of the collections.
	 * @param other <code>NaturalSetCollection</code> to intersect with
	 * @return The cardinality of the intersect of the two collections.
	 * @throws NaturalSetException if <code>this.domain().equals(other.domain())</code> does not evaluate to <code>true</code>.
	 */
	
	/**
	 * Constructs a collection containing the intersection of <code>a</code> and <code>b</code>. 
	 * @param a
	 * @param b
	 * @return A new <code>NaturalSetCollection</code> object
	 * @throws NaturalSetException if <code>a.domain().equals(b.domain())</code> does not evaluate to <code>true</code>.
	 */
	public static NaturalSetCollection intersect(NaturalSetCollection a, NaturalSetCollection b) throws NaturalSetException{
		NaturalSetCollection intersect = null;
		if(b.domain.equals(a.domain)){			
			intersect = new NaturalSetCollection(b.domain);
			ArrayList<NaturalSet> first = b.elements, second = a.elements;
			int firstPosition = 0, secondPosition = 0;
			int firstEnd = 	first.size() - 1, secondEnd = second.size() - 1;
	
			while(firstPosition < firstEnd && secondPosition < secondEnd) {
				int order = first.get(firstPosition).compareTo(second.get(secondPosition));			
				if(order == 0){
					intersect.elements.add(first.get(firstPosition));
					firstPosition++;
					secondPosition++;
				} else if (order > 0) {
					secondPosition++;
				} else {
					firstPosition++;
				}
			}
			
			if(firstPosition != firstEnd){
				first = a.elements;
				second = b.elements;
	
				secondPosition = firstPosition;
				firstPosition = secondEnd;
				secondEnd = firstEnd;
				firstEnd = firstPosition;
			} 
			
			NaturalSet last = first.get(firstEnd);
			NaturalSet current = second.get(secondPosition);
			
			int order = current.compareTo(last);	
			
			while(secondPosition <= secondEnd && order <= 0){
				if(order == 0){ intersect.elements.add(current); }
				secondPosition++;
			}
		} else {
			throw new NaturalSetException("NaturalSetCollections have to be defined over the same Domain to intersect");
		}
		
		return intersect;
	}

	/**
	 * Xor every {@link NaturalSet} in the collection with the corresponding <code>NaturalSet</code> in <code>other</code>.  
	 * @param other The <code>NaturalSetCollection</code> to Xor with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 * @throws NaturalSetException if <code>this.domain().equals(other.domain())</code> does not evaluate to <code>true</code> or the two sets are of different cardinality.
	 */
	public NaturalSetCollection lineByLineXor(NaturalSetCollection other) throws NaturalSetException {
		if(size() == other.size() && domain.equals(other.domain)){
			for (int i=0; i < size(); i++) {
				elements.get(i).xor(other.elements.get(i)); 
			}
			
		} else { 
			throw new NaturalSetException("Line by line XOR is only defined for sets of equal domain and cardinality");
		
		} return this;
	}	
	
	/**
	 * Xor every {@link NaturalSet} with <code>other</code>.  
	 * @param other The <code>NaturalSet</code> to Xor with.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 * @throws NaturalSetException if <code>this.domain().equals(other.domain())</code> does not evaluate to <code>true</code>.
	 */
	public NaturalSetCollection xor(NaturalSet other) throws NaturalSetException{
		if(domain.equals(other.domain)){
			for (int i=0; i < size(); i++) {
				elements.get(i).xor(other); 
			}
			
		} else {
			throw new NaturalSetException("XOR is only defined for sets of equal domain");
			
		} return this;
	}
		
	/**
	 * Print the collection as a binary matrix.
	 * @param w the <code>PrintWriter</code> to print to.
	 * @throws IOException
	 */
	public void print(PrintWriter w) throws IOException{
		for(NaturalSet naturalSet : elements){
			w.println(naturalSet.toBinaryString());
		}
	}	
	
	/**
	 * @param index
	 * @return The <code><NaturalSet/code> at <code>index</code> position.
	 */
	public NaturalSet get(int index){
		return elements.get(index);
	}
	
	/**
	 * Construct the transposed collection.
	 * @param c The collection to transpose.
	 * @return A new <code>NaturalSetCollection</code> object.
	 * @throws NaturalSetException If <code>c</code> is empty.
	 */
	public static NaturalSetCollection transpose(NaturalSetCollection c) throws NaturalSetException {
		NaturalSetCollection t = null;
		if(c.size() > 0){
			t = new NaturalSetCollection(new NaturalDomain(c.size()), c.domain.closureCardinality());
			for(int j=c.domain.min(); j<=c.domain.max(); j++){
				OpenBitSet tline = new OpenBitSet(c.size());
				for(int i=0; i<c.size(); i++){
					if(c.elements.get(i).contains(j)){tline.fastSet(i);}
				}
				t.add(new NaturalSet(t.domain, tline));
			}
		} else {
			throw new NaturalSetException("Can not transpose empty collection.");
		}
		return t;
	}	
	
	public Iterator<NaturalSet> iterator(){
		return new NaturalSetIterator(this);
	}
	
	public class NaturalSetIterator implements Iterator<NaturalSet> {
		private Iterator<NaturalSet> iterator;
		
		public NaturalSetIterator(NaturalSetCollection naturalSetSet){
			iterator = naturalSetSet.elements.iterator();
		}

		public boolean hasNext() {
			return iterator.hasNext();
		}

		public NaturalSet next() {
			return iterator.next();
		}

		public void remove() {
			iterator.remove();
		}
	}

	/**
	 * @return The {@link NaturalDomain} over which the {@link NaturalSet} objects in the collection are defined. 
	 */
	public NaturalDomain domain() {
		return domain;
	}		

	public String toString() {
		StringBuilder display = new StringBuilder();
		display.append("NaturalSetCollection {");
		display.append("\n\tNaturalSet domain: ");
		display.append(domain );
		display.append("\n\tsize: ");
		display.append(size() );
		display.append("\n}");
					
		return display.toString();
	}	
}
