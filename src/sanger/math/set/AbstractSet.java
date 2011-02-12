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

import java.util.Locale;

import org.apache.solr.util.OpenBitSet;


/**
 * Abstract implementation of a finite, set of natural numbers.
 * <p>The set can be stored internally in either the <strong>bounds</strong> or the <strong>binary</strong> mode. 
 * The bounds mode is specified by the <code>min</code> and <code>max</code> values, 
 * which is memory efficient but only supports continuous intervals. The binary mode stores an {@link OpenBitSet}, 
 * consuming more memory for large sets, but capable of encoding non continuous intervals.
 * When continuous, the set will automatically collapse to the bounds mode.
 * {@link SetTopology#EMPTY} sets are defined as a special case of the bounds mode where both <code>min</code> and <code>max</code> are set to {@link Integer#MIN_VALUE}.</p>
 * @author Lior Galanti
 */
public abstract class AbstractSet implements Cloneable{
	
	/**
	 * Enumerate the three valid set states.
	 */
	public enum SetTopology {
		/**
		 * Set contains no elements.
		 */
		EMPTY,
		
		/**
		 * Set contains a continuous set of natural numbers.
		 */
		CLOSED,
		/**
		 * Set contains a, possibly, non continuous set of natural numbers.
		 */
		FRAGMENTS; 
		/**
		 * return the lower case name of the state.
		 */
		public String toString() { return (this.name()).toLowerCase(Locale.ENGLISH); }
	};

	protected OpenBitSet map;
	protected int min;
	protected int max;
	
	/**
	 * Constructs an {@link SetTopology#EMPTY} set.
	 */
	protected AbstractSet(){
		this.map = null;
		this.min = Integer.MIN_VALUE;
		this.max = Integer.MIN_VALUE;
	}
		
	/**
	 * Constructs a {@link SetTopology#CLOSED} set of the interval <code>[min, max]</code>.
	 * @param min Smallest natural number in the set.
	 * @param max Largest natural number in the set.
	 */
	protected AbstractSet(int min, int max){
		this.map = null;
		this.min = min;
		this.max = max;
		if(min > max) { clear(); }
	}
	
	/**
	 * Constructs a set from <code>map</code>.
	 * @param map The {@link OpenBitSet} of natural numbers to be contained in the set.
	 */
	protected AbstractSet(OpenBitSet map){
		this.map = map;
		this.min = Integer.MIN_VALUE;
		this.max = Integer.MIN_VALUE;
	}
	
	/**
	 * Evaluates the {@link SetTopology} of the set.
	 * The value indicates that the following conditions hold:
	 * <ul>
	 * 	<li><code>{@link SetTopology#EMPTY}: min == Integer.MIN_VALUE && max == Integer.MIN_VALUE && map == null</code></li>
	 * 	<li><code>{@link SetTopology#CLOSED}: min > Integer.MIN_VALUE && max > Integer.MIN_VALUE && min >= max && map == null</code></li>
	 * 	<li><code>{@link SetTopology#FRAGMENTS}: min == Integer.MIN_VALUE && max == Integer.MIN_VALUE && map != null</code></li>
	 * </ul>
	 * @return The current topology
	 */
	public SetTopology topology(){		
		SetTopology result = SetTopology.EMPTY;			
		if(min != Integer.MIN_VALUE) result = SetTopology.CLOSED; 
		else if(map != null) result = SetTopology.FRAGMENTS; 
		return result;
	}

	/**
	 * Resets the set to a valid {@link SetTopology#EMPTY} state, regardless of its initial state.
	 */
	protected void clear(){
		min = Integer.MIN_VALUE;
		max = Integer.MIN_VALUE;
		map = null;
	}

	/**
	 * Convenience function for testing if the set is empty.
	 * @return True if the set is {@link SetTopology#EMPTY}.
	 */
	public boolean isEmpty(){
		return (topology() == SetTopology.EMPTY);
	}
	
	/**
	 * Convenience function for testing if the set is continuous.
	 * A set is continuous if it is either {@link SetTopology#CLOSED} or {@link SetTopology#EMPTY}.
	 * @return returns True if the set is continuous.
	 */
	public boolean isContinuous(){
		return (topology() == SetTopology.EMPTY || topology() == SetTopology.CLOSED);
	}

	/**
	 * Make this set contain the same elements as the other set.
	 * @param other The set to copy from.
	 */
	protected void copy(final AbstractSet other){
		this.min = other.min;
		this.max = other.max;
		this.map = other.map == null ? null : (OpenBitSet)other.map.clone(); 
	}
	
	/**
	 * Compute the cardinality of the set's closure.
	 * The closure is the smallest continuous set that contains the set. 
	 * @return Cardinality of the set's closure.
	 */
	public int closureCardinality(){
		return isEmpty() ? 0 : max() - min() + 1;
	}	
	
	/**
	 * Returns the smallest element in the set.
	 * @return The smallest element in the set.
	 */
	public abstract int min();

	/**
	 * Returns the largest element in the set.
	 * @return The largest element in the set.
	 */
	public abstract int max();
	
	/**
	 * Returns The cardinality of the set.
	 * @return The cardinality element in the set.
	 */
	public abstract int cardinality();

	/**
	 * Tests if the given element is contained in the set.
	 * @return True if the given element is contained in the set.
	 */
	public abstract boolean contains(int position);
	
	/**
	 * Changes the set to the binary mode.
	 */
	protected abstract void vectorize();
	
	/**
	 * Changes the set to the bounds mode, if it is continuous.
	 */
	protected abstract void normalize();	
	
	/**
	 * Returns the binary vector for the set, regardless of its mode.
	 * @return The {@link OpenBitSet} of the elements contained in the set.
	 */
	protected abstract OpenBitSet vector();
}
