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

/**
 * Infinit sites mutation event. Only biallelic sequences are supported, 
 * so a nucleotide can have either the ancestral or the mutated value.
 * Mutation events occure on the ancestral recombination graph and stored on {@link Edge edges}.
 * @author Lior Galanti
 *
 */
public class Mutation implements Comparable<Mutation> {	
	protected Genealogy genealogy;
	private int position;
	
	/**
	 * Constructs a new mutation event
	 * @param genealogy
	 * @param position
	 */
	protected Mutation(Genealogy genealogy, int position){
		this.genealogy = genealogy;
		this.position = position;
	}
	
	public int compareTo(Mutation mutation){
		return this.position - mutation.position;
	}
	
	/**
	 * @return The position where the <code>Mutation</code> occurs.
	 */
	public int position() {
		return position;
	}

	public String toString() { return String.valueOf(position); }
}
