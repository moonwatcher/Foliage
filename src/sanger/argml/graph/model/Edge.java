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
import java.util.LinkedList;

import sanger.math.set.NaturalSet;

/**
 * An ARG edge, directed forward in time. Edges have an active region 
 * which contains all the positions in the sequence with ancestral material
 * and a length. {@link Mutation Mutations} can also occure on the edge.
 * @author Lior Galanti
 */
public class Edge implements Comparable<Edge>, Iterable<Mutation>{
	
	protected Vertex source;
	protected Vertex target;
	protected double length;
	protected NaturalSet activeRegion;
	protected LinkedList<Mutation> mutations;

	/**
	 * Constructs a new transitional edge.
	 * This is part of the procedure for cropping or copying a {@link Genealogy genealogy}.
	 * Copying is preformed forward in time and creates an edge from 
	 * the <code>source</code> {@link Vertex vertex} in the new genealogy 
	 * to the <code>target</code> vertex in the original.
	 * The new edge's active region contains only coordinates and {@link Mutation mutations} 
	 * also contained in <code>region</code>.
	 * @param edge Edge in the original genealogy.
	 * @param source Source vertex in the original genealogy.
	 * @param region Region to copy when constructing the the new edge.
	 */
	protected Edge(Edge edge, Vertex source, NaturalSet region) {
		this.source = source;
		this.target = edge.target;
		this.mutations = new LinkedList<Mutation>();
		this.activeRegion = region;
		this.length = 1.0;
		
		source.addEdge(this);
		source.genealogy().edges.add(this);
		for(Mutation mutation : edge.mutations){
			if(region.contains(mutation.position())) {
				mutations.add(new Mutation(source.genealogy(), mutation.position()));
			}
		}
	}
	
	/**
	 * Constructs a new edge connecting <code>source</code> and <code>target</code> vertices.
	 * @param source Source {@link Vertex}.
	 * @param target Target {@link Vertex}. 
	 * @param region Active region of the edge.
	 */
	protected Edge(Vertex source, Vertex target, NaturalSet region) {
		this.source = source;
		this.target = target;
		this.activeRegion = region;		
		this.mutations = new LinkedList<Mutation>();
		this.length = 1.0;
		
		source.addEdge(this);		
		target.addEdge(this);
		source.genealogy().edges.add(this);
	}
	
	
	/**
	 * Mutate snp in given position while traversing the edge. 
	 * @param position
	 */
	public void mutate(int position){
		this.mutations.add(new Mutation(source.genealogy(),position));
	}
		
	/**
	 * @return True if there are mutation events on the edge.
	 */
	public boolean hasMutations(){
		return mutations.isEmpty() ? false : true;
	}

	/**
	 * @param region
	 * @return True if there are mutations in <code>region</code> on the edge.
	 */
	public boolean hasMutations(NaturalSet region){
		boolean result = false;
		for(Mutation mutation : mutations){
			if(region.contains(mutation.position())){
				result = true;
				break;
			}
		} 			
		return result;
	}
			
	public String toString(){
		StringBuilder display = new StringBuilder();
		display.append("Edge: ");
		display.append(source.getId());
		display.append(" > ");
		display.append(target.getId());
		display.append(" active on: ");
		display.append(activeRegion);
		if (!mutations.isEmpty()) { 
			display.append(" mutations: ");
			display.append(mutations);
		}
		return display.toString();
	}

	public Iterator<Mutation> iterator(){
		return mutations.iterator();
	}

	
	public int compareTo(Edge other){
		return activeRegion.compareTo(other.activeRegion);
	}
	
	/**
	 * @return the activeRegion
	 */
	public NaturalSet activeRegion() {
		return activeRegion;
	}
	

	/**
	 * @return the source
	 */
	public Vertex source() {
		return source;
	}
	

	/**
	 * @return the target
	 */
	public Vertex target() {
		return target;
	}	
	
	/**
	 * @return the length
	 */
	public double length() {
		return length;
	}
}
