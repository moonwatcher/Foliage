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

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.solr.util.OpenBitSet;

import sanger.argml.model.NonReoccuringDeque;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetException;
import sanger.math.set.NaturalSpace;

/**
 * An implementation of an acestral recombination graph vertex.
 * Vertices are multi furcated and have an implicit active region 
 * derived from the active regions of the edges that go through them.
 * The sum of all active regions on edges coming in and of those going out 
 * should be identical and is, in fact, the flux going through the vertex.
 * @author Lior Galanti
 */
public class Vertex implements Comparable<Vertex>, Iterable<Edge>{	
	private Genealogy genealogy;
	private int id;
	protected LinkedList<Edge> edges = new LinkedList<Edge>();	

	/**
	 * Constructs a new Vertex in <code>genealogy</code>.
	 * @param genealogy to add the new vertex to.
	 * @param id Id for the new vertex.
	 */
	protected Vertex(Genealogy genealogy, int id) {
		this.genealogy = genealogy;
		this.id = id;
		this.genealogy.vertices.add(this);
	}
	
	/**
	 * Constructs a trasitional vertex.
	 * @param genealogy {@link Genealogy} to relate the vertex to.
	 * @param other Vertex to copy.
	 * @throws NaturalSetException if the target genealogy is defined over a different {@link NaturalSpace} than the source.
	 */
	protected Vertex(Genealogy genealogy, Vertex other) throws NaturalSetException{
		this.genealogy = genealogy;
		this.id = other.id;
		this.genealogy.vertices.add(this);
		for(Edge e : other){			
			if(other.isSource(e)){
				NaturalSet ar = genealogy.snpDomain().project(e.activeRegion());
				if(!ar.isEmpty()){
					new Edge(e, this, ar);
				}
			} 
		}
	}
	
	/**
	 * Creates a new vertex and transferes all activity on positions bigger than <code>position</code> to it.
	 * New vertex will have copies of all the outgoing edges with activity in positions bigger than <code>position</code>.
	 * Edges coming out of the original vertex will be active only on positions smaller or equal to <code>position</code>.
	 * Mutations are devided between the two so that they remain on an edge active in their position.
	 * @param id Id for the new vertex.
	 * @param position Position to split the active regions.
	 * @return The new vertex.
	 */
	public Vertex split(int id, int position){
		Vertex v = null;
		v = new Vertex(genealogy, id);
		for ( Edge edge : this ) {
			if ( isSource(edge) ) {				
				NaturalSet region = edge.activeRegion.clone();
				region.removeBelow(position + 1);
				Edge e = new Edge(v, edge.target, region);					
				Iterator<Mutation> iterate = edge.iterator();
				while(iterate.hasNext()){
					Mutation m = iterate.next();
					if(region.contains(m.position())){
						iterate.remove();
						e.mutations.add(m);
					}
				}
				edge.activeRegion.removeAbove(position);
			}
		}
		return v;	
	}
	
	/**
	 * Merge <code>other</code> vertex into this one.
	 * All edges pointing at or from <code>other</code> are redirected to this vertex, 
	 * this vertex's id is set to <code>id</code> and <code>other</code> is removed from the genealogy. 
	 * @param id New id for this vertex after the merge is complete.
	 * @param other Vertex to merge into this vertex.
	 * @return <code><strong>this</strong></code> to facilitate chaining.
	 */
	public Vertex merge(int id, Vertex other) {
		for(Edge edge : other ){			
			if (other.isSource(edge)) edge.source = this;			
			if (other.isTarget(edge)) edge.target = this;
			this.addEdge(edge);							
		}
		
		other.edges.clear();
		genealogy.vertices.remove(other);		
		this.id = id;
		return this;
	}
	
	/**
	 * Construct a new vertex and connect it to this vertex as an ancestor.
	 * A new {@link Edge} is also created with the new vertex as a source and this vertex as target.
	 * The new edge will have the full active region of this vertex.
	 * @param id Id for the new Vertex.
	 * @return The new vertex.
	 */
	public Edge createAncestor(int id){
		return createAncestor(id, outActiveRegion());
	}

	/**
	 * Construct a new vertex and connect it to this vertex as an ancestor.
	 * A new {@link Edge} is also created with the new vertex as a source and this vertex as target.
	 * The new edge will have a <code>region</code> active region.
	 * @param id Id for the new Vertex.
	 * @param region active region for the edge from the ancestor.
	 * @return The new vertex.
	 */
	public Edge createAncestor(int id, NaturalSet region) {
		Vertex ancestor = new Vertex(genealogy, id);		
		return new Edge(ancestor, this, region);
	}
		
	/**
	 * Connect <code>ancestor</code> to this vertex.
	 * The new edge will point at this vertex as a target and at <code>ancestor</code> as a source.
	 * The new edge will have the full active region of this vertex. 
	 * @param ancestor Vertex to connect to.
	 */
	public void connectToAncestor(Vertex ancestor){
		new Edge(ancestor, this, outActiveRegion());	
	}

	/**
	 * Add an edge to the vertex.
	 * No validation preformed.
	 * @param edge
	 */
	protected void addEdge(Edge edge) {
		edges.add(edge);
	}
	
	
	/**
	 * Defined as the natural order relation on the vertex id. 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Vertex other){
		return this.id - other.id;
	}

	/**
	 * Two vertices are concidered equal if they have the same id.
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	public boolean equals(Object o) {
		boolean result = false;
		if (o instanceof Vertex) {
			Vertex other = (Vertex) o;
			if(other.id == this.id) { result = true; }
		}
		return result;
	}
	
	/**
	 * The vertex has only 1 incoming and 1 outgoing {@link Edge} active on <code>region</code>.
	 * @param region
	 * @return True if <code>(inDegree(region) == 1 && outDegree(region) == 1)</code> evaluates to <code>true</code>.
	 * @throws NaturalSetException if <code>region</code> is not defined over the same {@link NaturalDomain} as the host {@link Genealogy}.
	 */
	public boolean isDegenerate(NaturalSet region) throws NaturalSetException{
		return (inDegree(region) == 1 && outDegree(region) == 1);
	}
	
	/**
	 * The vertex is the grand most recent common ancestor.
	 * @return True if the vertex is the gmrca.
	 */
	public boolean isGmrca() {
		return (inDegree() == 0);
	}
		
	/**
	 * The vertex is a leaf node.
	 * The leaf nodes are the sequences in the present time.
	 * @return True if the vertex is a leaf node.
	 */
	public boolean isLeaf() {
		return (outDegree() == 0);
	}
		
	/**
	 * This vertex is <code>edge</code>'s target.
	 * @param edge {@link Edge} to test.
	 * @return True if this vertex is <code>edge</code>'s target.
	 */
	public boolean isTarget(Edge edge){
		return (edge.target() == this);
	}
	
	/**
	 * This vertex is <code>edge</code>'s source.
	 * @param edge {@link Edge} to test.
	 * @return True if this vertex is <code>edge</code>'s source.
	 */
	public boolean isSource(Edge edge){
		return (edge.source() == this);
	}
	
	/**
	 * Count the number of incoming {@link Edge} objects.
	 * @return The number of incoming {@link Edge} objects.
	 */
	public int inDegree(){
		int result = 0;
		for ( Edge edge : edges ) {
			if ( isTarget(edge) ) { result++; }
		} return result;
	}
	
	/**
	 * Count the number of incoming {@link Edge} objects on <code>region</code>.
	 * @param region Active region on which edges are counted.
	 * @return The number of incoming {@link Edge} objects on <code>region</code>.
	 * @throws NaturalSetException if <code>region</code> is not defined over the same {@link NaturalDomain} as the host {@link Genealogy}.
	 */
	public int inDegree(NaturalSet region) throws NaturalSetException{
		int result = 0;
		for ( Edge edge : edges ) {
			if ( isTarget(edge) && edge.activeRegion().intersectionCount(region) > 0) { result++; }
		} return result;
	}
	
	/**
	 * Count the number of outgoing {@link Edge} objects.
	 * @return The number of outgoing {@link Edge} objects.
	 */
	public int outDegree(){
		int result = 0;
		for ( Edge edge : edges ) {
			if ( isSource(edge) ) { result++; }
		} return result;
	}

	/**
	 * Count the number of outgoing {@link Edge} objects on <code>region</code>.
	 * @param region Active region on which edges are counted.
	 * @return The number of outgoing {@link Edge} objects on <code>region</code>.
	 * @throws NaturalSetException if <code>region</code> is not defined over the same {@link NaturalDomain} as the host {@link Genealogy}.
	 */	
	public int outDegree(NaturalSet region) throws NaturalSetException{
		int result = 0;
		for ( Edge edge : edges ) {
			if ( isSource(edge) && edge.activeRegion().intersectionCount(region) > 0) { result++; }
		} return result;
	}

	/**
	 * Union the active regions on all incoming edges.
	 * This is defined as the vertex's incoming active region
	 * and under valid conditions should be equal to {@link #outActiveRegion()}. 
	 * The method returns a new instance so changes to it will not effect the vertex.
	 * @return The union of active regions of all the vertex's incoming edges.
	 */
	public NaturalSet inActiveRegion(){
		NaturalSet region;
		if(inDegree() > 0){
			region = genealogy.snpDomain.createEmptyNaturalSet();
			for ( Edge edge : edges ) {
				try { if (isTarget(edge)) region.union(edge.activeRegion());} 
				catch (NaturalSetException e) {genealogy.env().log().printError(e);}
			}			
		} else {
			region = genealogy.snpDomain.createCompleteNaturalSet();
		} return region;
	}
	
	/**
	 * Union the active regions on all outgoing edges.
	 * This is defined as the vertex's active region
	 * and under valid conditions should be equal to {@link #inActiveRegion()}. 
	 * The method returns a new instance so changes to it will not effect the vertex.
	 * @return The union of active regions of all the vertex's outgoing edges.
	 */
	public NaturalSet outActiveRegion(){
		NaturalSet region;
		if(outDegree() > 0){
			region =  genealogy.snpDomain.createEmptyNaturalSet();
			for ( Edge edge : edges ) {
				try { if (isSource(edge)) region.union(edge.activeRegion());} 
				catch (NaturalSetException e) {genealogy.env().log().printError(e);}
			} 
		} else {
			region = genealogy.snpDomain.createCompleteNaturalSet();
		} return region;
	}

	/**
	 * Trace all mutation events on the lineage leading to this vertex from the {@link Genealogy#gmrca() GMRCA}.
	 * The returned value indicates which snp locations have mutated since the GMRCA.
	 * For the actual sequence at this vertex see {@link #haplotype(NaturalSet)}.
	 * @param region Active region on which to test for possible mutation events.
	 * @return All snp positions that had mutation events.
	 * @throws NaturalSetException if <code>region</code> is not defined over the same {@link NaturalDomain} as the host {@link Genealogy}.
	 */
	public NaturalSet mutations(NaturalSet region) throws NaturalSetException {
		HashMap<Vertex, NaturalSet> activity = new HashMap<Vertex, NaturalSet>();
		Deque<Vertex> tovisit = new ArrayDeque<Vertex>();
		NonReoccuringDeque<Vertex> verticesToVisit = new NonReoccuringDeque<Vertex>();

		activity.put(this, NaturalSet.intersect(inActiveRegion(), region));
		tovisit.add(this);
		
		Vertex current;
		NaturalSet currentActivity;
		
		while(!tovisit.isEmpty()){
			current = tovisit.poll();
			currentActivity = activity.get(current);
			
			for ( Edge edge : current ) {			
				NaturalSet incomingActivity = NaturalSet.intersect(currentActivity, edge.activeRegion());
				
				if ( current.isTarget(edge) && !incomingActivity.isEmpty()) {
					NaturalSet sourceActivity = activity.get(edge.source());	

					if ( sourceActivity == null ) {	// First visit
						activity.put(edge.source(), incomingActivity);
						tovisit.add(edge.source());

					} else {	// Consecutive visits	
						// Extend activity to include activity from new route.
						int compare = sourceActivity.unionAndCompare(incomingActivity);

						// Reschedule the source node for scanning
						if (!tovisit.contains(edge.source()) && compare > 0) {
							tovisit.add(edge.source());
						}
						
					}
				}
			}
		}			
		
		
		OpenBitSet result = new OpenBitSet(genealogy.snpDomain().closureCardinality());
		verticesToVisit.push(this);
		
		while (verticesToVisit.size() > 0) {
			current = verticesToVisit.poll();
			for (Edge edge : current) {
				activity.containsKey(current);
				NaturalSet active = NaturalSet.intersect(edge.activeRegion(), activity.get(current));
				if (current.isTarget(edge) && !active.isEmpty()) {
					for(Mutation mutation : edge){
						if(active.contains(mutation.position())){
							result.fastSet(genealogy.snpDomain().toRelativeCoordinate(mutation.position()));
						}
					}
					verticesToVisit.add(edge.source());
				}
			}
		}
		return new NaturalSet(genealogy.snpDomain, result);
	}
	
	/**
	 * Trace all mutation events on the lineage leading to this vertex from the {@link Genealogy#gmrca() GMRCA}.
	 * The returned value indicates which snp locations have mutated since the GMRCA.
	 * For the actual sequence at this vertex see {@link #haplotype()}.
	 * @return All snp positions that had mutation events.
	 * @throws NaturalSetException if <code>region</code> is not defined over the same {@link NaturalDomain} as the host {@link Genealogy}.
	 */
	public NaturalSet mutations() throws NaturalSetException {
		return mutations(genealogy.snpDomain.createCompleteNaturalSet());
	}
		
	/**
	 * The haplotype at this vertex.
	 * @param region Region to get the haplotype for.
	 * @return The haplotype.
	 * @throws NaturalSetException if <code>region</code> is not defined over the same {@link NaturalDomain} as the host {@link Genealogy}.
	 */
	public NaturalSet haplotype(NaturalSet region) throws NaturalSetException {
		NaturalSet hap = mutations(region);
		hap.xor(genealogy.ancestralSequence);
		return hap;
	}
	
	/**
	 * The haplotype at this vertex.
	 * @return The haplotype.
	 * @throws NaturalSetException if <code>region</code> is not defined over the same {@link NaturalDomain} as the host {@link Genealogy}.
	 */
	public NaturalSet haplotype() throws NaturalSetException {
		return haplotype(genealogy.snpDomain.createCompleteNaturalSet());
	}

	public String toString() {
		StringBuilder display = new StringBuilder();
		display.append("Vertex ");
		display.append(id);
		display.append(" {\n\tOut(");
		display.append(outDegree());
		display.append(")");
		display.append(outActiveRegion());
		display.append(":\n");

		for(Edge edge : edges) {
			if ( isSource(edge) ) {
				display.append("\t\t");
				display.append(edge);
				display.append(" ");
			}
		}

		display.append("\n\tIn(");
		display.append(inDegree());
		display.append(")");
		display.append(inActiveRegion());
		display.append(":\n");

		for(Edge edge : edges) {
			if (isTarget(edge) ) {
				display.append("\t\t");
				display.append(edge);
				display.append("\n");
			}
		}
		display.append("\n}"); 
					
		return display.toString();
	}
	
	public Iterator<Edge> iterator(){
		return edges.iterator();
	}

	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	/**
	 * @param id the id to set
	 */
	public void setId(int id) {
		this.id = id;
	}	

	public Genealogy genealogy(){
		return genealogy;
	}
	
}
