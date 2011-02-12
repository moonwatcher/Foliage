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
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedList;

import org.apache.solr.util.BitSetIterator;
import org.apache.solr.util.OpenBitSet;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.argml.io.TextOutput;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetCollection;
import sanger.math.set.NaturalSetException;

public class Genealogy extends Environmental implements Cloneable{	
	protected NaturalDomain snpDomain;
	protected NaturalDomain haplotypeDomain;
	protected NaturalSet ancestralSequence;	

	protected LinkedList<Vertex> vertices;
	protected LinkedList<Edge> edges;
	private LinkedList<Vertex> leaves;	
	protected Vertex gmrca;

	/**
	 * Construct an empty genealogy and populates a vertex for every haplotype.
	 * The ancestral sequence is set to 0's.
	 * @param env The runtime enviroment
	 * @param snpDomain The snp domain of the haplotypes in the sample.
	 * @param haplotypeDomain The haplotype domain of the sample.
	 */
	public Genealogy(Environment env, NaturalDomain snpDomain, NaturalDomain haplotypeDomain) {
		super(env);
		this.vertices = new LinkedList<Vertex>();
		this.edges = new LinkedList<Edge>();
		this.leaves = new LinkedList<Vertex>();
		this.gmrca = null;
		this.snpDomain = snpDomain; 
		this.haplotypeDomain = haplotypeDomain; 
		this.ancestralSequence = snpDomain.createEmptyNaturalSet();
		for(int i=haplotypeDomain.min(); i<=haplotypeDomain.max(); i++) {
			if(haplotypeDomain.contains(i)){
				leaves.add(new Vertex(this, i));
			}
		}
	}
	
	/**
	 * Construct an empty genealogy.
	 * @param env The runtime enviroment
	 * @param snpDomain The snp domain of the haplotypes in the sample.
	 * @param haplotypeDomain The haplotype domain of the sample.
	 * @param ancestral Ancestral sequence on the grand most recent common ancestor.
	 */
	protected Genealogy(Environment env, NaturalDomain snpDomain, NaturalDomain haplotypeDomain, NaturalSet ancestral) {
		super(env);
		this.vertices = new LinkedList<Vertex>();
		this.edges = new LinkedList<Edge>();
		this.leaves = null;
		this.gmrca = null;
		this.ancestralSequence = ancestral;
		this.snpDomain = snpDomain;
		this.haplotypeDomain = haplotypeDomain;
	}		

	
	/**
	 * Calculate all possible bipartitions for the tree.
	 * This procedure applies only when {@link #isTree()} return true.
	 * @return A collection of all possible bipartitions.
	 */
	public NaturalSetCollection biPartitions() {
		HashMap<Edge, NaturalSet> edgeToBiPartitionMap = null;
		NaturalSetCollection result = null;
		//	If L denotes the number of haplotypes and the root is multifurcated 
		//	then there are no identical bipartitions but there will be strictly less then L-2.
		//	If the root is bifurcated the L-2 edge in DFS order has to be the second edge coming into the root
		//	which generates the same bipartition as the other one coming into the root.
		//	Either way we always want strictly less then L-2.
		int pmax = 2 * haplotypeDomain.cardinality() - 3;
		
		if(isTree()){			
			edgeToBiPartitionMap = new HashMap<Edge, NaturalSet>();
			Vertex current = gmrca;
			
			while(current != null){
				Edge next = null, back = null;			
				for(Edge edge : current){
					if(current.isSource(edge) && !edgeToBiPartitionMap.containsKey(edge)){
						next = edge;
						break;
					} else {
						if(current.isTarget(edge)){
							back = edge;
						}
					}
				}
				
				if(!(current.outDegree() == 0)){
					if(next != null){
						current = next.target();
					} else {
						if(back != null && edgeToBiPartitionMap.size() < pmax){
							NaturalSet biPartition = haplotypeDomain.createEmptyNaturalSet();					
							for(Edge edge : current){
								if(current.isSource(edge)){
									try { biPartition.union(edgeToBiPartitionMap.get(edge)); } 
									catch (NaturalSetException e) { env().log().printError(e); }
									
								}
							}
							edgeToBiPartitionMap.put(back, biPartition);
							current = back.source();
						} else {
							current = null;
						}
					}
				} else {
					edgeToBiPartitionMap.put(back, new NaturalSet(haplotypeDomain, current.getId(), current.getId()));
					current = back.source();				
				}
			}			
			
			result = new NaturalSetCollection(haplotypeDomain);
			int inf = 1, sup = haplotypeDomain.cardinality() - 1;
			for(NaturalSet p : edgeToBiPartitionMap.values()){
				int card = p.cardinality();
				if( card > inf && card < sup){
					if(p.contains(p.domain().min())) p.inverse();
					try { result.add(p);} 
					catch (NaturalSetException e) { env().log().printError(e); }
				}
			}			
			result.sort();

			result.trimToSize();
		}
		return result;
	}
	
	/**
	 * Calculate all possible bipartitions for all local trees in the graph.
	 * The map stores a collection of bipartions for the local tree for every {@link #recombinationFreeRegions() recombination free regions} in the snp domain. 
	 * @return A map of collections of bipartitions.  
	 */
	public HashMap<NaturalSet, NaturalSetCollection> localTreesBiPartitions() {
		NaturalSetCollection frames = recombinationFreeRegions();
		HashMap<NaturalSet, NaturalSetCollection> localTreesBiPartitions = new HashMap<NaturalSet, NaturalSetCollection>(frames.size());

		for(NaturalSet frame : frames){
			if(frame.cardinality() > 0){
				try {			
					NaturalSetCollection biPartitions =  clipSnpDomain(this, frame).biPartitions();
					localTreesBiPartitions.put(frame, biPartitions);
					
				} catch (NaturalSetException e) { env().log().printError(e); }	
			}
		}
		return localTreesBiPartitions;
	}
	
	/**
	 * Count the recombinations on the graph for every snp.
	 * The array corresponds to the {@link #snpDomain() snp domain's} relative coordinates,
	 * use {@link NaturalDomain#toAbsoluteCoordinate(int)} to get absolute position in snp domain.
	 * A value in cell i counts the recombinations that occur between the position of snp i and i+1.   
	 * @return An array of integers with recombination count for every snp in the domain.
	 */
	public int[] recombinationRates(){
		int[] distribution = new int[snpDomain.closureCardinality()];
		for(Vertex vertex : vertices){
			if(vertex.inDegree() > 1){
				LinkedList<Edge> fragments = new LinkedList<Edge>();
				for(Edge edge : vertex){
					if(vertex.isTarget(edge)){
						if(!edge.activeRegion().isEmpty()){ // Do no concider empty regions
							fragments.add(edge);
						}
					}
				}
				Collections.sort(fragments);
				fragments.remove(fragments.size() - 1);
				for(Edge edge : fragments){
					try {
						int position = snpDomain.toRelativeCoordinate(edge.activeRegion().max());
						distribution[position]++;
					} catch (NaturalSetException e) { env().log().printError(e); }	
				}					
			}
		}
		return  distribution;		
	}
	
	/**
	 * Map locations where recombination events occur.
	 * Use {@link #recombinationRates()} for actual recombination counts.
	 * @return A set of all positions where recombination events occur.
	 */
	public NaturalSet recombinationPositions() {
		OpenBitSet recombinations = new OpenBitSet(snpDomain.closureCardinality());
		for(Vertex vertex : vertices){
			if(vertex.inDegree() > 1){
				for(Edge edge : vertex){
					if(vertex.isTarget(edge)){
						if(!edge.activeRegion().isEmpty()){ // Do no concider empty regions
							try {recombinations.fastSet(snpDomain.toRelativeCoordinate(edge.activeRegion().max()));} 
							catch (NaturalSetException e) { env().log().printError(e); }							
						}
					}
				}
			}
		}
		return  new NaturalSet(snpDomain, recombinations);
	}
		
	/**
	 * A collection of non overlapping sets where no recombination events occur.
	 * @return A collection of sets from the snp domain.
	 */
	public NaturalSetCollection recombinationFreeRegions(){
		NaturalSetCollection regions = new NaturalSetCollection(snpDomain);						
		BitSetIterator iterate = new BitSetIterator(recombinationPositions().toOpenBitSet());			
		int left = 0;
		int right = iterate.next();
		while(right > -1){
			try {regions.add(new NaturalSet(snpDomain, snpDomain.toAbsoluteCoordinate(left), snpDomain.toAbsoluteCoordinate(right)));} 
			catch (NaturalSetException e) { env().log().printError(e); }
			
			left = right + 1;
			right = iterate.next();
		}
		regions.trimToSize();
		regions.sort();

		return regions;
	}
	
	/**
	 * Clip the SNP domain.
	 * Constructs a new genealogy spanning only the positions contained in <code>region</code>.
	 * @param genealogy Genealogy to clip <code>region</code> from.
	 * @param region Set of position to crop.
	 * @return A new genealogy containing the clip.
	 * @throws NaturalSetException if <code>region</code> is not from the genealogy's snp domain.
	 */
	public static Genealogy clipSnpDomain(final Genealogy genealogy, final NaturalSet region) throws NaturalSetException {
		return clipSnpDomain(genealogy, region, genealogy.gmrca);
	}
		
	/**
	 * Clip the SNP domain from <code>vertex</code> to the present.
	 * Constructs a new genealogy from the sub graph rooted at <code>vertex</code>, 
	 * spanning only the positions contained in <code>region</code>.
	 * @param genealogy Genealogy to clip <code>region</code> from.
	 * @param region Set of position to crop.
	 * @param vertex Vertex to start clipping from.
	 * @return A new genealogy containing the clip.
	 * @throws NaturalSetException if <code>region</code> is not from the genealogy's snp domain.
	 */
	public static Genealogy clipSnpDomain(final Genealogy genealogy, final NaturalSet region, final Vertex vertex) throws NaturalSetException {
		Genealogy clip = null;
		if(genealogy.snpDomain.equals(region.domain())){
			NaturalDomain subDomain = genealogy.snpDomain.createSubDomain(region);
			HashMap<Vertex,Vertex> vertexMap = new HashMap<Vertex, Vertex>(subDomain.closureCardinality());
			Deque<Edge> toVisit = new ArrayDeque<Edge>(subDomain.closureCardinality());		
			clip = new Genealogy(genealogy.env(), subDomain, genealogy.haplotypeDomain, subDomain.project(genealogy.ancestralSequence));
			Vertex vertexToCopy = vertex;
			
			while(vertexToCopy.outDegree(region) == 1){
				for(Edge edge : vertexToCopy){
					if(vertexToCopy.isSource(edge) && edge.activeRegion().intersectionCount(region) > 0){
						vertexToCopy = edge.target();
						break;
					}
				}			
			}
			
			Vertex vertexCopy = new Vertex(clip, vertexToCopy);
			clip.gmrca = vertexCopy;		
			vertexMap.put(vertexToCopy, vertexCopy);
			toVisit.addAll(vertexCopy.edges);
			
			
			while(!toVisit.isEmpty()){
				Edge edgeToNextVertex = toVisit.poll();
				vertexToCopy = edgeToNextVertex.target();
				vertexCopy = vertexMap.get(vertexToCopy);
				LinkedList<Mutation> mutations = new LinkedList<Mutation>();
				
				if(vertexCopy == null) {
					if(!genealogy.env().flag("NoCollapse")) {
						while(vertexToCopy.isDegenerate(region)){
							for(Edge edge : vertexToCopy){
								if(vertexToCopy.isSource(edge) && (edge.activeRegion().intersectionCount(region) > 0)){
									for(Mutation mutation : edge){
										if(region.contains(mutation.position())) {
											mutations.add(new Mutation(clip, mutation.position()));
										}
									}
									vertexToCopy = edge.target();
									break;
								}
							}						
						}					
						vertexCopy = vertexMap.get(vertexToCopy);
					}
					
					if(vertexCopy == null){ 
						vertexCopy = new Vertex(clip, vertexToCopy); 
						vertexMap.put(vertexToCopy, vertexCopy);
						for(Edge edge : vertexCopy){
							if(vertexCopy.isSource(edge)){ toVisit.add(edge); }
						}					
					}
				}
				
				edgeToNextVertex.target = vertexCopy;
				edgeToNextVertex.mutations.addAll(mutations);
				vertexCopy.addEdge(edgeToNextVertex);
			}
			clip.sort();
			
			if(!genealogy.gmrca.equals(vertex)){
				OpenBitSet map = new OpenBitSet(genealogy.haplotypeDomain.closureCardinality());
				for(Vertex leaf : clip.leaves()){
					map.fastSet(leaf.getId());
				}
				clip.haplotypeDomain = genealogy.haplotypeDomain.createSubDomain(new NaturalSet(genealogy.haplotypeDomain, map));
			}
			
		} else {
			throw new NaturalSetException("Clipping region must be of the Genealogy's domain");
		
		}
		
		return clip;
	}
	
	public Genealogy clone(){
		Genealogy clone = null;
		try { 
			clone = Genealogy.clipSnpDomain(this, this.snpDomain().createCompleteNaturalSet());
			
		} catch (Exception e) {
			env().log().printError(e);
		}
		return clone;
	}
	
	/**
	 * Reconstruct the haplotypes on all the leaves by tracing mutation events on the lineage.
	 * @return A collection of haplotypes for all leaves.
	 */	
	public NaturalSetCollection haplotypes(){
		NaturalSetCollection haplotypes = null;
		try { haplotypes = haplotypes(snpDomain.createCompleteNaturalSet());}
		catch (NaturalSetException e) { env().log().printError(e);}
		return haplotypes;
	}
		
	/**
	 * Reconstruct the haplotypes in <code>region</code> on all the leaves by tracing mutation events on the lineage.
	 * @param region Set of SNPs to trace mutations for.
	 * @return A collection of haplotypes for all leaves.
	 * @throws NaturalSetException if region is not defined over the genealogy's domain.
	 */
	public NaturalSetCollection haplotypes(NaturalSet region) throws NaturalSetException {
		NaturalSetCollection hapmat = mutations(region);
		hapmat.xor(ancestralSequence);
		return hapmat;
	}
	
	/**
	 * Trace all mutation events on the lineage leading to all the leaves. 
	 * The returned value contains SNP locations that have mutated since the GMRCA.
	 * @return A collection of mutation sets for all leaves.
	 */
	public NaturalSetCollection mutations(){
		NaturalSetCollection mutations = null;
		try { mutations = mutations(snpDomain.createCompleteNaturalSet());}
		catch (NaturalSetException e) { env().log().printError(e);}
		return mutations;
	}
	
	/**
	 * Trace all mutation events in <code>region</code> on the lineage leading to all the leaves. 
	 * The returned value contains SNP locations that have mutated since the GMRCA. 
	 * @param region Set of SNPs to trace mutations for.
	 * @return A collection of mutation sets for all leaves.
	 * @throws NaturalSetException if region is not defined over the genealogy's domain.
	 */
	public NaturalSetCollection mutations(NaturalSet region) throws NaturalSetException {
		NaturalSetCollection haplotypes = new NaturalSetCollection(snpDomain, haplotypeDomain.cardinality());
		for (Vertex vertex : leaves()){
			haplotypes.add(vertex.mutations(region));
		}
		return haplotypes;
	}
		

	/**
	 * Sort the vertex and edge collections according to their natural order.
	 */
	public void sort(){
		Collections.sort(vertices);
		Collections.sort(edges);
		
		for(Edge edge : edges){
			Collections.sort(edge.mutations);
		}
		
		for(Vertex vertex : vertices){
			Collections.sort(vertex.edges);
		}
	}
	
	public String toString() {
		StringBuilder display = new StringBuilder();
		display.append("Genealogy {");
		display.append(" SNP:");
		display.append(snpDomain );
		display.append(", SEQ:");
		display.append(haplotypeDomain );					
		display.append(", V:");
		display.append(vertices.size());					
		display.append(", E:");
		display.append(edges.size());
		display.append(" }");					
		return display.toString();
	}


	public NaturalDomain snpDomain() {
		return snpDomain;
	}

	public NaturalDomain haplotypeDomain() {
		return haplotypeDomain;
	}	


	/**
	 * @return True if the genealogy has no recombination events.
	 */
	public boolean isTree(){
		return (vertices.size() == (edges.size() + 1));
	}

	/**
	 * Return a list of vertices with zero {@link Vertex#outDegree() out degree}.
	 * @return the leaf vertices
	 */
	public LinkedList<Vertex> leaves(){
		if(leaves == null) {
			leaves = new LinkedList<Vertex>();
			for(Vertex v : vertices){
				if(v.isLeaf()){ leaves.add(v); }
			}
			Collections.sort(leaves);
		}
		return leaves;
	}
	
	/**
	 * @return the vertices
	 */
	public LinkedList<Vertex> vertices() {
		return vertices;
	}

	/**
	 * @return the edges
	 */
	public LinkedList<Edge> edges() {
		return edges;
	}
	
	/**
	 * @return the ancestralSequence
	 */
	public NaturalSet ancestralSequence() {
		return ancestralSequence;
	}
	
	/**
	 * @param ancestralSequence the ancestralSequence to set
	 */
	public void setAncestralSequence(NaturalSet ancestralSequence) {
		this.ancestralSequence = ancestralSequence;
	}	

	/**
	 * @return the grand most recent common ancestor
	 */
	public Vertex gmrca() {
		return gmrca;
	}

	/**
	 * Report edges with empty active region to enviroment log.
	 */
	public void reportEmptyRegions(){
		int count = 0;
		for(Vertex vertex : vertices()){
			for(Edge edge : vertex){
				if(edge.activeRegion().cardinality() == 0){
					count++;
					env().log().printError("empty active region for edge: " + edge);
				}
			}
		}
		
		if(count > 0) { 
			env().log().printError("found " + count + " edges with empty active region");
		}
	}
	
	public void print(TextOutput printer) {
		printer.writer().println(toString());
		for ( Vertex vertex : vertices ) {
			printer.writer().println(vertex.toString());
		}
	}
	
	public void printDot2(TextOutput printer){
		printer.writer().println("digraph arg {");
		printer.writer().println("\tnode [fontsize=8, shape=circle, color=\"#" + env().stringProperty("CoalescenceColor") + "\", fontname=Sans, height=0.3]");
		printer.writer().println("\tedge [fontsize=8, color=\"#" + env().stringProperty("EdgeColor") + "\", fontname=Sans, arrowsize=\"0.5\" arrowtail=dot]");

		// Define MRCA as source
		printer.writer().println("\t{ node [peripheries=3, color=\"#" + env().stringProperty("MRCAColor") + "\"]; rank=source; " + gmrca().getId() + " [label=\"MRCA\"];}");
		
		// Define mutation nodes
		printer.writer().println("\t{ node [color=\"#" + env().stringProperty("MutationColor") + "\", shape=circle, peripheries=2];");
		
		for(Edge edge : edges){
			if(edge.hasMutations()){
				int mlb = (int)Math.round(Math.sqrt(edge.mutations.size() / Math.log10(snpDomain.cardinality())));
				int lb = 0;
				StringBuilder msb = new StringBuilder();
				for(Mutation mutation : edge){
					if(lb==0) msb.append("\\n");
					else msb.append(", ");
					msb.append(mutation.position());
					if(++lb >= mlb) lb=0;
				}
				msb.delete(0, 2);
				printer.writer().println("\t\tm" + edge.source().getId() + "_" + edge.target().getId() + " [label=\"" +  msb.toString() + "\"];");
			}
		}
		
		printer.writer().println("\t}");

		// Recombinations nodes
		printer.writer().println("\t{ node [color=\"#" + env().stringProperty("RecombinationColor") + "\", shape=circle, peripheries=2];");
		for(Vertex vertex : vertices){
			if(!vertex.isLeaf() && vertex.inDegree() > 1)
				printer.writer().println("\t\t" + vertex.getId() + " [label=\"\"];");
		}
		printer.writer().println("\t}");

		// Recombinations leafs
		printer.writer().println("\t{ node [color=\"#" + env().stringProperty("RecombinationColor") + "\", shape=circle, peripheries=2]; rank=sink; ");
		for(Vertex vertex : vertices){
			if(vertex.isLeaf() && vertex.inDegree() > 1)
				printer.writer().println("\t\t" + vertex.getId() + ";");
		}
		printer.writer().println("\t}");

		// Other leafs
		printer.writer().println("\t{ node [peripheries=2]; rank=sink;");
		for(Vertex vertex : vertices){
			if(vertex.isLeaf() && vertex.inDegree() < 2)
				printer.writer().println("\t\t" + vertex.getId() + ";");
		}
		printer.writer().println("\t}");

		// Coalescence
		printer.writer().println("\t{");
		for(Vertex vertex : vertices){
			if(!vertex.isLeaf() && !vertex.isGmrca() && vertex.inDegree() < 2)
				printer.writer().println("\t\t" + vertex.getId() + " [label=\"\"];");
		}
		printer.writer().println("\t}");
		
		for(Edge edge : edges){
			if(edge.hasMutations()){
				printer.writer().println("\t\t" + edge.source().getId() + "->m" + edge.source().getId() + "_" + edge.target().getId() + "[label=\"" + edge.activeRegion().toCompactString() + "\"]" + ";");
				printer.writer().println("\t\t" + "m" + edge.source().getId() + "_" + edge.target().getId() + "->" + edge.target().getId() + "[label=\"" + edge.activeRegion().toCompactString() + "\"]" + ";");				
			} else {
				printer.writer().println("\t\t" + edge.source().getId() + "->" + edge.target().getId() + "[label=\"" + edge.activeRegion().toCompactString() + "\"]" + ";");				
			}
		}
		printer.writer().println("}");
	}
	
	public void printDot1(TextOutput printer){
		printer.writer().println("digraph arg {");
//		printer.writer().println("\tranksep=0.5;");
//		printer.writer().println("\tratio=1;");
//		printer.writer().println("\tremincross=true;");
		printer.writer().println("\tnode [fontsize=8, shape=circle, color=\"#" + env().stringProperty("CoalescenceColor") + "\", fontname=Sans, height=0.3]");
		printer.writer().println("\tedge [fontsize=8, color=\"#" + env().stringProperty("EdgeColor") + "\", fontname=Sans, arrowsize=\"0.5\" arrowtail=dot]");

		// Define MRCA as source
		printer.writer().println("\t{ node [peripheries=3, color=\"#" + env().stringProperty("MRCAColor") + "\"]; rank=source; " + gmrca().getId() + " [label=\"MRCA\"];}");
		
		// Define mutation nodes
		printer.writer().println("\t{ node [color=\"#" + env().stringProperty("MutationColor") + "\", shape=circle, peripheries=2];");
		
		for(Edge edge : edges){
			if(edge.hasMutations()){
				int mlb = (int)Math.round(Math.sqrt(edge.mutations.size() / Math.log10(snpDomain.cardinality())));
				int lb = 0;
				StringBuilder msb = new StringBuilder();
				for(Mutation mutation : edge){
					if(lb==0) msb.append("\\n");
					else msb.append(", ");
					msb.append(mutation.position());
					if(++lb >= mlb) lb=0;
				}
				msb.delete(0, 2);
				printer.writer().println("\t\tm" + edge.source().getId() + "_" + edge.target().getId() + " [label=\"" +  msb.toString() + "\"];");
			}
		}
		
		printer.writer().println("\t}");

		// Recombinations nodes
		printer.writer().println("\t{ node [color=\"#" + env().stringProperty("RecombinationColor") + "\", shape=Mrecord, height=0.4];");
		for(Vertex vertex : vertices){
			if(vertex.inDegree() > 1){
				StringBuilder sb = new StringBuilder();
				for(Edge e : vertex.edges){
					if(vertex.isTarget(e)){
						sb.append("|<");
						sb.append(e.activeRegion().toCompactString());
						sb.append(">");
						sb.append(e.activeRegion().toCompactString());
					}
				}
				sb.delete(0, 1);
				printer.writer().println("\t\t" + vertex.getId() + " [label=\"" + sb.toString() + "\"];");
			}
		}
		printer.writer().println("\t}");

		// leafs
		printer.writer().println("\t{ node [peripheries=2]; rank=sink;");
		for(Vertex vertex : vertices){
			if(vertex.isLeaf()){
				if(vertex.inDegree() > 1){
					printer.writer().println("\t\tL" + vertex.getId() + "[label=\"" + vertex.getId() + "\"];");
				} else {
					printer.writer().println("\t\t" + vertex.getId() + "[label=\"" + vertex.getId() + "\"];");
				}
			}
		}
		printer.writer().println("\t}");

		// leafs
		for(Vertex vertex : vertices){
			if(vertex.isLeaf() && vertex.inDegree() > 1){
				printer.writer().println("\t\t" + vertex.getId() + "->L" + vertex.getId() + ";");
			}
		}

		// Coalescence
		printer.writer().println("\t{");
		for(Vertex vertex : vertices){
			if(!vertex.isGmrca() && vertex.inDegree() < 2 && !vertex.isLeaf())
				printer.writer().println("\t\t" + vertex.getId() + " [label=\"\"];");
		}
		printer.writer().println("\t}");
		
		for(Edge edge : edges){
			if(edge.target().inDegree() > 1){
				if(edge.hasMutations()){
					printer.writer().println("\t\t" + edge.source().getId() + "->m" + edge.source().getId() + "_" + edge.target().getId() + ";");
					printer.writer().println("\t\t" + "m" + edge.source().getId() + "_" + edge.target().getId() + "->\"" + edge.target().getId() + "\":\""+ edge.activeRegion().toCompactString() + "\";");				
				} else {
					printer.writer().println("\t\t" + edge.source().getId() + "->\"" + edge.target().getId() + "\":\""+ edge.activeRegion().toCompactString() + "\";");				
				}
				
			} else {
				if(edge.hasMutations()){
					printer.writer().println("\t\t" + edge.source().getId() + "->m" + edge.source().getId() + "_" + edge.target().getId() + ";");
					printer.writer().println("\t\t" + "m" + edge.source().getId() + "_" + edge.target().getId() + "->" + edge.target().getId() + ";");				
				} else {
					printer.writer().println("\t\t" + edge.source().getId() + "->" + edge.target().getId() + ";");				
				}
			}
		}
		printer.writer().println("}");
	}
}
