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

import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Iterator;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.argml.io.TextInput;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetException;

public class GenealogyFactory extends Environmental{
	protected TextInput input;
	protected NaturalDomain snpDomain = null;
	protected NaturalDomain haplotypeDomain = null;
	protected NaturalSet snpFilter = null;
	protected HashMap<Integer, Vertex> stubs = null;		
	protected Genealogy genealogy = null;
	protected boolean multifurcated;

	public GenealogyFactory(Environment env, TextInput input) throws FileNotFoundException{
		super(env);
		this.input = input;
		multifurcated = env.flag("MultiFurcate");
	}
	
	public void start(){
		genealogy = new Genealogy(input.env(), snpDomain, haplotypeDomain);
		stubs = new HashMap<Integer, Vertex>(haplotypeDomain.cardinality());
		for ( Vertex v : genealogy.vertices() ) {
			stubs.put(v.getId(), v);
		}			
	}
	
	public Genealogy finish(){
		Genealogy g = null;
		try {
			genealogy.sort();
			findGmrca();
			g = Genealogy.clipSnpDomain(genealogy, snpFilter);
			
			if(multifurcated){
				int i = 0;
				for ( Vertex vertex : g.vertices() ) {
					vertex.setId(i);
					i++;
				}
			}
			
		} catch (Exception e) {
			env().log().printError(e);
			
		} finally {
			stubs = null;
			genealogy = null;
		
		} return g;		
	}
		
	public void mutate(int sourceKey, int targetKey, int marker) throws NaturalSetException {
		Vertex target = stubs.get(targetKey);
		stubs.remove(targetKey);

		if(multifurcated){					
			if ( target.outDegree() == 1 ) {
				for ( Edge edge : target ) {
					if ( target.isSource(edge) ) {
						edge.mutate(marker);
						target.setId(sourceKey);
						stubs.put(sourceKey, target);
					}
				}
				
			} else {

				Edge edge = target.createAncestor(sourceKey, target.outActiveRegion());
				edge.mutate(marker);
				stubs.put(sourceKey, edge.source());
			}
			
		} else {
			Edge edge = target.createAncestor(sourceKey);
			edge.mutate(marker);
			stubs.put(sourceKey, edge.source());
		}
	}
	
	public void recombine(int leftKey, int rightKey, int childKey, int position) throws NaturalSetException {
		Vertex left, right;
		Vertex child = stubs.get(childKey);
		stubs.remove(childKey);

		if(multifurcated){
						
			if ( child.outDegree() == 1 ) {
				left = child;
				left.setId(leftKey);
				right = left.split(rightKey, position);
				
			} else {
				NaturalSet activity = child.outActiveRegion();
				left = child.createAncestor(leftKey, activity.copyUpTo(position)).source();
				activity.removeBelow(position + 1);
				right = child.createAncestor(rightKey, activity).source();	
			}
			
		} else {
			NaturalSet activity = child.outActiveRegion();			
			left = child.createAncestor(leftKey, activity.copyUpTo(position)).source();
			activity.removeBelow(position + 1);
			right = child.createAncestor(rightKey, activity).source();
			
		}

		stubs.put(leftKey, left);
		stubs.put(rightKey, right);
	}
	
	public void coalesce(int parentKey, int oneKey, int twoKey) throws NaturalSetException {
		Vertex one = stubs.get(oneKey);
		Vertex two = stubs.get(twoKey);
		stubs.remove(oneKey);
		stubs.remove(twoKey);
		
		if(multifurcated){
			if ( !one.isLeaf() && !two.isLeaf() ) { // niether nodes is a leaf
				Vertex parent = one.merge(parentKey, two);
				stubs.put(parentKey, parent);
				
			} else {
				if ( one.isLeaf() && two.isLeaf() ){ // both nodes are leafs
					Vertex parent = one.createAncestor(parentKey).source();
					two.connectToAncestor(parent);					
					stubs.put(parentKey, parent);
					
				} else { // One of the decendents is a leaf and the other is not
					Vertex leaf = one.isLeaf() ? one : two;
					Vertex internal = one.isLeaf() ? two : one;
					internal.setId(parentKey);
					leaf.connectToAncestor(internal);
					stubs.put(parentKey, internal);
				}
			}
			
		} else {
			Vertex parent = one.createAncestor(parentKey).source();
			two.connectToAncestor(parent);			
			stubs.put(parentKey, parent);	
		}
	}
	
	private void findGmrca() throws Exception{
		if(stubs.size() > 1){ // Need to simulate a root, multi furcating for now.
			int rootId = genealogy.vertices().get(genealogy.vertices().size() - 1).getId() + 1;
			Iterator<Vertex> iterate = stubs.values().iterator();
			genealogy.gmrca = iterate.next().createAncestor(rootId).source();
			while (iterate.hasNext()){
				iterate.next().connectToAncestor(genealogy.gmrca());
			}			
		} else {
			for ( Vertex vertex : stubs.values() ) {
				genealogy.gmrca = vertex;
			}
		}		
	}
	
	public void filterSnp(NaturalSet region) throws NaturalSetException{
		this.snpFilter = snpDomain.project(region);
	}

	
	/**
	 * @return the snpDomain
	 */
	public NaturalDomain getSnpDomain() {
		return snpDomain;
	}

	
	/**
	 * @return the haplotypeDomain
	 */
	public NaturalDomain getHaplotypeDomain() {
		return haplotypeDomain;
	}

			
}