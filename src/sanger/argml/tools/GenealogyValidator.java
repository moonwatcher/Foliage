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

package sanger.argml.tools;

import java.io.IOException;

import org.apache.solr.util.OpenBitSet;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.argml.graph.model.Genealogy;
import sanger.argml.graph.model.HaplotypeSet;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetCollection;
import sanger.math.set.NaturalSetException;

public class GenealogyValidator extends Environmental{
	private NaturalSetCollection errorSequences;
	private NaturalSet errorLocations;
	private NaturalSet ancestralSequence;
	private int errorCount;
	
	public GenealogyValidator(Environment env, Genealogy g, HaplotypeSet haplotypes) throws NaturalSetException {
		super(env);
		if(!haplotypes.hasMissing()){
			errorSequences = g.mutations();
			errorSequences.lineByLineXor(haplotypes.haplotypes());
			errorSequences = NaturalSetCollection.transpose(errorSequences);
			OpenBitSet error = new OpenBitSet(haplotypes.snpDomain().closureCardinality());
			OpenBitSet ancestor = new OpenBitSet(haplotypes.snpDomain().closureCardinality());
			
			int numberOfHaplotypes = errorSequences.domain().cardinality();
				
			for(int i = 0; i<=g.snpDomain().last(); i++) {
				NaturalSet line = errorSequences.get(i);
				int cardinality = line.cardinality();
				
				if(cardinality == numberOfHaplotypes){ ancestor.fastSet(i); } 
				else if (cardinality > 0) { error.fastSet(i); }
			}
			this.errorLocations = new NaturalSet(haplotypes.snpDomain(), error);
			this.ancestralSequence = new NaturalSet(haplotypes.snpDomain(), ancestor);
			this.errorCount = errorLocations.cardinality();
		} else {
			throw new NaturalSetException("Validating missing values not supported.");
		}
	}
	
	public boolean hasErrors(){
		return (errorCount > 0);
	}

	public NaturalSet ancestralSequence() {
		return ancestralSequence;
	}

	public NaturalSetCollection errorSequences() {
		return errorSequences;
	}

	public int errorCount() {
		return errorCount;
	}

	public NaturalSet errorVector() {
		return errorLocations;
	}	

	public void print() throws IOException {
		try {
			if(hasErrors()) {
				env().log().printError("errors found on " + errorCount + " columns!");
				env().log().println(errorLocations.toBinaryString());
				NaturalSetCollection.transpose(errorSequences).print(env().log().writer());			
			} else {
				env().log().printInfo("genealogy valid");
			}
		} catch (NaturalSetException e) {
			env().log().printError(e);
		}
	}
}
