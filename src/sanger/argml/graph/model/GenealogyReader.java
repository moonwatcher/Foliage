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

import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sanger.argml.environment.Environment;
import sanger.argml.format.IllegalFileFormatException;
import sanger.argml.io.TextInput;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetException;

public class GenealogyReader extends GenealogyFactory implements Iterable<Genealogy>{
	public enum MargaritaGenealogyReadState { START, ARGINFERENCE, ARGS, READY, WAITING, BUILDING, EOF; }
	protected MargaritaGenealogyReadState state = MargaritaGenealogyReadState.START;	
	protected static final Pattern coalescencePattern = Pattern.compile("^[0-9]+\\sco\\s([0-9]+)\\s([0-9]+)\\s([0-9]+)$");
	protected static final Pattern mutationPattern = Pattern.compile("^[0-9]+\\s(?:mu|er)\\s([0-9]+)\\s([0-9]+)\\s([0-9]+)$");
	protected static final Pattern recombinationPattern = Pattern.compile("^[0-9]+\\sre\\s([0-9]+)\\s([0-9]+)\\s([0-9]+)\\s([0-9]+)$");
	protected static final Pattern headerPattern = Pattern.compile("^([0-9]+)\\s([0-9]+)\\s(?:(?:[0-9\\.E-]+|NA)\\s){6}(?:[0-9\\.E-]+|NA)$");
	protected static final Pattern nextArgPattern = Pattern.compile("^ARG ([0-9]+)$");
	protected static final Pattern argInferencePattern = Pattern.compile("^%ARGINFERENCE$");
	protected static final Pattern argsPattern = Pattern.compile("^%ARGS");

	private NaturalDomain argDomain = null;
	private NaturalSet argFilter = null;
	private Integer nextIndex = -1;
	private Integer index = null;
	private Integer count = 0;
	
	
	public void start(){
		super.start();
		this.state = MargaritaGenealogyReadState.BUILDING;
	}
	
	public GenealogyReader(Environment env, TextInput input) throws NaturalSetException, IOException, IllegalFileFormatException{
		super(env, input);
		while(state != MargaritaGenealogyReadState.ARGS){
			processInstruction(input.reader().readLine());
		}
	}
	
	private void processInstruction(String instruction) throws NaturalSetException, IOException, IllegalFileFormatException{
		Matcher matcher = null;
		if(instruction != null || instruction == ""){
			if(state == MargaritaGenealogyReadState.START){
				matcher = argInferencePattern.matcher(instruction);
				if(matcher.matches()){
					state = MargaritaGenealogyReadState.ARGINFERENCE;
				}
				
			} else if (state == MargaritaGenealogyReadState.ARGINFERENCE) {
				matcher = headerPattern.matcher(instruction);
				if(matcher.matches()){
					if(count == 0){
						snpDomain = new NaturalDomain(Integer.parseInt(matcher.group(2))); 
						haplotypeDomain = new NaturalDomain(Integer.parseInt(matcher.group(1)));
						snpFilter = snpDomain.createCompleteNaturalSet();
					} count++;
				} else {
					matcher = argsPattern.matcher(instruction);
					if(matcher.matches()){
						argDomain = new NaturalDomain(count);
						argFilter = argDomain.createCompleteNaturalSet();
						state = MargaritaGenealogyReadState.ARGS;
					}
				}
				
			} else {
				matcher = coalescencePattern.matcher(instruction);
				if(!matcher.matches()){
					matcher = mutationPattern.matcher(instruction);
					if(!matcher.matches()){
						matcher = recombinationPattern.matcher(instruction);
						if(!matcher.matches()){
							matcher = nextArgPattern.matcher(instruction);
							if(matcher.matches()){
								index = nextIndex;
								nextIndex = Integer.parseInt(matcher.group(1));
								if(nextIndex > argFilter.max()) {
									state = MargaritaGenealogyReadState.EOF;
									input.close();
									
								} else {
									try {
										while(!argFilter.contains(nextIndex)){
											do {
												matcher = nextArgPattern.matcher(input.reader().readLine());
											} while (!matcher.matches());
											nextIndex = Integer.parseInt(matcher.group(1));
										}
										state = MargaritaGenealogyReadState.WAITING;
										
									} catch (Exception e) {
										throw new IllegalFileFormatException(e);
									}
								}									
							}
							
						} else {
							Integer childKey = Integer.decode(matcher.group(1));
							Integer leftKey = Integer.decode(matcher.group(2));
							Integer rightKey = Integer.decode(matcher.group(3));
							Integer position = Integer.decode(matcher.group(4));
							recombine(leftKey, rightKey, childKey, position);
						}
					} else {
						Integer targetKey = Integer.decode(matcher.group(1));
						Integer sourceKey = Integer.decode(matcher.group(2));
						Integer marker = Integer.decode(matcher.group(3));
						mutate(sourceKey, targetKey, marker);
					}
				} else {
					Integer oneKey = Integer.decode(matcher.group(1));
					Integer twoKey = Integer.decode(matcher.group(2));
					Integer parentKey = Integer.decode(matcher.group(3));
					coalesce(parentKey, oneKey, twoKey);
				}				
			}			
		} else {
			state = MargaritaGenealogyReadState.EOF;
			input.close();
			index++;
		}
	}

	public Iterator<Genealogy> iterator(){
		return new GenealogyIterator(this);
	}
	
	public class GenealogyIterator implements Iterator<Genealogy> {
		GenealogyReader factory;

		public GenealogyIterator(GenealogyReader factory) {
			this.factory = factory;
		}
		
		public boolean hasNext() {
			if(factory.state == MargaritaGenealogyReadState.ARGS){
				while(factory.state == MargaritaGenealogyReadState.ARGS){
					try { processInstruction(input.reader().readLine());
					} catch (Exception e) { input.env().log().printError(e);}					
				}
			}
			return (factory.state == MargaritaGenealogyReadState.WAITING);
		}

		public Genealogy next() {
			Genealogy result = null;
			if(hasNext()){
				result = factory.readNext();
			} else {
				throw new NoSuchElementException("No more ARGs in file.");
			}
			return result;
		}

		public void remove() {}
	}

	public Integer index() {
		return index;
	}

	private Genealogy readNext() {
		Genealogy g = null;
		try {
			start();			
			while(state == MargaritaGenealogyReadState.BUILDING){ processInstruction(input.reader().readLine()); }
			g = finish();

		} catch (Exception e) {
			input.env().log().printError(e);
			
		} return g;
	}	

	public void filterArg(NaturalSet region) throws NaturalSetException{
		this.argFilter = argDomain.project(region);
	}
	
	/**
	 * @return the argDomain
	 */
	public NaturalDomain getArgDomain() {
		return argDomain;
	}

}
