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

package sanger.argml.environment;

import java.util.Iterator;
import java.util.LinkedList;

public class Instruction implements Comparable<Instruction>, Iterable<Dependency>{
	protected String name;
	protected String help;
	protected String input;
	protected String output;
	protected int order;
	protected LinkedList<Dependency> dependencies = new LinkedList<Dependency>();
	
	public Instruction(String name, int order){
		this.name = name;
		this.order = order;
		this.help = "";
	}
	
	public void addDependency(Dependency dependency){
		dependencies.add(dependency);
	}

	public int compareTo(Instruction other){
		return this.order - other.order;
	}	
	
	public String toString(){
		StringBuilder display = new StringBuilder();
		display.append(name);
		for(Dependency d : this){
			if(d.getProperty().hasValue()){
				display.append(' ');
				display.append(d.getProperty().symbol);
				display.append(' ');
				display.append(d.getProperty().stringValue());
			}
		}
		return display.toString();		
	}
		
	public String help(){
		int pad = 0;
		for(Dependency d : this){ pad = Math.max(pad, d.getProperty().symbol.length()); }
		pad = Math.max(pad, name.length());
		StringBuilder display = new StringBuilder();
		display.append(' ');
		display.append(pad(name, pad + 2));
		display.append(help);
		display.append("\n");
		if(!dependencies.isEmpty()) {
			for(Dependency d : this){
				display.append("\n ");
				display.append(d.isOptional() ? "- " : "+ ");
				display.append(d.getProperty().help(pad));
			}
			display.append("\n");
		}

		return display.toString();	
	}

	public String description(int pad){
		StringBuilder display = new StringBuilder();
		display.append(' ');
		display.append(pad(name, pad));
		display.append(help);
		return display.toString();	
	}

	public Iterator<Dependency> iterator(){
		return dependencies.iterator();
	}
	
	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the order
	 */
	public int getOrder() {
		return order;
	}

	private String pad(String s, int pad){
		char [] buf = new char[pad + 3 - s.length()];
		for (int i=0; i<buf.length ; i++) buf[i] = ' ';
		return s + String.valueOf(buf);
	}
}
