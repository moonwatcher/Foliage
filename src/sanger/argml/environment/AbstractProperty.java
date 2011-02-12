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

import java.util.regex.Pattern;

public abstract class AbstractProperty implements Comparable<AbstractProperty>{
	protected Pattern pattern;
	protected String name;
	protected String symbol;
	protected String help;
	protected int order;
	
	public AbstractProperty( String name, String symbol, String help, int order){
		this.name = name;
		this.symbol = symbol;
		this.help = help;
		this.order = order;
	}
	
	public int compareTo(AbstractProperty other){
		return this.order - other.order;
	}
	
	public String getName() {
		return name;
	}
	
	public abstract String stringValue();
	
	public abstract void match(String command);
	
	public abstract String toString();
	
	public abstract boolean hasValue();

	public String help(int pad) {
		StringBuilder sb = new StringBuilder();
		sb.append('[');
		sb.append(symbol);
		sb.append("]");
		sb.append(padding(sb.toString(), pad));
		if(help != null)sb.append(help);
		return sb.toString();
	}

	private String padding(String s, int pad){
		char [] buf = new char[pad + 3 - s.length()];
		for (int i=0; i<buf.length ; i++) buf[i] = ' ';
		return String.valueOf(buf);
	}	
}
