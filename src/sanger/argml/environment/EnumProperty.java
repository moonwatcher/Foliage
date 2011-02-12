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

public class EnumProperty extends StringProperty implements Iterable<String> {
	private LinkedList<String> enumeration;
	private static String enumerateExp(LinkedList<String> enumeration){
		StringBuilder sb = new StringBuilder();
		for(String e : enumeration){
			sb.append('|');
			sb.append(e);
		}
		if(sb.length() > 0) sb.delete(0, 1);
		return sb.toString();
	}
	
	public EnumProperty(String name, String identifier, String value, LinkedList<String> enumeration, String help, int order) {
		super(name, identifier, value, help, order, enumerateExp(enumeration));
		this.enumeration = enumeration;
	}
	
	public Iterator<String> iterator(){
		return enumeration.iterator();
	}
	
	public String help(int pad){
		StringBuilder sb = new StringBuilder();
		sb.append(super.help(pad));
		sb.append("{ ");
		for(String literal : this){
			sb.append(literal);
			sb.append(", ");
		}
		sb.delete(sb.length() - 2, sb.length());
		sb.append(" }");
		return sb.toString();
	}
}
