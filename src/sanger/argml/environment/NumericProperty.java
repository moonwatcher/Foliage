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

import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class NumericProperty extends AbstractProperty {
	private Double value = null;
	
	public NumericProperty(String name, String symbol, Double value, String help, int order){
		super(name, symbol, help, order);
		pattern = Pattern.compile("[\\s]+(?:" + symbol + ")[\\s]+([0-9\\.]+)[\\s]+");
		this.value = value;
	}
	
	public void match(String command){
		Matcher matcher = pattern.matcher(command); 
		if(matcher.find()){ value = Double.parseDouble(matcher.group(1)); }	
	}

	public String stringValue(){
		return String.valueOf(value);
	}
	
	public Double getValue() {
		return value;
	}
	
	public String toString(){
		return name + " : " + String.valueOf(value);
	}
	
	public boolean hasValue(){
		return (value != null);
	}
}
