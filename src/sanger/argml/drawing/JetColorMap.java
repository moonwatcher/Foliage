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

package sanger.argml.drawing;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import sanger.argml.environment.Environment;


public class JetColorMap {
	private ArrayList<Color> scale;
	private int size;
	private static final Color offTheChartColor = new Color(0x00FFFFFF, true);
	private static final String offTheChartString = Integer.toHexString(0x00FFFFFF);
	
	public JetColorMap(Environment env){
		InputStream stream = null;
		BufferedReader reader = null;

		scale = new ArrayList<Color>();
		try {
			stream = this.getClass().getResourceAsStream("jcols");
			reader = new BufferedReader(new InputStreamReader(stream));
			
			String line = reader.readLine();
			while(line != null){
				scale.add(Color.decode(line));
				line = reader.readLine();
			}
			scale.trimToSize();
			size = scale.size() - 1; 
		
		} catch (IOException e) {
			env.log().printError(e);
			
		} finally{
			try {
				if(reader != null) reader.close();
				if(stream != null) stream.close();
				
			} catch (IOException e) {
				env.log().printError(e);
			}
		}
	}
	
	public String mapToHexColor(double x){
		String value = null;
		if(x <= 1.0){
			value = Integer.toHexString( scale.get((int)Math.round(x * size)).getRGB() & 0x00ffffff );
		} else {
			value = offTheChartString;
		}
		return value;
	}
	
	public Color mapToColor(double x){
		Color value = null;
		if(x<=1.0){
			value = scale.get((int)Math.round(x * size));
		} else {
			value = offTheChartColor;
		}
		return value;
	}

	
	public ArrayList<Color> scale() {
		return scale;
	}

	public int size() {
		return scale.size();
	}
 
}
