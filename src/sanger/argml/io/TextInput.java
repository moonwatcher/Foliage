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

package sanger.argml.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;


public class TextInput extends Environmental{
	protected File file;
	protected InputStream stream = null;
	protected InputStreamReader streamReader = null;
	protected BufferedReader reader = null;
		
	public TextInput(Environment env, String pathname) throws FileNotFoundException{
		super(env);
		this.file = new File(env.inbase(), pathname);
		stream = new FileInputStream(file);
		streamReader = new InputStreamReader(stream);
		reader = new BufferedReader(streamReader);
	}

	public TextInput(Environment env, InputStream stream){
		super(env);
		this.stream = stream;
		streamReader = new InputStreamReader(stream);
		reader = new BufferedReader(streamReader);
	}
	
	public void close() throws IOException {
		reader.close();
		streamReader.close();
		if(file!=null) stream.close();
	}
	
	public BufferedReader reader(){
		return reader;
	}
	
	public String toString(){
		return file!=null?file.toString():"stdout"; 
	}
}
