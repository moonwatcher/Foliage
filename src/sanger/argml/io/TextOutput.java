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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;


public class TextOutput extends Environmental{
	protected File file;
	protected FileWriter fileWriter;
	protected PrintWriter printWriter;
	
	protected TextOutput(Environment env){
		super(env);
	}
	
	public TextOutput(Environment env, OutputStream out){
		super(env);
		this.printWriter = new PrintWriter(out);
	}
	
	public TextOutput(Environment env, String pathname) throws IOException{
		super(env);
		file = new File(env.outbase(), pathname);
		if(file.getParentFile().exists() ? true : file.getParentFile().mkdirs()) {
			fileWriter = new FileWriter(file, false);
			printWriter = new PrintWriter(fileWriter, false);				
		}			
	}

	public void close() throws IOException {
		printWriter.flush();
		printWriter.close();
		if(file!=null) fileWriter.close();
	}
	
	public PrintWriter writer(){
		return printWriter;
	}

	public String toString(){
		return file!=null?file.toString():"stdout"; 
	}	
}
