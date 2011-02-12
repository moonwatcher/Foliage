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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;


public class XmlOutput extends Environmental{
	protected static final XMLOutputFactory factory = XMLOutputFactory.newInstance();
	public static final String LINE_SEPARATOR = System.getProperty("line.separator");
	
	protected File file;
	protected XMLStreamWriter writer = null;
	protected OutputStream stream = null;
	
	public XmlOutput(Environment env, OutputStream out) throws XMLStreamException{
		super(env);
		writer = factory.createXMLStreamWriter(System.out, "utf-8");
	}
	
	public XmlOutput(Environment env, String pathname) throws XMLStreamException, FileNotFoundException{
		super(env);
		this.file = new File(env.outbase(), pathname);
		if(file.getParentFile().exists() ? true : file.getParentFile().mkdirs()) {
			stream = new FileOutputStream(file);
			writer = factory.createXMLStreamWriter(stream, "utf-8");			
		}					
	}
	
	public void close() throws XMLStreamException, IOException {
		writer.flush();
		writer.close();
		if(file!=null) stream.close();
	}

	public XMLStreamWriter writer() {
		return writer;
	}

	public String toString(){
		return file!=null?file.toString():"stdout"; 
	}
	
}
