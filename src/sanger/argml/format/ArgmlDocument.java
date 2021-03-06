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

package sanger.argml.format;

import javax.xml.stream.XMLStreamException;

import sanger.argml.io.XmlOutput;

public class ArgmlDocument {
	private XmlOutput out;
	
	public ArgmlDocument(XmlOutput out){
		this.out = out;
	}

	public void writeStartDocument() throws XMLStreamException{
		out.writer().writeStartDocument();
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
		out.writer().writeStartElement("argml");
		out.writer().writeCharacters(XmlOutput.LINE_SEPARATOR);
	}
	
	public void writeEndDocument() throws XMLStreamException{
		out.writer().writeEndElement();
		out.writer().writeEndDocument();
	}	
}
