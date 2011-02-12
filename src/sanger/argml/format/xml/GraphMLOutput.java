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

package sanger.argml.format.xml;

import java.io.FileNotFoundException;
import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.argml.graph.model.Edge;
import sanger.argml.graph.model.Genealogy;
import sanger.argml.graph.model.Mutation;
import sanger.argml.graph.model.Vertex;
import sanger.argml.io.XmlOutput;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetException;

public class GraphMLOutput extends Environmental{
	private XmlOutput out;
	private enum HaplotypeVerbosity {ROOT, LEAF, ALL, NONE}
	
	private boolean nodeActivity;
	private HaplotypeVerbosity haplotypeVerbosity;
	
	private static final String GRAPHML_ARGML_EDGE_EXTENSION_KEY = "d0";
	private static final String GRAPHML_ARGML_NODE_EXTENSION_KEY = "d1";
	
	private static final String ARGML_HAPLOTYPE_QNAME = "haplotype";
	private static final String XML_ID_QNAME = "id";
	private static final String LINE_SEPARATOR = System.getProperty("line.separator");
	
	private static final String ARGML_NS="http://sanger.ac.uk/xsd/argml.xsd";
	private static final String ARGML_NS_PREFIX = "s";
	private static final String ARGML_TOPOLOGY_QNAME = "topology";
	private static final String ARGML_ACTIVITY_QNAME = "activity";
	private static final String ARGML_MAP_QNAME = "map";
	private static final String ARGNL_RIGHT_QNAME = "right";
	private static final String ARGML_LEFT_QNAME = "left";
	private static final String ARGML_MUTATION_QNAME = "mutation";
	private static final String ARGML_MARKER_QNAME = "marker";

	private static final String GRAPHML_NS = "http://graphml.graphdrawing.org/xmlns";
	private static final String GRAPHML_GRAPHML_QNAME = "graphml";
	private static final String GRAPHML_GRAPH_QNAME = "graph";
	private static final String GRAPHML_NODE_QNAME = "node";
	private static final String GRAPHML_EDGE_QNAME = "edge";
	private static final String GRAPHML_DATA_QNAME = "data";
	private static final String GRAPHML_PARSE_OUTDEGREE_QNAME = "parse.outdegree";
	private static final String GRAPHML_PARSE_INDEGREE_QNAME = "parse.indegree";
	private static final String GRAPHML_KEY_QNAME = "key";
	private static final String GRAPHML_TARGET_QNAME = "target";
	private static final String GRAPHML_SOURCE_QNAME = "source";

	public GraphMLOutput(Environment env, XmlOutput out){
		super(env);
		this.out =  out;
		this.nodeActivity = env.flag("NodeActivity");
		this.haplotypeVerbosity = Enum.valueOf(HaplotypeVerbosity.class, env.stringProperty("HaplotypeVerbosity").toUpperCase());
	}

	public void writeStart() throws XMLStreamException, FileNotFoundException{
		out.writer().writeStartDocument();
		out.writer().writeCharacters(LINE_SEPARATOR);
		out.writer().setDefaultNamespace(GRAPHML_NS);
		out.writer().writeStartElement(GRAPHML_NS, GRAPHML_GRAPHML_QNAME);
		out.writer().writeDefaultNamespace(GRAPHML_NS);
		out.writer().writeNamespace(ARGML_NS_PREFIX,ARGML_NS);
	}

	public void writeEnd()throws XMLStreamException, IOException{
		out.writer().writeEndElement();
		out.writer().writeEndDocument();
	}
	
	public void write(Genealogy genealogy) throws XMLStreamException, NaturalSetException{
		out.writer().writeStartElement(GRAPHML_NS, GRAPHML_GRAPH_QNAME);
			out.writer().writeAttribute("edgedefault", "directed");
			out.writer().writeAttribute("parse.nodes", Integer.toString(genealogy.vertices().size()));
			out.writer().writeAttribute("parse.edges", Integer.toString(genealogy.edges().size()));
			out.writer().writeCharacters(LINE_SEPARATOR);
			
			for ( Vertex vertex : genealogy.vertices()) {
				writeNode(vertex);
			}

			for ( Edge e : genealogy.edges()) {
				writeEdge(e);
			}
		out.writer().writeEndElement();
	}
	
	private void writeNode(Vertex v) throws XMLStreamException, NaturalSetException{
		out.writer().writeStartElement(GRAPHML_NS, GRAPHML_NODE_QNAME);
		out.writer().writeAttribute(XML_ID_QNAME, 'n' + Integer.toString(v.getId()));
		out.writer().writeAttribute(GRAPHML_PARSE_INDEGREE_QNAME, Integer.toString(v.inDegree()));
		out.writer().writeAttribute(GRAPHML_PARSE_OUTDEGREE_QNAME, Integer.toString(v.outDegree()));

			if(extendGraphMLNode(v)){
				out.writer().writeStartElement(GRAPHML_NS, GRAPHML_DATA_QNAME);
				out.writer().writeAttribute(GRAPHML_KEY_QNAME, GRAPHML_ARGML_NODE_EXTENSION_KEY);
					if(nodeActivity) { writeInterval(v.outActiveRegion()); }
					if(includeHaplotype(v)) { writeHaplotype(v); }
				out.writer().writeEndElement();
			}
			
		out.writer().writeEndElement();
		out.writer().writeCharacters(LINE_SEPARATOR);
	}

	private void writeHaplotype(Vertex v) throws XMLStreamException, NaturalSetException{	
		out.writer().writeStartElement(ARGML_NS, ARGML_HAPLOTYPE_QNAME);
			out.writer().writeCharacters(v.haplotype().toBinaryString());
		out.writer().writeEndElement();
	}

	private void writeEdge(Edge e) throws XMLStreamException{
		out.writer().writeStartElement(GRAPHML_NS, GRAPHML_EDGE_QNAME);
			out.writer().writeAttribute(GRAPHML_SOURCE_QNAME, 'n' + Integer.toString(e.source().getId()));
			out.writer().writeAttribute(GRAPHML_TARGET_QNAME, 'n' + Integer.toString(e.target().getId()));			
			out.writer().writeStartElement(GRAPHML_NS, GRAPHML_DATA_QNAME);
				out.writer().writeAttribute(GRAPHML_KEY_QNAME, GRAPHML_ARGML_EDGE_EXTENSION_KEY);
				writeInterval(e.activeRegion());
				writeMutations(e);
			out.writer().writeEndElement();		
		out.writer().writeEndElement();
		out.writer().writeCharacters(LINE_SEPARATOR);
	}

	private void writeInterval(NaturalSet i) throws XMLStreamException{		
		out.writer().writeStartElement(ARGML_NS, ARGML_ACTIVITY_QNAME);
			if(i.isEmpty() || !i.isContinuous()){ out.writer().writeAttribute(ARGML_TOPOLOGY_QNAME, i.topology().toString()); }
			if(!i.isEmpty()){
				out.writer().writeAttribute(ARGML_LEFT_QNAME, Integer.toString(i.min()));
				out.writer().writeAttribute(ARGNL_RIGHT_QNAME, Integer.toString(i.max()));
				if(!i.isContinuous()){
					out.writer().writeStartElement(ARGML_NS, ARGML_MAP_QNAME);
						out.writer().writeCharacters(i.toBinaryString());
					out.writer().writeEndElement();
				}
			}
		out.writer().writeEndElement();		
	}
	
	private void writeMutations(Edge e) throws XMLStreamException{
		for(Mutation m : e){
			out.writer().writeStartElement(ARGML_NS, ARGML_MUTATION_QNAME);
			out.writer().writeAttribute(ARGML_MARKER_QNAME, String.valueOf(m.position()));
			out.writer().writeEndElement();
		}
	}

	
	public boolean extendGraphMLNode(final Vertex vertex) {
		return 
			nodeActivity || 
			(vertex.isLeaf() && haplotypeVerbosity == HaplotypeVerbosity.LEAF) || 
			(vertex.isGmrca() && haplotypeVerbosity == HaplotypeVerbosity.ROOT) ||
			includeHaplotype(vertex);
	}
	
	public boolean includeHaplotype(final Vertex vertex) {
		return 
			(vertex.isLeaf() && haplotypeVerbosity == HaplotypeVerbosity.LEAF) || 
			(vertex.isGmrca() && haplotypeVerbosity == HaplotypeVerbosity.ROOT) ||
			haplotypeVerbosity == HaplotypeVerbosity.ALL;			
	}	
}
