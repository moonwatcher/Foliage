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

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;

import sanger.argml.io.TextOutput;
import sanger.argml.io.XmlInput;


public class Environment implements Iterable<AbstractProperty> {
	public static String banner = "Foliage v0.1PR ARG manipulation library.\nLior Galanti lg8@sanger.ac.uk 2007.\nWellcome Trust Sanger Institute.\n\n";
	public static String usage = " Usage:\t\tfoliage COMMAND [OPTIONS]...\n Help:\t\tfoliage help --topic COMMAND for command syntax.\n Example:\tfoliage help --topic dot\n\n";
	public static String legend = "\n + : mandatory\n - : optional\n";
	public static String commandTitle = " Available commands:\n\n";

	private int instructionLength = 0;
	
	private Pattern commandExp = Pattern.compile("^[\\s]+([^\\s-][^\\s]*)[\\s]+");
	private HashMap<String, AbstractProperty> properties = new HashMap<String, AbstractProperty>();
	private HashMap<String, Instruction> commands = new HashMap<String, Instruction>();
	private Instruction instruction;
	private File outbase;
	private File inbase;
	private LogOutput log;
	
	public Environment(String command) throws XMLStreamException, IOException, UnknownCommandException{
		readConfig();
		Matcher m = commandExp.matcher(command);
		if(m.find()) {
			instruction = commands.get(m.group(1));
			if(instruction == null) instruction = commands.get("help");
		} else instruction = commands.get("help");
		
		for(Dependency d : instruction){
			d.getProperty().match(command);
			if(!d.isOptional() && !d.getProperty().hasValue()){
				throw new UnknownCommandException("error: missing mandatory parameters!", instruction);
			}
		}
		
		outbase = new File(stringProperty("OutputBase"));
		inbase = new File(stringProperty("InputBase"));
		log = stringPropertyExist("Log") ? new LogOutput(this, stringProperty("Log")) : new LogOutput(this, System.err);
		
		if(!instruction.getName().equals("help")){
			log.printInfo("started on " + (new Date(System.currentTimeMillis())).toString());
			log.printInfo("Instruction: " + instruction.name);		
			log.printInfo(this.toString());
			log.printInfo("Command: " + instruction);
		}
		
	}
	
	private void readConfig() throws XMLStreamException{
		XmlInput i = new XmlInput(this, this.getClass().getResourceAsStream("/sanger/argml/cli/config.xml"));
		int order = 0;

		String state = "";
		String name = null;
		String symbol = null;
		String value = null;
		String help = null;
		LinkedList<String> enumeration = null;
		Instruction c = null;
		
		for(int event=i.parser().next(); event!=XMLStreamConstants.END_DOCUMENT; event=i.parser().next()) {
			switch (event) {
				case XMLStreamConstants.START_ELEMENT:
					if(i.parser().getLocalName().equals("properties")) {
						state = "properties";
						order = 0;
						
					} else if (i.parser().getLocalName().equals("command")) {
						state = "command";
						order = 0;
						
					} else if (state.equals("command")) {
						if (i.parser().getLocalName().equals("instruction")){
							c = new Instruction(i.parser().getAttributeValue("", "name"), order++);

						} else if (i.parser().getLocalName().equals("depend")) {
							c.addDependency(
								new Dependency(
									properties.get(i.parser().getAttributeValue("", "name")),
									Boolean.parseBoolean(i.parser().getAttributeValue("", "optional"))
								)
							);

						} else if (i.parser().getLocalName().equals("help")) {
							c.help = formatHelp(i.parser().getElementText());

						} else if (i.parser().getLocalName().equals("input")) {
							c.input = i.parser().getElementText();
						
						} else if (i.parser().getLocalName().equals("output")) {
							c.output = i.parser().getElementText();
						}
						
						
					} else if (state.equals("properties")) {
						if (	i.parser().getLocalName().equals("string")	|| 
								i.parser().getLocalName().equals("integer")	|| 
								i.parser().getLocalName().equals("boolean")	||
								i.parser().getLocalName().equals("double")	||
								i.parser().getLocalName().equals("enum")	) {

							name = i.parser().getAttributeValue("", "name");
							symbol = i.parser().getAttributeValue("", "symbol");
							value = i.parser().getAttributeValue("", "default");
							
						} else if (i.parser().getLocalName().equals("restriction")) {
							enumeration = new LinkedList<String>();
							
						} else if (i.parser().getLocalName().equals("help")) {
							help = formatHelp(i.parser().getElementText());
							
						} else if (i.parser().getLocalName().equals("enumeration")) {
							enumeration.add(i.parser().getAttributeValue("", "value"));
						}
					}
				break;
				
				case XMLStreamConstants.END_ELEMENT:
					if (state.equals("properties")) {
						if (i.parser().getLocalName().equals("string")) {
							properties.put(name, new StringProperty(name, symbol, value, help, order++));
							name = null; symbol = null; value = null; help = null;
							
						} else if (i.parser().getLocalName().equals("integer")) {
							properties.put(name, new NumericProperty(name, symbol, value == null ? null : Math.floor(Double.parseDouble(value)), help, order++));
							name = null; symbol = null; value = null; help = null;
							
						} else if (i.parser().getLocalName().equals("double")) {
							properties.put(name, new NumericProperty(name, symbol, value == null ? null : Double.parseDouble(value), help, order++));
							name = null; symbol = null; value = null; help = null;
							
						} else if (i.parser().getLocalName().equals("boolean")) {
							properties.put(name, new Flag(name, symbol, help, order++));
							name = null; symbol = null; value = null; help = null;
							
						} else if (i.parser().getLocalName().equals("enum")) {
							properties.put(name, new EnumProperty(name, symbol, value, enumeration, help, order++));		
							name = null; symbol = null; value = null; help = null; enumeration = null;
							
						}
						
					} else if(state.equals("command")){
						if (i.parser().getLocalName().equals("instruction")){
							commands.put(c.getName(), c);
							instructionLength = Math.max(instructionLength, c.getName().length());
							c = null;
						}
					}

				break;	
			}
		}
		try { i.close(); } catch (Exception e) { e.printStackTrace(); }
	}

	
	public boolean stringPropertyExist(String name){
		return (stringProperty(name) != null);
	}
	
	public String stringProperty(String name){
		String result = null;
		AbstractProperty property = properties.get(name);
		if(property != null && property instanceof StringProperty) result = ((StringProperty)property).getValue();
		return result;
	}
	
	public boolean numericPropertyExist(String name){
		return (numericProperty(name) != null);
	}
	
	public Double numericProperty(String name){
		Double result = null;
		AbstractProperty property = properties.get(name);		
		if(property != null && property instanceof NumericProperty) result = ((NumericProperty)property).getValue();
		return result;
	}	
	
	public Integer integerProperty(String name){
		Double dresult = null;
		Integer result = null;
		AbstractProperty property = properties.get(name);		
		if(property != null && property instanceof NumericProperty) {
			dresult = ((NumericProperty)property).getValue();
			if(dresult!=null) result = (int)Math.floor(dresult);
		}
		return result;
	}	
		
	public Boolean flag(String name){
		Boolean result = null;
		AbstractProperty property = properties.get(name);		
		if(property != null && property instanceof Flag) result = ((Flag)property).getValue();
		return result;
	}
		
	
	public Iterator<AbstractProperty> iterator(){
		LinkedList<AbstractProperty> order = new LinkedList<AbstractProperty>();
		order.addAll(properties.values());
		Collections.sort(order);
		return order.iterator();
	}

	public LogOutput log() {
		return log;
	}
	
	public File outbase() {
		return outbase;
	}
	
	public File inbase() {
		return inbase;
	}
		
	public String toString(){
		StringBuilder display = new StringBuilder();
		display.append("Properties {\n");
		for(Dependency d : instruction){
			if(d.getProperty().hasValue()){
				display.append("\t");
				display.append(d.getProperty());
				display.append("\n");
			}
		}
		display.append("}");
		return display.toString();	
	}

	public void printCommands(TextOutput out){
		StringBuilder display = new StringBuilder();
		display.append(banner);
		display.append(usage);
		display.append(commandTitle);
		LinkedList<Instruction> sc = new LinkedList<Instruction>(commands.values());
		Collections.sort(sc);
		for(Instruction i : sc){
			display.append(i.description(instructionLength));
			display.append("\n");			
		}
		display.append('\n');
		
		out.writer().print(display.toString());
	}

	public void printTopic(TextOutput out, String topic){
		StringBuilder display = new StringBuilder();
		display.append(banner);
		if(commands.containsKey(topic)) {
			display.append(commands.get(topic).help());
			display.append(legend);
			
		} else {
			display.append("No such topic.\n");
		}
		
		out.writer().print(display.toString());
	}

	public void close() throws IOException{
		if(log != null){
			if(!instruction.getName().equals("help")){
				log.printTotalBenchmark();
				log.printInfo("closing log\n");
			}
			log.close();
		}
	}

	
	/**
	 * @return the instruction
	 */
	public Instruction instruction() {
		return instruction;
	}

	private String formatHelp(String s){
		s = s.replaceAll("[\\s]+", " ");
		s = s.replace(". ", ".\n");
		s = s.trim();
		return s;
	}
	
}
