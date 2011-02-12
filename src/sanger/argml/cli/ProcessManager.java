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

package sanger.argml.cli;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.imageio.stream.ImageOutputStream;
import javax.xml.stream.XMLStreamException;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.argml.environment.UnknownCommandException;
import sanger.argml.format.IllegalFileFormatException;
import sanger.argml.graph.model.GenealogyReader;
import sanger.argml.graph.model.HaplotypeReader;
import sanger.argml.graph.model.HaplotypeSet;
import sanger.argml.graph.model.Statistics;
import sanger.argml.graph.model.StatisticsFactory;
import sanger.argml.io.TextInput;
import sanger.argml.io.TextOutput;
import sanger.argml.io.XmlInput;
import sanger.argml.io.XmlOutput;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import sanger.math.set.NaturalSetException;

public class ProcessManager extends Environmental{
		
	public ProcessManager(String cmd) 
	throws UnknownCommandException, IOException, IllegalFileFormatException, NaturalSetException, XMLStreamException{
		super(new Environment(cmd));
	}
	
	public TextOutput textOutput() throws IOException{
		return env().stringPropertyExist("Output") ? new TextOutput(env(), env().stringProperty("Output")) : new TextOutput(env(), System.out);
	}
	
	public TextInput textInput() throws FileNotFoundException{
		return env().stringPropertyExist("Input") ? new TextInput(env(), env().stringProperty("Input")) : new TextInput(env(), System.in);
	}
	
	
	public XmlInput xmlInput() throws FileNotFoundException, XMLStreamException{
		return env().stringPropertyExist("Input") ? new XmlInput(env(), env().stringProperty("Input")) : new XmlInput(env(), System.in);
		
	}
	
	public XmlOutput xmlOutput() throws FileNotFoundException, XMLStreamException{
		return env().stringPropertyExist("Output") ? new XmlOutput(env(), env().stringProperty("Output")) : new XmlOutput(env(), System.out);
	}
	
	
	public ImageOutputStream imageOutputStream() throws IOException{
		ImageOutputStream out = null;
		if(env().stringPropertyExist("Output")){
			File file = new File(env().outbase(), env().stringProperty("Output"));
			if(!file.getParentFile().exists()) file.getParentFile().mkdirs();			
			out = ImageIO.createImageOutputStream(file);
		} else {
			out = ImageIO.createImageOutputStream(System.out);
		}

		return out;
	}	
	
	public HaplotypeSet createHaplotypeSet(TextInput input) 
		throws FileNotFoundException, NaturalSetException, IOException {
	
		HaplotypeReader i = new HaplotypeReader(env(), input);
		HaplotypeSet result = i.readHaplotypeSet();
		i.close();
		result = HaplotypeSet.clipSnpDomain(result, filterSNP(result.snpDomain()));
		return result;
	}

	public GenealogyReader createGenealogyFactory(TextInput input) 
		throws IOException, NaturalSetException, IllegalFileFormatException {
		
		GenealogyReader f = new GenealogyReader(env(), input);		
			
		f.filterArg(filterARG(f.getArgDomain()));
		f.filterSnp(filterSNP(f.getSnpDomain()));
		return f;
	}	

	public StatisticsFactory createStatisticsFactory(XmlInput input)
		throws NaturalSetException, XMLStreamException, FileNotFoundException {	
		return createStatisticsFactory(input, null);
	}
	
	public StatisticsFactory createStatisticsFactory(XmlInput input, String filter)
		throws NaturalSetException, XMLStreamException, FileNotFoundException {		
		StatisticsFactory f = new StatisticsFactory(input, filter);
		if(env().numericPropertyExist("LowerSnp")) f.setLowerSnpFilter(env().integerProperty("LowerSnp"));
		if(env().numericPropertyExist("UpperSnp")) f.setUpperSnpFilter(env().integerProperty("UpperSnp"));		
		return f;
	}
	
	public Statistics createEmptyStatistics(TextInput input) 
		throws NaturalSetException, IOException {
		
		Statistics result = new Statistics(env(), env().stringProperty("Name"), new HaplotypeReader(env(), input));
		result = Statistics.clip(result, filterSNP(result.snpDomain()));
		return result;
	}

	
	public NaturalSet filterSNP(NaturalDomain domain){
		return new NaturalSet( domain, 
				env().numericPropertyExist("LowerSnp") ? env().numericProperty("LowerSnp").intValue() : domain.min(),
				env().numericPropertyExist("UpperSnp") ? env().numericProperty("UpperSnp").intValue() : domain.max());		
	}

	public NaturalSet filterARG(NaturalDomain domain){
		return new NaturalSet( domain, 
				env().numericPropertyExist("LowerArg") ? env().numericProperty("LowerArg").intValue() : domain.min(),
				env().numericPropertyExist("UpperArg") ? env().numericProperty("UpperArg").intValue() : domain.max());		
	}

}
