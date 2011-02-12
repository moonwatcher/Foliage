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
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.regex.Pattern;

import javax.imageio.stream.ImageOutputStream;

import org.apache.oro.io.Perl5FilenameFilter;

import sanger.argml.drawing.RSquarePlotPainter;
import sanger.argml.drawing.TreeTraversalDistancePainter;
import sanger.argml.format.ArgmlDocument;
import sanger.argml.format.xml.GraphMLOutput;
import sanger.argml.graph.model.CoordinateTranslator;
import sanger.argml.graph.model.Genealogy;
import sanger.argml.graph.model.GenealogyReader;
import sanger.argml.graph.model.HaplotypeSet;
import sanger.argml.graph.model.Statistics;
import sanger.argml.graph.model.StatisticsFactory;
import sanger.argml.io.TextInput;
import sanger.argml.io.TextOutput;
import sanger.argml.io.XmlInput;
import sanger.argml.io.XmlOutput;
import sanger.argml.tools.GenealogyValidator;
import sanger.margarita.ArgBuilderForUnphasedData;
import sanger.margarita.InputParser;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSet;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;

public class ArgmlProcessor {
	protected static final DecimalFormat doubleToString = new DecimalFormat("0.###E0");

	private static String formatDouble(double x){
		String result = "0";
		if(x!=0.0){
			if(	x == Double.POSITIVE_INFINITY || x == Double.NEGATIVE_INFINITY ) {
				result = String.valueOf(x);
			} else {
				result = doubleToString.format(x);
			}
		}
		return result;
	}

	private static String normalizeCommand(String[] args){
		StringBuilder sb = new StringBuilder();
		sb.append(' ');
		for(int i=0; i<args.length; i++) { sb.append(args[i]); sb.append(' '); }
		return sb.toString();
	}
	
	public static void main(String[] args) {
		// Force foliage to use wstx implementation of stax
		System.setProperty("javax.xml.stream.XMLInputFactory", "com.ctc.wstx.stax.WstxInputFactory");
		System.setProperty("javax.xml.stream.XMLOutputFactory", "com.ctc.wstx.stax.WstxOutputFactory");
		System.setProperty("javax.xml.stream.XMLEventFactory", "com.ctc.wstx.stax.WstxEventFactory");
		
		String command = normalizeCommand(args);
		ProcessManager p = null;
		try {
			p = new ProcessManager(command);
			process(p);
			
		} catch (Exception e) {
			System.out.println(e.getMessage());
			
		} finally {
			try { if(p != null) { p.env().close(); } } 
			catch(IOException e){
				System.out.println("Failed to close log file.");
				e.printStackTrace();	
			}
		}			
	}

	public static void process(ProcessManager p) {
		try {
			help(p);
			haplotypeFromPhase(p);
			permuteHaplotypes(p);
			validatingGenealogy(p);
			validatingSubGenealogy(p);
			calculateStatistics(p);
			summarizeStatistics(p);
			pairwiseTreeTraversalDistance(p);
			pairwiseTreeTraversalDistanceP(p);
			ld(p);
			mixedPlot(p);
			graphML(p);
			debug(p);
			dotOutput(p);
			margaritaArgOutput(p);
			collectStatistics(p);
			filterStatistics(p);

		} catch (Exception e) {
			p.env().log().printError(e);
		}
	}	
	
	private static void help(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("help")){
			TextOutput out = null;
			
			try {
				out = p.textOutput();
				if(p.env().stringPropertyExist("Topic")){
					p.env().printTopic(out, p.env().stringProperty("Topic"));					
				} else {
					p.env().printCommands(out);					
				}
				
			} finally {
				if(out!=null) out.close();
			}			
		}
	}
	
	private static void dotOutput(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("dot")){
			TextInput in = null;
			TextOutput out = null;
			
			try {
				in = p.textInput();
				out = p.textOutput();

				GenealogyReader f = p.createGenealogyFactory(in);				
				for(Genealogy g : f) {
					if(p.env().integerProperty("Type") == 1){
						g.printDot1(out);
					} else {
						g.printDot2(out);
					}
					p.env().log().printBenchmark("writing genealogy " + f.index());
				}
				
			} finally {
				if(in!=null) in.close();
				if(out!=null) out.close();
			}
		}
	}
	
	private static void debug(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("debug")){
			TextInput in = null;
			TextOutput out = null;
			
			try {
				in = p.textInput();
				out = p.textOutput();

				GenealogyReader gi = p.createGenealogyFactory(in);				
				for(Genealogy g : gi) {
					out.writer().println("ARG " + gi.index());
					g.print(out);
					out.writer().println();
					p.env().log().printBenchmark("writing genealogy " + gi.index());
				}
				
			} finally {
				if(in!=null) in.close();
				if(out!=null) out.close();
			}
		}
	}
	
	private static void graphML(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("graphml")){
			TextInput in = null;
			XmlOutput out = null;
			
			try {
				in = p.textInput();
				out = p.xmlOutput();
				
				GenealogyReader f = p.createGenealogyFactory(in);
				GraphMLOutput gml = new GraphMLOutput(p.env(), out);

				gml.writeStart();
				for(Genealogy g : f) {					
					gml.write(g);
					p.env().log().printBenchmark("writing genealogy " + f.index());
				}				
				gml.writeEnd();
				
			} finally {
				if(in!=null) in.close();
				if(out!=null) out.close();
			}
		}
	}
	
	private static void haplotypeFromPhase(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("from-phase")){
			TextInput fpin = null;
			TextInput fpout = null;
			TextOutput out = null;
			
			try {
				HaplotypeSet hs = null;
				fpin = new TextInput(p.env(), p.env().stringProperty("FastPHASEInput"));
				
				if(p.env().stringPropertyExist("FastPHASEOutput")){
					fpout = new TextInput(p.env(), p.env().stringProperty("FastPHASEOutput"));
					hs = HaplotypeSet.readFastPHASE(p.env(), fpin, fpout);
				} else {
					hs = HaplotypeSet.readFastPHASE(p.env(), fpin);
				}

				out = p.textOutput();
				hs.write(out);
				p.env().log().printBenchmark("writing haplotypes file " + out);					
				
			} finally { 
				if(fpin!=null) fpin.close(); 
				if(fpout!=null) fpout.close(); 
				if(out!=null) out.close(); 
			}			
		}
	}

	private static void permuteHaplotypes(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("subsample-haplotype")){
			TextInput in = null;
			TextOutput out = null;

			try {
				in = p.textInput();
				HaplotypeSet hs = p.createHaplotypeSet(in);
				
				int subsample = p.env().numericProperty("SubSampleSize").intValue();
				int sample =  p.env().numericProperty("SampleSize").intValue();

				for (int i=1; i<=sample; i++){
					NaturalDomain rd = hs.haplotypeDomain().createRandomSubDomain(subsample);
					HaplotypeSet rhs = HaplotypeSet.clipHaplotypeDomain(hs, rd);
					
					out = new TextOutput(p.env(), i + "_" + p.env().stringProperty("Output"));
					rhs.write(out);
					p.env().log().printBenchmark("writing haplotype sub sample file " + out + " : " + rhs.haplotypeDomain());
					out.close();
					
					if(p.env().flag("NonOverlappingSubSample")){
						NaturalDomain ord = hs.haplotypeDomain().createRandomSubDomain(rd, subsample);
						HaplotypeSet orhs = HaplotypeSet.clipHaplotypeDomain(hs, ord);
						
						out = new TextOutput(p.env(), i + "_f_" + p.env().stringProperty("Output"));
						orhs.write(out);
						p.env().log().printBenchmark("writing non overlapping haplotype sub sample file " + out + " : " + orhs.haplotypeDomain());
						out.close();
					}
				}									
			} finally { 
				if(out!=null) out.close(); 
				if(in!=null) in.close(); 
			}			
		}
	}	
	
	private static void pairwiseTreeTraversalDistance(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("pttd")){
			XmlInput in = null;
			ImageOutputStream out = null;
			int order = 0;

			try {
				in = p.xmlInput();
				out = p.imageOutputStream();
				Statistics one = null, two = null;
				TreeTraversalDistancePainter pttd = null;
				
				StatisticsFactory f = p.createStatisticsFactory(in);
			
				for(Statistics s : f){
					if(s.name().equals(p.env().stringProperty("Name"))){
						one = s;
					}
					
					if(p.env().stringPropertyExist("Other") && s.name().equals(p.env().stringProperty("Other"))){ 
						two = s;
					}
				}
				
				
				if(one!=null){
					pttd = new TreeTraversalDistancePainter(p.env(), one);
					pttd.paintDiagram();
					if(p.env().flag("RecombinationDensityDistribution")) pttd.paintRecombinationDensity(order++);					
					if(p.env().flag("SnpDensityDistribution")) pttd.paintSnpDensity(order++);
					pttd.paintLegend();
				}
				
				if(two!=null){
					TreeTraversalDistancePainter opttd = new TreeTraversalDistancePainter(p.env(), pttd, two);
					opttd.paintCompareDiagram();
					if(p.env().flag("RecombinationDensityDistribution")) opttd.paintRecombinationDensity(order++);
					if(p.env().flag("SnpDensityDistribution")) opttd.paintSnpDensity(order++);
				}
				
				pttd.write(out);
								
			} finally { 
				if(in!=null) in.close();
				if(out!=null) out.close();
			}
		}
	}	

	private static void pairwiseTreeTraversalDistanceP(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("pttdp")){
			XmlInput in = null;	
			ImageOutputStream out = null;
			
			try {
				Normal normal = new Normal(0.0, 1.0, new DRand());

				in = p.xmlInput();
				int tile = p.env().integerProperty("Tile");
				Statistics one = null, two = null;
				ArrayList<double[]> population = new ArrayList<double[]>();
				double[] x = null;
				
				
				StatisticsFactory factory = p.createStatisticsFactory(in, null);//p.env().stringProperty("Name"));
				Pattern outgroupFilter = Pattern.compile(p.env().stringProperty("Other"));
	
				for(Statistics s : factory){
					if(outgroupFilter.matcher(s.name()).matches()) {
						if(two == null) {
							two = s;
							
						} else {
							two.addStatistics(s);
						}
						p.env().log().printBenchmark("outgroup: " + s);

					} else {
						if(one == null) {
							one = s;
							
						} else {
							one.addStatistics(s);
							CoordinateTranslator translate = new CoordinateTranslator(p.env(), s.snpDomain(), s.basePairDomain(), s.markerPositions());					
							population.add(translate.tiledDensity(tile, s.recombinationRates()));						
						}
						p.env().log().printBenchmark("population: " + s);
					}
				}
				
				CoordinateTranslator translate = new CoordinateTranslator(p.env(), two.snpDomain(), two.basePairDomain(), two.markerPositions());					
				x  = translate.tiledDensity(tile, two.recombinationRates());

				int regions = population.get(0).length;
				int populationSize = population.size();
	
				double[] meanOfSqures = new double[regions];

				double[] mean = new double[regions];
				double[] sd = new double[regions];
				double[] zvalue = new double[regions];
				double[] pvalue = new double[regions];
				int[] qvalue = new int[regions];
				int[] bqvalue = new int[regions];
				boolean[] significant = new boolean[regions];
				
				for(int i=0; i<regions; i++) {
					//	Sum all values
					for(int j=0; j<population.size(); j++) {
						mean[i] += population.get(j)[i];
						meanOfSqures[i] += Math.pow(population.get(j)[i], 2);
					}
					
					//	normalize values
					mean[i] /= populationSize;
					meanOfSqures[i] /= populationSize;
					sd[i] = Math.sqrt(meanOfSqures[i] - Math.pow(mean[i], 2));
				}
				
				for(int i=0; i<regions; i++) {
					qvalue[i] = 0;
					if(x[i] > mean[i]){
						for(int j=0; j<population.size(); j++) {
							if(population.get(j)[i] >= x[i]) qvalue[i]++;
						}	
					} else {
						for(int j=0; j<population.size(); j++) {
							if(population.get(j)[i] <= x[i]) qvalue[i]++;
						}						
					}
					
					bqvalue[i] = 0;
					double otherside = 2 * mean[i] - x[i];
					if(x[i] > mean[i]){
						for(int j=0; j<population.size(); j++) {
							if(population.get(j)[i] <= otherside) bqvalue[i]++;
						}	
					} else {
						for(int j=0; j<population.size(); j++) {
							if(population.get(j)[i] >= otherside) bqvalue[i]++;
						}	
					}
					
				}
				
				for(int i=0; i<regions; i++){
					if(sd[i] == 0){
						significant[i] = (mean[i] != x[i]);
						zvalue[i] = Double.POSITIVE_INFINITY;
						
						if(significant[i]) 
							pvalue[i] = 1.0 / (double)populationSize;
						else 
							pvalue[i] = 0.0;
								
					} else {
						zvalue[i] = ( x[i] - mean[i] ) / sd[i];
						double cdf = normal.cdf(zvalue[i]);
						pvalue[i] = Math.min(cdf, 1.0 - cdf);
						
						//	Significance testing
						//	Setting who is significant
						//	instead of checking the P-Value you can check the q or bq here
						significant[i] = pvalue[i] < 0.01;						
					}

					StringBuilder sb = new StringBuilder();
					sb.append("tile: ");
					sb.append(i);
					sb.append(" z: ");
					sb.append(formatDouble(zvalue[i]));
					sb.append(" p: ");
					sb.append(formatDouble(pvalue[i]));
					sb.append(" q: ");
					sb.append(qvalue[i]);
					sb.append(" bq: ");
					sb.append(bqvalue[i]);

					sb.append(" x: ");
					sb.append(formatDouble(x[i]));
					sb.append(" m: ");
					sb.append(formatDouble(mean[i]));
					sb.append(" sd: ");
					sb.append(formatDouble(sd[i]));
					
					p.env().log().printInfo(sb.toString());
				}

				out = p.imageOutputStream();
				TreeTraversalDistancePainter pttd = null;
				int order = 0;
				
				if(one!=null){
					pttd = new TreeTraversalDistancePainter(p.env(), one);
					pttd.paintDiagram();
					if(p.env().flag("RecombinationDensityDistribution")) pttd.paintRecombinationDensity(order++);					
					if(p.env().flag("SnpDensityDistribution")) pttd.paintSnpDensity(order++);
					pttd.paintLegend();
				}
				
				if(two!=null){
					TreeTraversalDistancePainter opttd = new TreeTraversalDistancePainter(p.env(), pttd, two);
					opttd.paintCompareDiagram();
					if(p.env().flag("RecombinationDensityDistribution")) opttd.paintRecombinationDensity(order++, significant);
					if(p.env().flag("SnpDensityDistribution")) opttd.paintSnpDensity(order++);
				}
				
				
				pttd.write(out);
				p.env().log().printBenchmark("write " + out);
								
			} finally {
				if(out!=null) out.close();
				if(in!=null) in.close();
			}
		}
	}	
		
	private static void ld(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("ld")){
			TextInput in = null;
			TextInput oin = null;
			ImageOutputStream out = null;
			int order = 0;

			try {
				in = p.textInput();
				HaplotypeSet hs = p.createHaplotypeSet(in);
				out = p.imageOutputStream();
				RSquarePlotPainter ldp = new RSquarePlotPainter(p.env(), hs);
				ldp.paintDiagram();
				if(p.env().flag("MinorAlleleDensityDistribution")) ldp.paintMinorAlleleDensity(order++);
				if(p.env().flag("SnpDensityDistribution")) ldp.paintSnpDensity(order++);

				ldp.drawLegend();
				
				if(p.env().stringPropertyExist("Other")) {
					oin = new TextInput(p.env(), p.env().stringProperty("Other"));
					HaplotypeSet ohs = p.createHaplotypeSet(oin);
					RSquarePlotPainter oldp = new RSquarePlotPainter(p.env(), ldp, ohs);
					oldp.paintCompareDiagram();
					if(p.env().flag("MinorAlleleDensityDistribution")) oldp.paintMinorAlleleDensity(order++);
					if(p.env().flag("SnpDensityDistribution")) oldp.paintSnpDensity(order++);

				}
				
				ldp.write(out);
								
			} finally { 
				if(in!=null) in.close();
				if(oin!=null) in.close();
				if(out!=null) out.close();
			}
		}
	}	

	private static void mixedPlot(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("pttd-ld")){
			XmlInput pttdin = null;
			TextInput ldin = null;
			ImageOutputStream out = null;
			int order = 0;

			try {
				ldin = new TextInput(p.env(), p.env().stringProperty("Other"));
				HaplotypeSet hs = p.createHaplotypeSet(ldin);
				
				pttdin = p.xmlInput();
				out = p.imageOutputStream();
				Statistics one = null;
				
				StatisticsFactory f = p.createStatisticsFactory(pttdin);
			
				for(Statistics s : f){
					if(s.name().equals(p.env().stringProperty("Name"))) 
						one = s;
				}
				
				TreeTraversalDistancePainter pttd = new TreeTraversalDistancePainter(p.env(), one);
				pttd.paintDiagram();
				if(p.env().flag("RecombinationDensityDistribution")) pttd.paintRecombinationDensity(order++);
				if(p.env().flag("SnpDensityDistribution")) pttd.paintSnpDensity(order++);
				pttd.paintLegend();
				RSquarePlotPainter ldp = new RSquarePlotPainter(p.env(), pttd, hs);
				ldp.paintCompareDiagram();
				if(p.env().flag("MinorAlleleDensityDistribution")) ldp.paintMinorAlleleDensity(order++);
				if(p.env().flag("SnpDensityDistribution")) ldp.paintSnpDensity(order++);

				pttd.write(out);
								
			} finally { 
				if(ldin!=null) ldin.close();
				if(pttdin!=null) pttdin.close();
				if(out!=null) out.close();
			}
		}
	}	
	
	private static void summarizeStatistics(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("ssummarize")){
			XmlInput in = null;
			XmlOutput out = null;
			ArgmlDocument doc = null;
			try {
				Statistics summary = null;
				in = p.xmlInput();
				StatisticsFactory f = p.createStatisticsFactory(in, p.env().stringProperty("Pattern"));
				out = p.xmlOutput();
				doc = new ArgmlDocument(out);
				doc.writeStartDocument();
				for(Statistics s : f){
					if(summary==null) summary = s;
					else summary.addStatistics(s);
					p.env().log().printBenchmark("adding " + s);
				}
				summary.setName(p.env().stringProperty("Name"));
				summary.writeElement(out);
				doc.writeEndDocument();
				p.env().log().printBenchmark("writing " + out);
			} finally {
				if(in!=null) in.close();
				if(out!=null) out.close();
			}			
		}		
	}	
	
	private static void filterStatistics(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("sfilter")){
			XmlInput in = null;
			XmlOutput out = null;
			ArgmlDocument doc = null;
			try {
				in = p.xmlInput();
				StatisticsFactory f = p.createStatisticsFactory(in, p.env().stringProperty("Pattern"));
				out = p.xmlOutput();
				doc = new ArgmlDocument(out);
				doc.writeStartDocument();
				for(Statistics s : f){
					s.writeElement(out);
					p.env().log().printBenchmark("adding " + s);
				}				
				in.close();
				
				if(p.env().stringPropertyExist("Other")){
					in = new XmlInput(p.env(), p.env().stringProperty("Other"));
					f = p.createStatisticsFactory(in, p.env().stringProperty("Pattern"));
					for(Statistics s : f){
						s.writeElement(out);
						p.env().log().printBenchmark("adding " + s);
					}				
				}
				doc.writeEndDocument();
				p.env().log().printBenchmark("writing " + out);
			} finally {
				if(in!=null) in.close();
				if(out!=null) out.close();
			}			
		}		
	}	
	
	private static void collectStatistics(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("scollect")){
			XmlInput in = null;
			XmlOutput out = null;
			ArgmlDocument doc = null;
			try {
				out = p.xmlOutput();
				doc = new ArgmlDocument(out);
				doc.writeStartDocument();

				FilenameFilter filter = new Perl5FilenameFilter(p.env().stringProperty("Pattern"));
				File[] files = p.env().inbase().listFiles(filter);
				for(File file : files){
					in = new XmlInput(p.env(), file.getName());
					StatisticsFactory f = p.createStatisticsFactory(in);					
					for(Statistics s : f){
						s.writeElement(out);
						p.env().log().printBenchmark("adding " + s);
					}				
					in.close();
				}				
				doc.writeEndDocument();
				p.env().log().printBenchmark("writing " + out);				
			} finally {
				if(in!=null) in.close();
				if(out!=null) out.close();
			}			
		}		
	}
			
	private static void calculateStatistics(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("statistic")){
			TextInput in = null;
			TextInput hf = null;
			XmlOutput out = null;
			
			try {
				hf = new TextInput(p.env(), p.env().stringProperty("Haplotypes"));
				in = p.textInput();
				out = p.xmlOutput();
				
				Statistics s = p.createEmptyStatistics(hf);
				GenealogyReader f = p.createGenealogyFactory(in);				
				p.env().log().printBenchmark("initialized " + (p.env().flag("MultiFurcate") ? "multifurcating " : "bifurcating") + " genealogy reader");
				
				for(Genealogy g : f) {				
					p.env().log().printBenchmark("reading genealogy " +  f.index() + " : " + g);
					s.addGenealogy(g);
					p.env().log().printBenchmark("calculating statistics for genealogy " +  f.index());						
				}

				s.writeDocument(out);
				p.env().log().printBenchmark("writing " + out);
				
			} finally {
				if(in!=null) in.close();
				if(hf!=null) hf.close();
				if(out!=null) out.close();
			}			
		}
	}
	
	private static void validatingGenealogy(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("validate")){
			TextInput in = null;
			TextInput hf = null;

			try {
				in = p.textInput();
				hf = new TextInput(p.env(), p.env().stringProperty("Haplotypes"));
				
				GenealogyReader f = p.createGenealogyFactory(in);
				HaplotypeSet hs = p.createHaplotypeSet(hf);			
				p.env().log().printBenchmark("initialized " + (p.env().flag("MultiFurcate") ? " multi furcated " : "") + "genealogy factory");
				
				for(Genealogy g : f) {					
					p.env().log().printBenchmark("build: " + g.toString());

					if(p.env().flag("ReportEmptyRegions")) {
						g.reportEmptyRegions();
						p.env().log().printBenchmark("check for edges with empty active region");
					}

					GenealogyValidator v = new GenealogyValidator(p.env(), g, hs);
					v.print();
					p.env().log().printBenchmark("validated " + f.index());
				}

			} finally {
				if(in!=null) in.close();
				if(hf!=null) hf.close();
			}
		}
	}
	
	private static void validatingSubGenealogy(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("svalidate")){
			TextInput in = null;
			TextInput hf = null;

			try {
				in = p.textInput();
				hf = new TextInput(p.env(), p.env().stringProperty("Haplotypes"));

				GenealogyReader f = p.createGenealogyFactory(in);
				HaplotypeSet hs = p.createHaplotypeSet(hf);				
				p.env().log().printBenchmark("initialized " + (p.env().flag("MultiFurcate") ? " multi furcated " : "") + "genealogy factory");
				
				for(Genealogy g : f) {
					p.env().log().printBenchmark("build: " + g.toString());
					
					for(int i=g.snpDomain().min(); i<=g.snpDomain().max(); i++){
						for(int j=i; j<=g.snpDomain().max(); j++){
							NaturalSet zone = new NaturalSet(g.snpDomain(), i, j);
							Genealogy sg = Genealogy.clipSnpDomain(g, zone);
							
							if(p.env().flag("ReportEmptyRegions")) {
								sg.reportEmptyRegions();
								p.env().log().printBenchmark("check for edges with empty active region");
							}
							
							HaplotypeSet shs = HaplotypeSet.clipSnpDomain(hs, zone);					
							GenealogyValidator v = new GenealogyValidator(p.env(), sg, shs);
							v.print();
							p.env().log().printBenchmark("validation");
						}
					}
				}

			} finally {
				if(in!=null) in.close();
				if(hf!=null) hf.close();
			}
		}
	}

	private static void margaritaArgOutput(ProcessManager p) throws Exception{
		if(p.env().instruction().getName().equals("margarita")){
			TextInput in = null;
			TextOutput out = null;
			
			try {
				in = p.textInput();
				out = p.textOutput();
				
				out.writer().println("Margarita 250707");
				final InputParser ip = new InputParser();
				ip.parseFile(p.env(), in);
				
				final ArgBuilderForUnphasedData ab = new ArgBuilderForUnphasedData(p.env(), out);
				ab.buildArgs((int)Math.floor(p.env().numericProperty("SampleSize")),ip);
				ab.printArgs();
				
			} finally {
				if(in!=null) in.close();
				if(out!=null) out.close();
			}			
		}
	}
}
