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
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;

import sanger.argml.io.TextOutput;

public class LogOutput extends TextOutput{
	private final static long DAY = 86400000;
	private final static long HOUR = 3600000;
	private final static long MINUTE = 60000;
	private final static long SECOND = 1000;
	private long clock;
	private long totalClock;
	
	public LogOutput(Environment env, OutputStream out){
		super(env, out);
		open();
	}
	
	public LogOutput(Environment env, String filename) throws IOException{
		super(env);
		file = new File(env.outbase(), filename);
		if(file.getParentFile().exists() ? true : file.getParentFile().mkdirs()) {
			fileWriter = new FileWriter(file, true);
			printWriter = new PrintWriter(fileWriter, true);				
		}
		open();
	}
	
	private void open(){
		totalClock = System.currentTimeMillis();
		clock = System.currentTimeMillis();
	}

	public void close() throws IOException{
		super.close();
	}
		
	public void println(String line){
		printWriter.println(line);
	}
	
	public void print(String line){
		printWriter.print(line);
	}
	
	public void reset(){
		clock = System.currentTimeMillis();
	}
	
	public void printBenchmark(String comment){
		long interval = System.currentTimeMillis() - clock;
		printWriter.println("INFO\t[" + normalizeTime(interval) + "]\t" + comment );
		clock = System.currentTimeMillis();
	}
	
	public void printTotalBenchmark(){
		long interval = System.currentTimeMillis() - totalClock;
		printWriter.println("INFO\t[" + normalizeTime(interval) + "]\ttotal runtime" );
	}
	
	public PrintWriter writer() {
		return printWriter;
	}
	
	public void printError(String message){
		printWriter.println("ERROR\t" + message);
	}

	public void printError(Exception e){
		printWriter.println("ERROR\t" + e.getMessage());
		e.printStackTrace(printWriter);
	}
	
	public void printInfo(String message){
		printWriter.println("INFO\t" + message);
	}
	
	private String normalizeTime(long ms){
		StringBuilder sb = new StringBuilder();
		String time;
		
		time = Integer.toString((int)Math.floor(ms/DAY));
		if(time.length() == 1) { sb.append('0'); }   
		sb.append(time);
		ms%=DAY;
		
		sb.append(':');
		time = Integer.toString((int)Math.floor(ms/HOUR));
		if(time.length() == 1) { sb.append('0'); }  
		sb.append(time);
		ms%=HOUR;

		sb.append(':');
		time = Integer.toString((int)Math.floor(ms/MINUTE));
		if(time.length() == 1) { sb.append('0'); }  
		sb.append(time);
		ms%=MINUTE;

		sb.append(':');
		time = Integer.toString((int)Math.floor(ms/SECOND));
		if(time.length() == 1) { sb.append('0'); }  
		sb.append(time);
		ms%=SECOND;

		sb.append('.');		
		sb.append(ms);
		
		return sb.toString();
	}

	public void write(){}
}
