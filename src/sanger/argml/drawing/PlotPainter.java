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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

import javax.imageio.stream.ImageOutputStream;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;
import sanger.argml.io.ImageOutput;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSetException;

/**
 * Abstract implementation of a heat plot painter.
 * <p>Can be used to plot any pairwise snp statistic scaled to chromosome coordinates.
 * The painter facilitates comparing two data sets or two statistics for the same data set 
 * as long as they are defined over the same snp and chromosome coordinates. 
 * The main observation here is that every pixel in the final output covers a fixed number 
 * of bases in chromosome coordinates given by the <code>tile</code> variable.</p>
 * <p>moving from the snp to the chromosome coordinate system is done by calculating 
 * an average of all pairwise correlation overlapping the pixel, weighted by the fraction 
 * of the pixel the corresponding snp frame is covering. The snp frame is defined as the interval, 
 * in chromosome coordinate, starting in the middle of the gap between a snp and the one preceding it 
 * and ending with in the middle of the gap between a snp and the one following. 
 * More formally, if <code>s(i)</code> is the chromosome coordinate of snp i, then its fram is: 
 * <code>[(s(i-1) + s(i)) / 2, (s(i) + s(i+1)) / 2]</code>.</p>
 * <p>Another facility provided by the painter is plotting distributions aligned with the same coordinate system,
 * and the snp density at each tile, which is useful for normalizing values by the snp density.</p>
 * @author Lior Galanti
 *
 */
/**
 * @author lg8
 *
 */
public abstract class PlotPainter extends Environmental{
	protected static final DecimalFormat doubleToString = new DecimalFormat("0.00");
	
	//protected CoordinateTranslator translate; 
	
	protected BufferedImage image;
	protected Graphics2D stage;
		
	protected JetColorMap colorMap;
	protected double tile;
	protected Integer gridSpacing;
	protected int distributionRange;
	protected int axisRange;
	protected int padding;
	protected int scaleRange;
	
	protected ArrayList<Integer> customgrid;

	protected Color lineColor;
	protected Color backgroundColor;
	protected Color gridColor;
	protected Color notchColor;
	protected Color legendColor;
	protected Color positiveTileColor;
	protected Color negativeTileColor;
	protected Color barColor;
	protected Color significantBarColor;
	protected Color meanColor;

	protected String font;
	protected int fontsize;
	
	protected int width;
	protected int height;
	protected int xPx = 0;
	protected int yPx = 0;
	
	protected ArrayList<Frame> frames;
	protected ArrayList<Tile> tiles;

	
	/**
	 * Implemented in a derived class to return the snp domain.
	 * @return
	 */
	protected abstract NaturalDomain snpDomain();

	/**
	 * Implemented in a derived class to return the chromosome domain.
	 * @return
	 */
	protected abstract NaturalDomain basePairDomain();
	
	/**
	 * Implemented in a derived class to return the chromosome coordinate of a snp.
	 * @param i the snp coordinate.
	 * @return snp i's chromosome coordinate.
	 * @throws NaturalSetException
	 */
	protected abstract int coordinate(int i) throws NaturalSetException;
		
	/**
	 * Implemented in a derived class to return the pairwise correlation statistic of the two snps.
	 * @param i snp corrdinate for the first snp
	 * @param j snp coordinate for the second snp
	 * @return the pairwise correlation statistic of the two snps.
	 * @throws NaturalSetException
	 */
	protected abstract double pairwiseSnpCorrelation(int i, int j) throws NaturalSetException;
	
	public PlotPainter(Environment env) throws NaturalSetException {
		super(env);
	}

	/**
	 * Used in mixed mode plotters to initialize a new plotter that will draw on another plotter's canvas.
	 * @param env
	 * @param base
	 * @throws NaturalSetException
	 */
	public PlotPainter(Environment env, PlotPainter base) throws NaturalSetException {
		super(env);
		this.image = base.image;
		this.stage = base.stage;
	}
	
	protected int distributions() {
		int result = 0;
		if(env().flag("MinorAlleleDensityDistribution")) result++;
		if(env().flag("RecombinationDensityDistribution")) result++;
		if(env().flag("SnpDensityDistribution")) result++;
			
		if(env().stringPropertyExist("Other")) {
			if(env().flag("MinorAlleleDensityDistribution")) result++;
			if(env().flag("RecombinationDensityDistribution")) result++;
			if(env().flag("SnpDensityDistribution")) result++;			
		}
		
		return result;
	}
			
	protected void initialize() throws NaturalSetException{
		readProperties();
		buildTiles();		
	}
	
	protected void readProperties(){
		this.colorMap = new JetColorMap(env());
		
		this.tile = env().numericProperty("Tile");
		this.gridSpacing = env().integerProperty("GridSpacing");		
		this.distributionRange = env().integerProperty("DistributionRange");
		this.axisRange = env().integerProperty("AxisRange");
		this.padding = env().integerProperty("Padding");
		this.scaleRange = env().integerProperty("ScaleRange");
		
		xPx = padding;
		yPx = padding;
		
		// Colors
		this.lineColor = new Color(Long.decode(env().stringProperty("LineColor")).intValue(), true);
		this.backgroundColor = new Color(Long.decode(env().stringProperty("BackgroundColor")).intValue(), true);
		this.gridColor = new Color(Long.decode(env().stringProperty("GridColor")).intValue(), true);
		this.notchColor = new Color(Long.decode(env().stringProperty("NotchColor")).intValue(), true);
		this.legendColor = new Color(Long.decode(env().stringProperty("LegendColor")).intValue(), true);
		this.positiveTileColor = new Color(Long.decode(env().stringProperty("PositiveTileColor")).intValue(), true);
		this.negativeTileColor = new Color(Long.decode(env().stringProperty("NegativeTileColor")).intValue(), true);
		this.barColor = new Color(Long.decode(env().stringProperty("BarColor")).intValue(), true);		
		this.significantBarColor = new Color(Long.decode(env().stringProperty("SignificantBarColor")).intValue(), true);		
		this.meanColor = new Color(Long.decode(env().stringProperty("MeanColor")).intValue(), true);

		this.font = env().stringProperty("Font");
		this.fontsize = env().integerProperty("FontSize");
		
		if(env().stringPropertyExist("CustomGridLine")){
			String[] list = env().stringProperty("CustomGridLine").split(",");
			customgrid = new ArrayList<Integer>(list.length);
			for(int i=0; i<list.length; i++){
				customgrid.add(Integer.parseInt(list[i]));
			}
			
			Collections.sort(customgrid);
		}
	}
	
	protected void buildTiles() throws NaturalSetException {
		int min = snpDomain().min();
		int max = snpDomain().max();
		
		frames = new ArrayList<Frame>(snpDomain().closureCardinality());
		tiles = new ArrayList<Tile>((int)Math.ceil(basePairDomain().closureCardinality() / tile));
		
		//	build frames
		frames.add(new Frame(min, coordinate(min), coordinate(min+1)));
		
		for(int i=min+1; i<max;i++){
			frames.add(new Frame(i, coordinate(i-1), coordinate(i), coordinate(i+1))); 
		}
		
		frames.add(new Frame(max, coordinate(max), coordinate(max-1)));
		
		//	build tiles
		int count = 0;
		double mark = basePairDomain().min();
		Frame left = frames.get(0);
		
		for(Frame frame : frames){
			
			//	while the mark is still to the left of the mark we need to add another tile
			while(frame.bpright > mark + tile){
				tiles.add(new Tile(count, left, frame, mark, mark + tile));
				count++;
				mark += tile;				
				if(mark > left.bpright) left = frame;
			}
		}
		
		// Complete the last, incomplete, tile
		if(mark < basePairDomain().max()){
			tiles.add(new Tile(count, left, frames.get(frames.size()-1), mark, basePairDomain().max()));
		}
	}	
	
	protected void buildStage(){
		width = xPx + tiles.size() + 4 * padding + 2 * axisRange + scaleRange + 30 * fontsize;
		height = yPx + tiles.size() + distributions() * (distributionRange + padding ) + 2 * padding;
		
		image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		stage = image.createGraphics();
		stage.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		stage.setFont(new Font(font, Font.PLAIN, fontsize));
		Graphics2D background =  (Graphics2D)stage.create();
		background.setColor(backgroundColor);
		background.fillRect(0, 0, image.getWidth(), image.getHeight());
	}
			
	public void paintDiagram() throws NaturalSetException {
		if(stage == null) buildStage();
		paintPlot(xPx, yPx);
		if(gridSpacing != null) { 
			paintGrid(xPx, yPx);
			paintMarkerAxis(tiles.size() + padding + xPx, yPx);
			paintSnpAxis(tiles.size() + 3 * padding + axisRange + xPx, yPx);
		}
		paintHeatScale(xPx + tiles.size() + axisRange * 2 + 4 * padding, yPx);
	}

	
	public void paintCompareDiagram() throws NaturalSetException{
		if(stage == null) buildStage();
		AffineTransform t = AffineTransform.getScaleInstance(-1,1);
		t.concatenate(AffineTransform.getTranslateInstance(-tiles.size(), 0));
		t.concatenate(AffineTransform.getQuadrantRotateInstance(1, tiles.size() / 2.0, tiles.size() / 2.0));
		stage.setTransform(t);

		paintPlot(yPx, xPx);
		if(gridSpacing != null) { 
			paintGrid(yPx, xPx);
		}
		stage.setTransform(AffineTransform.getScaleInstance(1, 1));
	}
	
	public void paintDistribution(int order, double[] distribution, boolean mean) throws NaturalSetException{
		paintDistribution(order, distribution, mean, null);
	}

	public void paintDistribution(int order, double[] distribution, boolean mean, boolean[] significant) throws NaturalSetException{
		paintDistribution(xPx, yPx + tiles.size() + (order * (distributionRange + padding)) + padding, distribution, mean, significant);
	}

	public void paintDistribution(int x, int y, double[] distribution, boolean mean) throws NaturalSetException{
		paintDistribution(x, y, distribution, mean, null);
	}

	public void paintDistribution(int x, int y, double[] distribution, boolean mean, boolean[] significant) throws NaturalSetException{
		double vmax = 0.0;
		double vmean = 0.0;
		
		for(double v : distribution){
			vmax = Math.max(vmax, v);
			vmean += v;			
		}	vmean /= distribution.length;
		
		Graphics2D lines = (Graphics2D)stage.create();
		lines.translate(x, y);
		lines.setColor(lineColor);
				
		Graphics2D plot = (Graphics2D)stage.create();	
		plot.translate(x, y + 1);
		plot.setColor(barColor);

		Graphics2D values = (Graphics2D)stage.create();
		values.translate(x + tiles.size() + 4 * padding, y + stage.getFont().getSize() / 2);
		values.setColor(lineColor);
		
		lines.drawLine(0, 0, tiles.size() + 3 * padding, 0);
				
		lines.setColor(lineColor);
		lines.drawLine(tiles.size() + padding, 1, tiles.size() + padding, distributionRange);
		lines.drawLine(tiles.size() + padding, distributionRange, tiles.size() + 3 * padding, distributionRange);
		values.drawString(doubleToString.format(vmax) + " max", 0, distributionRange);
		
		if(mean){
			int meanLine = (int)Math.round(vmean/vmax*distributionRange);
			lines.setColor(meanColor);
			lines.drawLine(0, meanLine, tiles.size() + 3 * padding, meanLine);
			values.setColor(meanColor);
			values.drawString(doubleToString.format(vmean) + " mean", 0, meanLine);
		}
		
		lines.setColor(positiveTileColor);
		
		boolean odd = true;
		int block = 0;
		int i = 0;

		for(Tile p : tiles){
			if(significant!=null) {
				if(significant[i]) plot.setColor(significantBarColor);
				else plot.setColor(barColor);
			}
			lines.fillRect(p.position, 1, 1, distributionRange);
			int bar = (int)Math.round((distribution[i] / vmax) * (double)distributionRange);
			plot.fillRect(p.position, 0, 1, bar);
			block += p.width();
			if(block > gridSpacing){
				block=0;
				odd = odd ? false : true;
				lines.setColor(odd ? positiveTileColor : negativeTileColor);
			}
			i++;
		}
	}	
	
	public void paintPlot(int x, int y) throws NaturalSetException{
		if(stage == null) buildStage();
		Graphics2D canvas = (Graphics2D)stage.create();
		canvas.translate(x, y);
		
		for(Tile p : tiles){
			for(int j=0; j<=p.position; j++){
				Tile op = tiles.get(j);
				canvas.setColor(colorMap.mapToColor(p.correlation(op)));
				canvas.fillRect(p.position, op.position, 1, 1);
			}
		}
	}
			
	public void paintGrid(int x, int y) throws NaturalSetException{
		if(stage == null) buildStage();
		int block = 0;
		Graphics2D grid = (Graphics2D)stage.create();
		grid.translate(x, y);
		grid.setColor(gridColor);
		
		int i = 0;
		for(Tile p : tiles){
			block+=p.width();
			if(block > gridSpacing){
				block=0;
				grid.drawLine(p.position, p.position, tiles.size()-1, p.position);
				grid.drawLine(p.position, p.position, p.position, 0);
			}
			
			if(customgrid!=null && i<customgrid.size() && customgrid.get(i) < p.bpright && customgrid.get(i) > p.bpleft){
				grid.drawLine(p.position, p.position, tiles.size(), p.position);
				grid.drawLine(p.position, p.position, p.position, 0);
				i++;
			}
		}
	}
	
	public void paintSnpAxis(int x, int y) throws NaturalSetException{
		if(stage == null) buildStage();
		int block = 0;
		boolean even = true;

		Graphics2D lines = (Graphics2D)stage.create();
		lines.translate(x, y);
		lines.setColor(lineColor);
		
		Graphics2D values = (Graphics2D)stage.create();
		values.setColor(lineColor);
		values.translate(x + 3 * padding, y + stage.getFont().getSize() / 2);
		
		lines.drawLine(0, 0, 0, tiles.size() - 1);
		
		Tile first = tiles.get(0);
		Tile last = tiles.get(tiles.size() -1);
		
		lines.drawLine(1, first.position, 2 * padding, first.position);
		values.drawString(String.valueOf(first.minSnp()), 0, first.position);
		
		for(Tile p : tiles){
			block+=p.width();
			if(block > gridSpacing){
				block=0;
				even = even ? false : true;
				if(even) lines.drawLine(1, p.position, padding, p.position);
				else lines.drawLine(1, p.position, 2 * padding, p.position);
				values.drawString(String.valueOf(p.meanSnp()), 0, p.position);
			}
		}
		
		lines.drawLine(1, last.position, 2 * padding, last.position);
	}
	
	public void paintMarkerAxis(int x, int y) throws NaturalSetException{
		if(stage == null) buildStage();
		int block = 0;
		boolean even = true;

		Graphics2D lines = (Graphics2D)stage.create();
		lines.translate(x, y);
		lines.setColor(lineColor);
		
		Graphics2D values = (Graphics2D)stage.create();
		values.setColor(lineColor);
		values.translate(x + 3 * padding, y + stage.getFont().getSize() / 2);
		
		lines.drawLine(0, 0, 0, tiles.size() - 1);
		
		Tile first = tiles.get(0);
		Tile last = tiles.get(tiles.size() -1);
		
		lines.drawLine(1, first.position, 2 * padding, first.position);
		values.drawString(String.valueOf(first.minBp()), 0, first.position);
		
		int i =0;
		for(Tile p : tiles){
			block+=p.width();
			if(block > gridSpacing){
				block=0;
				even = even ? false : true;
				if(even) lines.drawLine(1, p.position, padding, p.position);
				else lines.drawLine(1, p.position, 2 * padding, p.position);
				values.drawString(String.valueOf(p.meanBp()), 0, p.position);
				
			}
			
			if(customgrid!=null && i<customgrid.size() && customgrid.get(i) < p.bpright && customgrid.get(i) > p.bpleft){
				lines.drawLine(1, p.position, 2 * padding, p.position);
				values.drawString(String.valueOf(p.meanBp()), 0, p.position);
				i++;
			}
		}
		
		lines.drawLine(1, last.position, 2 * padding, last.position);
//		values.drawString(String.valueOf(last.maxBp()), 0, last.position);
	}
	
	public void paintHeatScale(int x, int y){
		if(stage == null) buildStage();
		int slices = 256;
		int sliceSize = (int)Math.round( slices / colorMap.size());
		AffineTransform origin = AffineTransform.getTranslateInstance(x, y);
		
		// frame
		Graphics2D paint = (Graphics2D)stage.create();
		
		paint.setTransform(origin);
		paint.setColor(legendColor);
		paint.fillRect(0, 0, scaleRange + 4 * padding + stage.getFont().getSize() * 4, slices + 2 + 4 * padding);
		paint.setColor(notchColor);
		paint.fillRect(2 * padding, 2 * padding, scaleRange + 2, sliceSize * colorMap.size() + 2);
					
		// heat scale
		paint.setTransform(origin);
		paint.translate(2 * padding + 1, 2 * padding + 1);
		for(Color color : colorMap.scale()){
			paint.setColor(color);
			paint.fillRect(0, 0, scaleRange, sliceSize);
			paint.translate(0, sliceSize);
		}		

		// notch
		paint.setTransform(origin);
		paint.translate(2 * padding + scaleRange, sliceSize * 16 + 1 + 2 * padding);
		paint.setColor(notchColor);
		for(int i=16; i<colorMap.size(); i+=16){
			if(i % 32 == 0) paint.drawLine(0, 0, -(2 * padding), 0);
			else paint.drawLine(0, 0, -padding, 0);			
			paint.translate(0, sliceSize * 16);
		}

		// notch values
		paint.setTransform(origin);
		paint.setColor(notchColor);
		paint.translate(3 * padding + scaleRange, 2 * padding + paint.getFont().getSize() / 2);		

		for(int i=0; i<=colorMap.size(); i+=16){
			if(i % 32 == 0) paint.drawString(String.valueOf(((double)i / (double)colorMap.size())), 0, 0);			
			paint.translate(0, sliceSize * 16);
		}		
	}
	
	public void write(ImageOutputStream out) throws IOException, NaturalSetException{
		ImageOutput imageOutput = new ImageOutput(env());
	    imageOutput.paint(image, out);
	}

	
	protected double coverage(int i) throws NaturalSetException{
		return coordinate(i+1) - coordinate(i);
	}
	
	protected double overlappingCoverage(int i, double min, double max) throws NaturalSetException{
		return	((min < coordinate(i+1)) && (max > coordinate(i))) ? (Math.min(coordinate(i+1), max) - Math.max(coordinate(i), min)) : 0.0;
	}
	
	protected double overlappingCoverageFraction(int i, double min, double max) throws NaturalSetException{
		return overlappingCoverage(i, min, max) / coverage(i);
	}
	
	public double[] tiledAverage(double tile, double[] parameter) throws NaturalSetException {		
		double[] average = new double[(int)Math.ceil(basePairDomain().closureCardinality() / (double)tile)];
		int left = snpDomain().min();
		int count = 0;
		double mark = basePairDomain().min();

		for(Integer i : snpDomain()){
			while(coordinate(i) > mark + tile){
				for(int k=left; k<i; k++){
					average[count] +=  parameter[snpDomain().toRelativeCoordinate(k)] / (i-left);
				}
				count++;
				mark += tile;
				left = i-1;
			}
		}
		
		// Complete the last, incomplete, tile
		if(mark < basePairDomain().max()){
			for(int k=left; k<snpDomain().max(); k++){
				average[count] +=  parameter[snpDomain().toRelativeCoordinate(k)] / (snpDomain().max()-left);
			}
		}

		return  average;
	}
		
	public double[] tiledDensity(double tile, double[] parameter) throws NaturalSetException {		
		double[] density = new double[(int)Math.ceil(basePairDomain().closureCardinality() / (double)tile)];
		int left = snpDomain().min();
		int count = 0;
		double mark = basePairDomain().min();

		for(Integer i : snpDomain()){
			while(coordinate(i) > mark + tile){
				for(int k=left; k<i; k++){
					
					density[count] += overlappingCoverageFraction(k, mark, mark + tile) * parameter[snpDomain().toRelativeCoordinate(k)];
				}
				count++;
				mark += tile;
				left = i-1;
			}
		}
		
		// Complete the last, incomplete, tile
		if(mark < basePairDomain().max()){
			for(int k=left; k<snpDomain().max(); k++){
				density[count] += overlappingCoverageFraction(k, mark, basePairDomain().max()) * parameter[snpDomain().toRelativeCoordinate(k)];
			}
		}

		return  density;
	}

	public double[] tiledSnpDensity(double tile) throws NaturalSetException {		
		double[] density = new double[(int)Math.ceil(basePairDomain().closureCardinality() / (double)tile)];
		int left = snpDomain().min();
		int count = 0;
		double mark = basePairDomain().min();

		for(Integer i : snpDomain()){
			while(coordinate(i) > mark + tile){
				for(int k=left; k<i; k++){
					density[count] += overlappingCoverageFraction(k, mark, mark + tile);
				}
				count++;
				mark += tile;
				left = i-1;
			}
		}
		
		// Complete the last, incomplete, tile
		if(mark < basePairDomain().max()){
			for(int k=left; k<snpDomain().max(); k++){
				density[count] += overlappingCoverageFraction(k, mark, basePairDomain().max());
			}
		}

		return  density;
	}
	
	
	protected Frame frameAt(int i) throws NaturalSetException{
		return frames.get(snpDomain().toRelativeCoordinate(i));
	}
		
	protected class Tile {
		protected int position;
		
		protected Frame fleft;
		protected Frame fright;
		
		protected double bpleft;
		protected double bpright;
				
		protected Tile(int position, Frame fleft, Frame fright, double bpleft, double bpright){
			this.position = position;
			this.fleft = fleft;
			this.fright = fright;
			this.bpleft = bpleft;
			this.bpright = bpright;
		}
		
		protected double width(){
			return bpright - bpleft;
		}

		protected double correlation(Tile other) throws NaturalSetException{
			double c = 0.0;
			for(int i=this.fleft.index; i<=this.fright.index; i++){
				for(int j=other.fleft.index; j<=other.fright.index; j++){
					c += frameAt(i).covered(this.bpleft, this.bpright) * frameAt(j).covered(other.bpleft, other.bpright) * pairwiseSnpCorrelation(i, j); 
				}
			}
			c /= (this.width() * other.width());
			return c;
		}
		
		public String toString(){
			return position + " { BP:[" + bpleft + ", " + bpright + "][" + width() +"] SNP[" + fleft.index + ", "+ fright.index + "] }"; 
		}		
	
		public int minSnp(){
			return fleft.index;
		}
		
		public int maxSnp(){
			return fright.index;
		}
		
		public int meanSnp(){
			return (int)Math.round((fright.index + fleft.index) / 2.0);
		}
						
		public int minBp(){
			return (int)Math.round(bpleft);
		}

		public int maxBp(){
			return (int)Math.round(bpright);
		}
		
		public int meanBp(){
			return (int)Math.round((bpright + bpleft) / 2.0);
		}		
	}
	
	protected class Frame {
		protected int index;
		protected double bpleft;
		protected double bpright;

		protected Frame(int index, int prev, int position, int next){
			this.index = index;
			bpleft = (prev + position) / 2.0;
			bpright = (next + position) / 2.0;
		}
		
		protected Frame(int index, int edge, int internal){
			this.index = index;
			if(internal > edge) { //first
				bpleft = edge;
				bpright = (internal + edge) / 2.0;

			} else { // last
				bpleft = (internal + edge) / 2.0;
				bpright = edge;
			}
		}
		
		protected double covered(double min, double max){
			return ((min < bpright) && (max > bpleft)) ? (Math.min(bpright, max) - Math.max(bpleft, min)) : 0.0;
		}
				
		protected double width(){
			return bpright - bpleft;
		}

		public String toString(){
			return index + " : [" + bpleft + ", " + bpright + "] [" + width() +"]"; 
		}
	
	}
}
