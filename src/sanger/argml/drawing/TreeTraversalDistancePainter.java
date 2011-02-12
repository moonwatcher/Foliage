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
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.Date;

import sanger.argml.environment.Environment;
import sanger.argml.graph.model.Statistics;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSetException;

public class TreeTraversalDistancePainter extends PlotPainter{
			
	private int hotspotFactor;
	private Color hotspotColor;

	private Statistics s;
		
	public TreeTraversalDistancePainter(Environment env, Statistics s) throws NaturalSetException {
		super(env);
		this.s = s;
		initialize();		
	}
	
	public TreeTraversalDistancePainter(Environment env, PlotPainter base, Statistics s) throws NaturalSetException {
		super(env, base);
		this.s = s;
		initialize();
	}
			
	protected NaturalDomain snpDomain() { 
		return s.snpDomain(); 
	}

	protected NaturalDomain basePairDomain() { 
		return s.basePairDomain(); 
	}
	
	protected int coordinate(int position) throws NaturalSetException { 
		return s.coordinate(position); 
	}
	
	protected double pairwiseSnpCorrelation(int i, int j) throws NaturalSetException {
		return s.localTreeCorrelation(i, j);
	}
	
	public void paintRecombinationDensity(int x, int y, boolean[] significant) throws NaturalSetException{
		double[] r = tiledDensity(tile, s.recombinationRates());
		paintDistribution(x, y, r, true, significant);
		paintRecombinationHotspot(x, y, r);		
	}
	
	private void paintRecombinationHotspot(int x, int y, double[] distribution) throws NaturalSetException{
		double vmax = 0.0;
		double vmean = 0.0;
		
		for(double v : distribution){
			vmax = Math.max(vmax, v);
			vmean += v;			
		}	vmean /= distribution.length;
		
		Graphics2D lines = (Graphics2D)stage.create();
		lines.translate(x, y);
		lines.setColor(hotspotColor);

		Graphics2D values = (Graphics2D)stage.create();
		values.translate(x + tiles.size() + 4 * padding, y + stage.getFont().getSize() / 2);
		values.setColor(hotspotColor);
		
		int hotLine = (int)Math.round(vmean/vmax*distributionRange) * hotspotFactor;
		lines.drawLine(0, hotLine, tiles.size() + 3 * padding, hotLine);
		
		values.drawString(doubleToString.format(hotspotFactor * vmean) + " hot", 0, hotLine);		
	}	

	public void paintRecombinationDensity(int order, boolean[] significant) throws NaturalSetException{
		paintRecombinationDensity(xPx, yPx + tiles.size() + (order * (distributionRange + padding)) + padding, significant);
	}

	public void paintRecombinationDensity(int order) throws NaturalSetException{
		paintRecombinationDensity(order, null);
	}

		
	public void paintSnpDensity(int x, int y) throws NaturalSetException{
		double[] r = tiledSnpDensity(tile);

		paintDistribution(x, y, r, true);
	}
	
	public void paintSnpDensity(int order) throws NaturalSetException{
		paintSnpDensity(xPx, yPx + tiles.size() + (order * (distributionRange + padding)) + padding);
	}

	public void paintLegend(){
		paintLegend(xPx + scaleRange + 9 * padding + stage.getFont().getSize() * 4 + tiles.size() + axisRange * 2, yPx);
	}
	
	public void paintLegend(int x, int y){
		AffineTransform origin = AffineTransform.getTranslateInstance(x, y);
		
		// frame
		Graphics2D paint = (Graphics2D)stage.create();
		
		paint.setTransform(origin);

		paint.setTransform(origin);
		paint.setColor(lineColor);
		int lineSpace = (int)Math.ceil(stage.getFont().getSize() * 1.1);
		paint.translate(padding, padding + lineSpace);
		paint.drawString("Top:  " + env().stringProperty("Name"), 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Bottom:  " + env().stringProperty("Other"), 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Tile (pixel):  " + env().integerProperty("Tile") +" bp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Min bp:  " + s.basePairDomain().min() +" bp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Max bp:  " + s.basePairDomain().max() +" bp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Min SNP:  " + s.snpDomain().min(), 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Max SNP:  " + s.snpDomain().max(), 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("SNP density:  " + s.basePairDomain().closureCardinality() / s.snpDomain().closureCardinality() +" bp/snp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Grid spacing :  " + env().integerProperty("GridSpacing") + " bp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Time :  " + new Date(System.currentTimeMillis()), 0, 0);
	}
	
	protected void readProperties(){
		super.readProperties();
		this.hotspotFactor = env().integerProperty("HotspotFactor");
		
		// Colors
		this.hotspotColor = new Color(Long.decode(env().stringProperty("HotspotColor")).intValue(), true);
	}

}
