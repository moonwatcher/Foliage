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

import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.Date;

import sanger.argml.environment.Environment;
import sanger.argml.graph.model.HaplotypeSet;
import sanger.math.set.NaturalDomain;
import sanger.math.set.NaturalSetException;

public class RSquarePlotPainter extends PlotPainter{
	private HaplotypeSet hs;
		
	public RSquarePlotPainter(Environment env, HaplotypeSet hs) throws NaturalSetException {
		super(env);
		this.hs = hs;
		initialize();
	}

	public RSquarePlotPainter(Environment env, PlotPainter base, HaplotypeSet hs) throws NaturalSetException {
		super(env, base);
		this.hs = hs;
		initialize();
	}

	protected NaturalDomain snpDomain(){ 
		return hs.snpDomain(); 
	}

	protected NaturalDomain basePairDomain(){
		return hs.basePairDomain();
	}
	
	protected int coordinate(int i) throws NaturalSetException {
		return hs.coordinate(i);
	}
	
	protected double pairwiseSnpCorrelation(int i, int j) throws NaturalSetException{
		return hs.rsquare(i, j);
	}
		
	public void drawLegend(){drawLegend(xPx + scaleRange + 9 * padding + stage.getFont().getSize() * 4 + tiles.size() + axisRange * 2, yPx);};
	
	public void drawLegend(int x, int y){
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
		paint.drawString("Min bp:  " + hs.basePairDomain().min() +" bp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Max bp:  " + hs.basePairDomain().max() +" bp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Min SNP:  " + hs.snpDomain().min(), 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Max SNP:  " + hs.snpDomain().max(), 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("SNP density:  " + hs.basePairDomain().closureCardinality() / hs.snpDomain().closureCardinality() +" bp/snp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Grid spacing :  " + env().integerProperty("GridSpacing") + " bp", 0, 0);
		paint.translate(0, lineSpace);
		paint.drawString("Time :  " + new Date(System.currentTimeMillis()), 0, 0);
	}
	
	public void paintMinorAlleleDensity(int x, int y) throws NaturalSetException{
		double[] r = tiledDensity(tile, hs.minorAlleleFrequencies());
		double[] s = tiledSnpDensity(tile);
		for(int i=0; i<r.length; i++){
			r[i] /= s[i];
		}

		paintDistribution(x, y, r, true);
	}
	
	public void paintMinorAlleleDensity(int order) throws NaturalSetException{
		paintMinorAlleleDensity(xPx, yPx + tiles.size() + (order * (distributionRange + padding)) + padding);
	}

	public void paintSnpDensity(int x, int y) throws NaturalSetException{
		double[] r = tiledSnpDensity(tile);

		paintDistribution(x, y, r, true);
	}
	
	public void paintSnpDensity(int order) throws NaturalSetException{
		paintSnpDensity(xPx, yPx + tiles.size() + (order * (distributionRange + padding)) + padding);
	}

}
