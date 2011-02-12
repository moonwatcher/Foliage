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

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Iterator;

import javax.imageio.ImageIO;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;

import sanger.argml.environment.Environment;
import sanger.argml.environment.Environmental;


public class ImageOutput extends Environmental{
	protected ImageWriter writer;

	
	public ImageOutput(Environment env){
		super(env);
		
		Iterator<ImageWriter> writers = ImageIO.getImageWritersByFormatName("png");
		this.writer = (ImageWriter)writers.next();	
	}

	public ImageOutput(Environment env, String format){
		super(env);
		
		Iterator<ImageWriter> writers = ImageIO.getImageWritersByFormatName(format);
		this.writer = (ImageWriter)writers.next();
	}
	
	public void paint(BufferedImage image, ImageOutputStream stream) throws IOException {
	    writer.setOutput(stream);	
	    writer.write(image);
	    stream.flush();
	}
}
