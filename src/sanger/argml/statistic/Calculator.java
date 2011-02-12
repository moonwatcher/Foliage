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

package sanger.argml.statistic;

import sanger.argml.environment.Environment;


public class Calculator {
	
	public static double PearsonCorrelationCoefficient(Environment env, double[] x, double[] y) {
		double correlation = -100;
		double EPS=0.001;
		
		try {
			double sum_sq_x = 0;
			double sum_sq_y = 0;
			double sum_coproduct = 0;
			double mean_x = x[0];
			double mean_y = y[0];
			for(int i=1;i<x.length;i++){
			    double sweep = i / (i + 1.0);
			    double delta_x = x[i] - mean_x;
			    double delta_y = y[i] - mean_y;
			    sum_sq_x += delta_x * delta_x * sweep;
			    sum_sq_y += delta_y * delta_y * sweep;
			    sum_coproduct += delta_x * delta_y * sweep;
			    mean_x += delta_x / (i + 1.0);
			    mean_y += delta_y / (i + 1.0);				
			}
			double pop_sd_x = Math.sqrt(sum_sq_x / x.length);
			double pop_sd_y = Math.sqrt(sum_sq_y / y.length);
			if (pop_sd_x > EPS && pop_sd_y > EPS) {
				double cov_x_y = sum_coproduct / x.length;
				correlation = cov_x_y / (pop_sd_x * pop_sd_y);
			}
		} catch (Exception e) {
			correlation = -100;
			env.log().printError(e);
		}
				
		return correlation;
	}

	public static double max(double[] x){
	    double max = x[0];
	    for (int i=1; i<x.length; i++) {
	        if (x[i] > max) max = x[i];
	    }
	    return max;
	}
	
	public static double mean(double[] x) {
	    double sum = 0;
	    for (int i=0; i<x.length; i++) sum += x[i];
	    return sum / x.length;
	}	

	public static double mean(int[] x) {
	    double sum = 0;
	    for (int i=0; i<x.length; i++) sum += x[i];
	    return sum / x.length;
	}
	
	public static double sum(double[] x){
	    double sum = 0;
	    for (int i=0; i<x.length; i++) sum += x[i];
	    return sum;
	}

	public static int sum(int[] x){
	    int sum = 0;
	    for (int i=0; i<x.length; i++) sum += x[i];
	    return sum;
	}

}
