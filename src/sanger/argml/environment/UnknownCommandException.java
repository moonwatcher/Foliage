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

public class UnknownCommandException extends Exception {

	public UnknownCommandException(String message, Instruction command) {
		super(message + "\nsee syntax description.\n" + command + "\n\n" + command.help());
	}

	public UnknownCommandException() {
		// TODO Auto-generated constructor stub
	}

	public UnknownCommandException(String message) {
		super(message);
		// TODO Auto-generated constructor stub
	}

	public UnknownCommandException(Throwable cause) {
		super(cause);
		// TODO Auto-generated constructor stub
	}

	public UnknownCommandException(String message, Throwable cause) {
		super(message, cause);
		// TODO Auto-generated constructor stub
	}

}
