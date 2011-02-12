/*
 * ArgStructure.java
 *
 * Copyright (c) 2007 Genome Research Ltd.
 * Author: Mark Minichiello
 *
 * THIS SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
 * OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * This code is free software; you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation.
 *
 * Any redistribution or derivation in whole or in part including any substantial portion
 * of this code must include this copyright and permission notice.
 *
 */

package sanger.margarita;

public class ArgStructure {
    
    public enum Type {Mu, Co, Re;}
    public int time, child1, child2, parent1, parent2, location;
    public Type t;
    
    /** Creates a new instance of ArgStructure */
    ArgStructure(final int time, final Type t, final int child1, final int child2, final int parent1, final int parent2, final int location){
        this.time = time; this.t = t; this.child1 = child1; this.child2 = child2; this.parent1 = parent1; this.parent2 = parent2; this.location = location;
    }
    
    /**
     * Prints the ArgStructure to the screen.
     */
    public final String toString(){
        switch(this.t){
            case Mu : return time + " mu " + child1 + " " + parent1 + " " + location;
            case Co : return time + " co " + child1 + " " + child2 + " " + parent1;
            case Re : return time + " re " + child1 + " " + parent1 + " " + parent2 + " " + location;
        }
        return null;
    }
}
