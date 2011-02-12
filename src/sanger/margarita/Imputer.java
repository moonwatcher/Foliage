/*
 * Imputer.java
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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.LinkedList;

public class Imputer {
    
    int numargs, numsequences, numcases, numcontrols, nummarkers;
    double[] markerlocations;
    byte[][] inputsequences;
    byte[][][] leafsequences;
    LinkedList<ArgStructure>[] args;
    
    /** Creates a new instance of Imputer */
    public Imputer(final LinkedList<ArgStructure>[] args, final InputParser ip) {
        this.args = args;
        numargs = args.length;
        numsequences = ip.getNumberOfSequences();
        numcases = ip.getNumberOfCases();
        numcontrols = ip.getNumberOfControls();
        nummarkers = ip.getNumberOfMarkers();
        markerlocations = ip.getMarkerLocations();
        inputsequences = ip.getInputSequences();
        leafsequences = calculateImputations();
    }
    
    private final byte[][][] calculateImputations(){    
        Hashtable<Integer,byte[]> currentsequences = new Hashtable<Integer,byte[]>();
        final byte[][][] leafsequences = new byte[numargs][][];
        for (int whicharg = 0; whicharg<numargs; whicharg++){
            // Set the currentsequences to be the leaf sequences.
            currentsequences.clear();
            for (int seq = numsequences; --seq>=0;){
                final byte[] sequence = new byte[nummarkers];
                System.arraycopy(inputsequences[seq],0,sequence,0,nummarkers);
                currentsequences.put(seq,sequence);
            }
            // Make upward pass of the ARG.
            final ArgStructure[] arg = args[whicharg].toArray(new ArgStructure[args[whicharg].size()]);
            for (ArgStructure struct : arg){
                switch(struct.t){
                    case Mu : {
                        final byte[] sequence = currentsequences.remove(struct.child1);
                        switch (sequence[struct.location]){
                            case 0 : {
                                sequence[struct.location] = 1;
                                break;
                            } case 1 : {
                                sequence[struct.location] = 0;
                                break;
                            }
                        }
                        currentsequences.put(struct.parent1,sequence);
                        break;
                    } case Co : {
                        final byte[] sequence1 = currentsequences.remove(struct.child1); // Recycle this one.
                        final byte[] sequence2 = currentsequences.remove(struct.child2);
                        for (int marker = nummarkers; --marker>=0;){
                            if (sequence1[marker]>1) sequence1[marker] = sequence2[marker];
                            else if (sequence1[marker]!=sequence2[marker] && sequence2[marker]<=1) System.out.println("ARG is not correct ERROR 1.");
                        }
                        currentsequences.put(struct.parent1,sequence1);
                        break;
                    } case Re : {
                        final byte[] sequence = currentsequences.remove(struct.child1);
                        final byte[] rightsequence = new byte[nummarkers];
                        System.arraycopy(sequence,struct.location+1,rightsequence,struct.location+1,nummarkers-struct.location-1);
                        for (int marker = 0; marker<=struct.location; marker++) rightsequence[marker] = 4; // 4 means undefined.
                        for (int marker = struct.location+1; marker<nummarkers; marker++) sequence[marker] = 4;
                        currentsequences.put(struct.parent1,sequence);
                        currentsequences.put(struct.parent2,rightsequence);
                    }
                }
            }
            // Make downward pass of the arg.
            for (int str = arg.length; --str>=0;){
                switch(arg[str].t){
                    case Mu : {
                        final byte[] sequence = currentsequences.remove(arg[str].parent1);
                        if (sequence[arg[str].location]==0) sequence[arg[str].location] = 1;
                        else sequence[arg[str].location] = 0;
                        currentsequences.put(arg[str].child1,sequence);
                        break;
                    } case Co : {
                        final byte[] sequence1 = currentsequences.remove(arg[str].parent1);
                        final byte[] sequence2 = new byte[nummarkers];
                        System.arraycopy(sequence1,0,sequence2,0,nummarkers);
                        currentsequences.put(arg[str].child1,sequence1);
                        currentsequences.put(arg[str].child2,sequence2);
                        break;
                    } case Re : {
                        final byte[] sequence1 = currentsequences.remove(arg[str].parent1);
                        final byte[] sequence2 = currentsequences.remove(arg[str].parent2);
                        for (int marker = arg[str].location+1; marker<nummarkers; marker++)
                            sequence1[marker] = sequence2[marker];
                        currentsequences.put(arg[str].child1,sequence1);
                    }
                }
            }
            // Currentsequences contains the leaf sequences. Store these leaf sequences.
            leafsequences[whicharg] = new byte[numsequences][nummarkers];
            for (int seq = numsequences; --seq>=0;){
                System.arraycopy(currentsequences.get(seq),0,leafsequences[whicharg][seq],0,nummarkers);
            }
        }
        return leafsequences;
    }
    
    
    /**
     * Returns the leaf sequences as imputed by the ARGs.
     *
     * @return  The imputed leaf sequences.
     */
    public final byte[][][] getImputations(){
        return leafsequences;
    }
    
    /**
     * Outputs the imputed leaf sequences to a file.
     *
     * @param the file to which the imputations are output.
     */
    public final void outputImputations(final String outputfilename){
        final int numargs = args.length;
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(outputfilename));
            for (int whicharg = 0; whicharg<numargs; whicharg++){
                // Output the leafsequences.
                out.write("%IMPUTATION\n");
                out.write(numcases + " " + numcontrols + " " + nummarkers + "\n");
                for (double location : markerlocations) out.write(location + "\n");
                for (int seq = 0; seq<numsequences; seq++){
                    for (int mar = 0; mar<nummarkers; mar++){
                        if (leafsequences[whicharg][seq][mar]==0) out.write("0");
                        else if (leafsequences[whicharg][seq][mar]==1) out.write("1");
                        else if (leafsequences[whicharg][seq][mar]==2) out.write("U"); // Very unusual to get this...
                        else if (leafsequences[whicharg][seq][mar]==3) out.write("M"); // ...or this.
                        else System.err.println("Error.");
                    }
                    out.write("\n");
                }
            }
            out.close();
            // Output the consensus imputation.
            final byte[][][] genotypevote = new byte[3][numsequences/2][nummarkers];
            for (int arg = 0; arg<numargs; arg++){
                for (int seq = 0; seq<numsequences; seq+=2){
                    for (int mar = 0; mar<nummarkers; mar++){
                        if (leafsequences[arg][seq][mar]==2) genotypevote[1][seq/2][mar]++; // Het. Very unusual.
                        else if (leafsequences[arg][seq][mar]+leafsequences[arg][seq+1][mar]==2) genotypevote[2][seq/2][mar]++; // Hom 1.
                        else if (leafsequences[arg][seq][mar]+leafsequences[arg][seq+1][mar]==1) genotypevote[1][seq/2][mar]++; // Het.
                        else if (leafsequences[arg][seq][mar]+leafsequences[arg][seq+1][mar]==0) genotypevote[0][seq/2][mar]++; // Hom 0.
                    }
                }
            }
            final char[][] consensus = new char[numsequences][nummarkers];
            for (int ind = 0; ind<numsequences/2; ind++){
                for (int mar = 0; mar<nummarkers; mar++){
                    if (genotypevote[2][ind][mar]>=genotypevote[1][ind][mar] && genotypevote[2][ind][mar]>=genotypevote[0][ind][mar]){
                        consensus[ind*2][mar] = '1';
                        consensus[ind*2+1][mar] = '1';
                    } else if (genotypevote[1][ind][mar]>=genotypevote[2][ind][mar] && genotypevote[1][ind][mar]>=genotypevote[0][ind][mar]){
                        consensus[ind*2][mar] = 'U';
                        consensus[ind*2+1][mar] = 'U';
                    } else if (genotypevote[0][ind][mar]>=genotypevote[2][ind][mar] && genotypevote[0][ind][mar]>=genotypevote[1][ind][mar]){
                        consensus[ind*2][mar] = '0';
                        consensus[ind*2+1][mar] = '0';
                    } else System.err.println("Error.");
                }
            }
            out = new BufferedWriter(new FileWriter(outputfilename + ".consensusimputation"));
            out.write(numcases + " " + numcontrols + " " + nummarkers + "\n");
            for (double location : markerlocations ) out.write(location + "\n");
            for (int seq = 0; seq<numsequences; seq++){
                for (int mar = 0; mar<nummarkers; mar++)
                    out.write(consensus[seq][mar]);
                out.write("\n");
            }
            out.close();
        } catch (IOException e){System.err.println("IO Error " + e);}
    }
}