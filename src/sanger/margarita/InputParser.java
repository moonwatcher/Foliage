/*
 * InputParser.java
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

import java.io.IOException;
import java.util.LinkedList;

import sanger.argml.environment.Environment;
import sanger.argml.io.TextInput;

public class InputParser {
    
    // Global variables.
    private int numcases, numcontrols, numsequences, nummarkers, numunphasedcharacters; // Properties of the data.
    private double[] markerlocations; // Locations of the markers.
    private byte[][] inputsequences; // The input sequences with alleles as 0 and 1, unphased as 2 and missing data 3.
    private short[][] allelecounts; // Number of 0, 1 and missing alleles.
    private LinkedList<int[]> possiblemutations; // Mutations that can be made immediately.
    
    // ************************************************ INITIALISATION METHODS ************************************************
    
    /** Creates a new instance of InputParser */
    public InputParser() {}
    
    /**
     * Reads in a file that is in Margarita format.
     * The data structures are initialised that will be acted upon
     * by ArgBuilder.
     *
     * @param  filename The file that contains the marker
     *  locations and sequences in Margarita format.
     */
    public final void parseFile(final Environment env, final TextInput input){
        try {            
            // Read in and parse the first line of the input.
            final String[] s = input.reader().readLine().split(" ");
            // Check the first line is OK.
            if (s.length!=3) env.log().printError("First line of file is not correct. Have you given the number of cases, controls, and markers?");
            // Initialise the data structures that we can populate now and when reading in the marker locations.
            numcases = Integer.parseInt(s[0]);
            numcontrols = Integer.parseInt(s[1]);
            numsequences = numcases+numcontrols;
            nummarkers = Integer.parseInt(s[2]);
            markerlocations = new double[nummarkers];
            // Read in the marker locations.
            try {
                for (int marker = 0; marker<nummarkers; marker++) markerlocations[marker] = Double.parseDouble(input.reader().readLine());
            } catch (Exception e) {
            	env.log().printError("Your positions of markers are not correct. Are there enough of them?\n" + e);
            }
            for (int marker = nummarkers-1; marker>=1;){
                if (markerlocations[marker]<=markerlocations[--marker]){
                	env.log().printError("Positions of markers are not in increasing order. Reformat input file.");
                    break;
                }
            }
            // Initialise the data structures that we populate when reading in the sequences.
            inputsequences = new byte[numsequences][nummarkers];
            allelecounts = new short[3][nummarkers];
            // Read in the sequences.
            try {
                final char[] buffer = new char[nummarkers]; // A buffer for the sequences to be read into.
                for (int seq = 0; seq<numsequences; seq++){
                	input.reader().read(buffer,0,nummarkers); // Read in the next line, this is faster than in.readLine();
                	input.reader().read(); // And skip the newline character.
                    for (int marker = nummarkers; --marker>=0;){
                        switch(buffer[marker]){
                            case '0' : {
                                allelecounts[0][marker]++;
                                break;
                            } case '1' : {
                                allelecounts[1][marker]++;
                                inputsequences[seq][marker] = 1;
                                break;
                            } case '2' : {
                            } case 'U' : {
                                if (seq%2==1 && inputsequences[seq-1][marker]<2)
                                	env.log().printError("You are not writing unphased sequences in the correct way. Reformat input file.");
                                allelecounts[2][marker]++;
                                inputsequences[seq][marker] = 2; // 2 means unphased.
                                numunphasedcharacters++;
                                break;
                            } case '3' : {
                            } case 'M' : {
                                allelecounts[2][marker]++;
                                inputsequences[seq][marker] = 3; // 3 means missing.
                                numunphasedcharacters++;
                                break;
                            } default : {
                            	env.log().printError("Unrecognised character " + buffer[marker] + " in sequence " + seq + ", marker " + marker + ". Correct this.");
                            }
                        }
                    }
                }
            } catch (Exception e) {
            	env.log().printError("Your sequences are not correct. Are there enough of them? Do they have the correct number of markers?\n" + e);
            }
            // Now check whether all positions are phase resolvable. If there is a position in a sequence that is unphased, when all the others
            // are missing or are unphased, then this will cause a crash. Need to randomly set one unphased genotype to 0/1.
            for (int mar = 0; mar<nummarkers; mar++){
                if (allelecounts[0][mar]==0 && allelecounts[1][mar]==0 && allelecounts[2][mar]>0){
                    boolean ok = false;
                    for (int seq = 0; seq<numsequences && !ok; seq+=2){
                        if (inputsequences[seq][mar]==2){
                            inputsequences[seq][mar] = 0;
                            inputsequences[seq+1][mar] = 1;
                            allelecounts[2][mar]-=2;
                            allelecounts[0][mar]++;
                            allelecounts[1][mar]++;
                            ok = true; // Then there is a sequence with an unphased genotype.
                        }
                    }
                    if (!ok){
                        // Then all the data is missing at this position.
                        if (inputsequences[0][mar]!=3 && inputsequences[1][mar]!=3) env.log().printError("Bad column of missing data in input file.");
                        inputsequences[0][mar] = 0;
                        inputsequences[1][mar] = 0;
                        allelecounts[2][mar]-=2;
                        allelecounts[0][mar]+=2;
                    }
                }
            }
            // Initialise the data structures that we populate once the file has been read.
            possiblemutations = new LinkedList<int[]>();
            for (int marker = nummarkers; --marker>=0;){
                // Check whether we can mutate this marker.
                if (allelecounts[2][marker]==0){
                    if (allelecounts[1][marker]==1) { // Then we can mutate this.
                        for (int seq = numsequences; --seq>=0;){ // Find the sequence with the singleton.
                            if (inputsequences[seq][marker]==1) {
                                possiblemutations.add(new int[]{seq,marker,0});
                                break;
                            }
                        }
                    } else if (allelecounts[0][marker]==1) { // Then we can mutate this.
                        for (int seq = numsequences; --seq>=0;){ // Find the sequence with the singleton.
                            if (inputsequences[seq][marker]==0) {
                                possiblemutations.add(new int[]{seq,marker,1});
                                break;
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {env.log().printError("Error while reading file.\n" + e);}
    }
    
    // ************************************************ INTERFACE ************************************************
    
    /**
     * How many case haplotype sequences are there in the data?
     *
     * @return  the number of case haplotype sequences.
     */
    public final int getNumberOfCases(){
        return numcases;
    }
    
    /**
     * How many control haplotype sequences are there in the data?
     *
     * @return  the number of control haplotype sequences.
     */
    public final int getNumberOfControls(){
        return numcontrols;
    }    
    
    /**
     * How many haplotype sequences are there in the data?
     *
     * @return  the number of haplotype sequences.
     */
    public final int getNumberOfSequences(){
        return numsequences;
    }
    
    /**
     * How many markers are typed in the data?
     *
     * @return  the number of markers.
     */
    public final int getNumberOfMarkers(){
        return nummarkers;
    }
    
    /**
     * How many unphased or missing characters are there in the data?
     *
     * @return  the number of unphased or missing characters.
     */
    public final int getNumberOfUnphasedCharacters(){
        return numunphasedcharacters;
    }
    
    /**
     * What are the positions of the markers?
     *
     * @return  the positions of the markers.
     */
    public final double[] getMarkerLocations(){
        return markerlocations;
    }
    
    /**
     * What are the input sequences?
     *
     * @return  the input sequences, where the first dimension is the sequence index,
     * and the second dimension is the marker index.
     * 0 and 1 are the alleles, 2 means unphased and 3 means missing data.
     */
    public final byte[][] getInputSequences(){
        return inputsequences;
    }
    
    /**
     * How many of each allele at each marker?
     *
     * @return  the number of alleles at each markers. [0][marker] is
     * the number of 0 alleles at marker, similarly for [1][marker]
     * and [2][marker] is the number of missing data at that marker.
     */
    public final short[][] cloneAlleleCounts(){
        final short[][] allelecountsclone = new short[3][nummarkers];
        for (int allele = 3; --allele>=0;) System.arraycopy(allelecounts[allele],0,allelecountsclone[allele],0,nummarkers);
        return allelecountsclone;
    }
    
    /**
     * Which mutations are immediately possible?
     *
     * @return  the number of alleles at each markers. [0][marker] is
     * the number of 0 alleles at marker, similarly for [1][marker]
     * and [2][marker] is the number of missing data at that marker.
     */
    public final LinkedList<int[]> clonePossibleMutations(){
        final LinkedList<int[]> possiblemutationsclone = new LinkedList<int[]>();
        for (int[] mut : possiblemutations) possiblemutationsclone.add(new int[]{mut[0],mut[1],mut[2]});
        return possiblemutationsclone;
    }
}
