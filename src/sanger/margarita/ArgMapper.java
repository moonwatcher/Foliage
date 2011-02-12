/*
 * ArgMapper.java
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

import java.util.BitSet;
import java.util.Collections;
import java.util.LinkedList;

import JSci.maths.statistics.ChiSqrDistribution;

public class ArgMapper {
    
    // Global variables.
    private int numargs, numcases, numcontrols, numsequences, numbipartitions, nummarkers;
    private double casefreq, controlfreq;
    private double[] markerlocations;
    private ArgStructure[][] args;
    private byte[][] inputsequences;
    final private ChiSqrDistribution chidist = new ChiSqrDistribution(1);
    
    /**
     * Creates a new instance of ArgMapper
     *
     * @param  args  Inferred ARGs.
     * @param  ip    The InputParser containing the sequence and
     * other data from which the ARGs are to inferred.
     */
    public ArgMapper(final LinkedList<ArgStructure>[] args, final InputParser ip) {
        // Store the ARGs.
        numargs = args.length;
        this.args = new ArgStructure[numargs][];
        for (int whicharg = numargs; --whicharg>=0;)
            this.args[whicharg] = args[whicharg].toArray(new ArgStructure[args[whicharg].size()]);
        
        // Extract the useful information from the InputParser.
        numsequences = ip.getNumberOfSequences();
        numbipartitions = numsequences-3;
        numcases = ip.getNumberOfCases();
        numcontrols = ip.getNumberOfControls();
        casefreq = (double)numcases/(double)numsequences;
        controlfreq = (double)numcontrols/(double)numsequences;
        nummarkers = ip.getNumberOfMarkers();
        markerlocations = ip.getMarkerLocations();
        inputsequences = ip.getInputSequences(); // These are used for the Chi-Square test.
    }
 
    
    /**
     * Maps disease loci. Only performs permutations until num cutoff are found
     * that have mapping score exceeding that for the true case control
     * configuration. Set cutoff equal to numpermutations to peform all permutations.
     *
     * @param numpermutations   The maximum number of permtuations to perform.
     * @param cutoff            How many acceptances before stopping.
     */
    public final void mapWithSmartPermutations(final int numpermutations, final int cutoff){
        
        // Create permutations for individuals.
        final LinkedList<Boolean> permlist = new LinkedList<Boolean>();
        final BitSet truecasecontrols = new BitSet(numsequences);
        for (int leaf = numcases; --leaf>=0;) truecasecontrols.set(leaf);
        for (int leaf = numcases/2; --leaf>=0;) permlist.add(true);
        for (int leaf = numcontrols/2; --leaf>=0;) permlist.add(false);
        final BitSet[] permutations = new BitSet[numpermutations];
        for (int perm = numpermutations, leaf = 0; --perm>=0; leaf = 0){
            Collections.shuffle(permlist);
            permutations[perm] = new BitSet(numsequences);
            for (boolean cc : permlist){
                if (cc) permutations[perm].set(leaf++); else leaf++;
                if (cc) permutations[perm].set(leaf++); else leaf++;
            }
        }
        
        // Do the mapping.
        int[][][] marginaltrees;
        double[] mappingscore = new double[nummarkers];
        double[] argpvalue = new double[nummarkers];
        double permscore;
        int permcount;
        Object[] interpretation;
        BitSet chromosomesunder;
        
        System.out.println("%INTERPRETATION");
        System.out.println("MARKER POSITION ARG CUT_CHISQ_SCORE CUT_CHISQ_PVALUE CHROMOSOMES_UNDER_CUT FREQ_CASES FREQ_CONTROLS");
        for (int marker = 0; marker<nummarkers; marker++){
            mappingscore[marker] = 0.0;
            marginaltrees = getMarginalTreesForMarker(marker);
            for (int whicharg = 0; whicharg<numargs; whicharg++){
                interpretation = treeMapWithInterpretation(marginaltrees[whicharg],truecasecontrols);
                chromosomesunder = (BitSet)interpretation[1];
                System.out.print(marker + " " + markerlocations[marker] + " " + whicharg + " " + (Double)interpretation[0] + " " + (1.0-chidist.cumulative((Double)interpretation[0])) + " ");
                // This prints out the chromosomes under the best cut.
                for (int seq = 0; seq<numsequences; seq++){
                    if (chromosomesunder.get(seq)) System.out.print("1");
                    else System.out.print("0");
                }
                System.out.println(" " + (Double)interpretation[2] + " " + (Double)interpretation[3]);
                mappingscore[marker]+=(Double)interpretation[0];
            }
            permcount = 0;
            for (int perm = 0;;){
                permscore = 0.0;
                for (int whicharg = numargs; --whicharg>=0;){
                    if ((permscore+=treeMap(marginaltrees[whicharg],permutations[perm])) >= mappingscore[marker]){
                        permcount++;
                        break;
                    }
                }
                perm++;
                if (permcount>=cutoff || perm==numpermutations) {
                    mappingscore[marker]/=(double)numargs;
                    argpvalue[marker] = (double)permcount/(double)perm;
                    break;
                }
            }
        }
        
        System.out.println("%MAPPING");
        System.out.println("MARKER POSITION ARG_MAP_SCORE PERM_P-VALUE CHI_P-VALUE");
        // Output the mapping results.
        for (int marker = 0; marker<nummarkers; marker++){
            System.out.println(marker + " " + markerlocations[marker] + " " + mappingscore[marker] + " " + argpvalue[marker] + " " + chiSquare(marker,truecasecontrols));
        }
    }
    
    private final int[][][] getMarginalTreesForMarker(final int marker){
        final int[][][] marginaltrees = new int[numargs][3][numsequences-1];
        for (int whicharg = numargs; --whicharg>=0;){
            int parent = numsequences;
            int node = 0;
            final int[] argtreemap = new int[args[whicharg][args[whicharg].length-1].parent1+1];
            for (int leaf = parent; --leaf>=0;) argtreemap[leaf] = leaf;
            for (ArgStructure struct : args[whicharg]){
                switch (struct.t){
                    case Mu : {
                        if (argtreemap[struct.child1]!=-1) argtreemap[struct.parent1] = argtreemap[struct.child1];
                        else argtreemap[struct.parent1] = -1;
                        break;
                    } case Re : {
                        if (argtreemap[struct.child1]!=-1){
                            if (marker<=struct.location) {
                                argtreemap[struct.parent1] = argtreemap[struct.child1];
                                argtreemap[struct.parent2] = -1;
                            } else {
                                argtreemap[struct.parent1] = -1;
                                argtreemap[struct.parent2] = argtreemap[struct.child1];
                            }
                        } else {
                            argtreemap[struct.parent1] = -1;
                            argtreemap[struct.parent2] = -1;
                        }
                        break;
                    } case Co : {
                        if (argtreemap[struct.child1]!=-1){
                            if (argtreemap[struct.child2]!=-1){
                                marginaltrees[whicharg][0][node] = argtreemap[struct.child1];
                                marginaltrees[whicharg][1][node] = argtreemap[struct.child2];
                                marginaltrees[whicharg][2][node++] = parent;
                                argtreemap[struct.parent1] = parent++;
                            } else {
                                argtreemap[struct.parent1] = argtreemap[struct.child1];
                            }
                        } else if (argtreemap[struct.child2]!=-1){
                            argtreemap[struct.parent1] = argtreemap[struct.child2];
                        } else {
                            argtreemap[struct.parent1] = -1;
                        }
                    }
                }
            }
        }
        return marginaltrees;
    }
    
    private final double treeMap(final int[][] marginaltree, final BitSet casecontrols){
        // Return the chi-tree value for the best bipartition.
        final int[] cases = new int[numsequences+numbipartitions]; // This counts the number of cases below each edge.
        final int[] controls = new int[cases.length]; // And this the number of controls.
        // Populate the case/controls bipartition arrays with the number of cases and controls in each bipartition.
        for (int leaf = numsequences; --leaf>=0;){
            if (casecontrols.get(leaf)) cases[leaf]=1;
            else controls[leaf]=1;
        }
        // Now calculate the score for each of the bipartitions, remembering the best score.
        double score;
        double bestscore = 0.0;
        double e0, e1, e2, e3; // Expected values.
        double t1, t2, t3, t4; // (O-E) values.
        int mutcount; // The number of mutants caused by a mutation bipartitioning the tree here.
        int notmutcount;
        for (int j = 0; j<numbipartitions; j++){
            final int parent = marginaltree[2][j];
            cases[parent] = cases[marginaltree[0][j]] + cases[marginaltree[1][j]];
            controls[parent] = controls[marginaltree[0][j]] + controls[marginaltree[1][j]];
            mutcount = cases[parent] + controls[parent];
            notmutcount = numsequences-mutcount;
            e0 = notmutcount*controlfreq;   //marker=0, case=0
            e1 = notmutcount*casefreq;      //marker=0, case=1
            e2 = mutcount*controlfreq;      //marker=1, case=0
            e3 = mutcount*casefreq;         //marker=1, case=1
            t1 = numcontrols-controls[parent]-e0;
            t2 = numcases-cases[parent]-e1;
            t3 = controls[parent]-e2;
            t4 = cases[parent]-e3;
            if ((score = t1*t1/e0 +
                    t2*t2/e1 +
                    t3*t3/e2 +
                    t4*t4/e3) > bestscore) {
                bestscore = score;
            }
        }
        return bestscore;
    }
    
    private final Object[] treeMapWithInterpretation(final int[][] marginaltree, final BitSet casecontrols){
        // Return the chi-tree value for the best bipartition.
        final int[] cases = new int[numsequences+numbipartitions]; // This counts the number of cases below each edge.
        final int[] controls = new int[cases.length]; // And this the number of controls.
        final BitSet[] undercut = new BitSet[cases.length];
        // Populate the case/controls bipartition arrays with the number of cases and controls in each bipartition.
        for (int leaf = numsequences; --leaf>=0;){
            if (casecontrols.get(leaf)) cases[leaf]=1;
            else controls[leaf]=1;
            undercut[leaf] = new BitSet(numsequences);
            undercut[leaf].set(leaf);
        }
        // Now calculate the score for each of the bipartitions, remembering the best score.
        double score;
        double bestscore = 0.0;
        int bestparent = 0;
        double e0, e1, e2, e3; // Expected values.
        double t1, t2, t3, t4; // (O-E) values.
        int mutcount; // The number of mutants caused by a mutation bipartitioning the tree here.
        int notmutcount;
        for (int j = 0; j<numbipartitions; j++){
            final int parent = marginaltree[2][j];
            cases[parent] = cases[marginaltree[0][j]] + cases[marginaltree[1][j]];
            controls[parent] = controls[marginaltree[0][j]] + controls[marginaltree[1][j]];
            undercut[parent] = (BitSet)(undercut[marginaltree[0][j]].clone());
            undercut[parent].or(undercut[marginaltree[1][j]]);
            mutcount = cases[parent] + controls[parent];
            notmutcount = numsequences-mutcount;
            e0 = notmutcount*controlfreq;   //marker=0, case=0
            e1 = notmutcount*casefreq;      //marker=0, case=1
            e2 = mutcount*controlfreq;      //marker=1, case=0
            e3 = mutcount*casefreq;         //marker=1, case=1
            t1 = numcontrols-controls[parent]-e0;
            t2 = numcases-cases[parent]-e1;
            t3 = controls[parent]-e2;
            t4 = cases[parent]-e3;
            if ((score = t1*t1/e0 +
                    t2*t2/e1 +
                    t3*t3/e2 +
                    t4*t4/e3) > bestscore) {
                bestscore = score;
                bestparent = parent;
            }
        }
        return new Object[]{bestscore,undercut[bestparent],(double)cases[bestparent]/(double)numcases,(double)controls[bestparent]/(double)numcontrols};
    }
    
    private final double chiSquare(final int marker, final BitSet casecontrols){
        int o1 = 0;
        int o2 = 0;
        int o3 = 0;
        int o4 = 0;
        int gotcases = numcases;
        int gotcontrols = numcontrols;
        for (int seq = numsequences; --seq>=0;){ // Loop over the case sequences.
            if (casecontrols.get(seq)){
                switch (inputsequences[seq][marker]){
                    case 0 : {o1++; break;}
                    case 1 : {o2++; break;}
                    case 2 : {if (seq%2==0) o1++; else o2++; break;}
                    case 3 : gotcases--;
                }
            } else {
                switch (inputsequences[seq][marker]){
                    case 0 : {o3++; break;}
                    case 1 : {o4++; break;}
                    case 2 : {if (seq%2==0) o3++; else o4++; break;}
                    case 3 : gotcontrols--;
                }
            }
        }
        // Calculate the expected values.
        double t1 = (double)(o1+o3);
        double t2 = (double)(o2+o4);
        final double e1 = t1*((double)gotcases)/((double)(gotcases+gotcontrols));
        final double e2 = t2*((double)gotcases)/((double)(gotcases+gotcontrols));
        final double e3 = t1*((double)gotcontrols)/((double)(gotcases+gotcontrols));
        final double e4 = t2*((double)gotcontrols)/((double)(gotcases+gotcontrols));
        // Now calculate chisquare.
        t1 = ((double)o1-e1);
        t2 = ((double)o2-e2);
        final double t3 = ((double)o3-e3);
        final double t4 = ((double)o4-e4);
        return 1.0-chidist.cumulative(t1*t1/e1+t2*t2/e2+t3*t3/e3+t4*t4/e4);
    }
}

