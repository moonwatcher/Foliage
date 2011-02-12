/*
 * ArgBuilderForUnphasedData.java
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

import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.Random;

import sanger.argml.environment.Environment;
import sanger.argml.io.TextOutput;

public class ArgBuilderForUnphasedData {
    // lg8 Add a PrintOutput abstract output 
	private TextOutput output;
	private Environment env;
	
    private int numsequences, numcases, numcontrols, nummarkers, nmm1; // Number of sequences and markers.
    private double[] markerlocations;
    private byte[][] inputsequences; //
    private short[][] allelecounts; // The number of currentsequences with alleles as 0, 1 and missing.
    private double[][] distancematrix; // The distances between all pairs of SNPs.
    private double longestpossible; // The longest possible length of a shared segment that does not span the whole distance.
    private LinkedList<int[]> possiblemutations; // The list of possible mutations.
    private Hashtable<Integer,int[]> currentsequences; // The sequences of the live edges in the ARG.
    private LinkedList<ArgStructure>[] args; // All the ARGs.
    private LinkedList<ArgStructure> currentarg; // The ARG currently under construction.
    private int numcoalescences, numrecombinations, numgeneconversions, nummutations, time, nextparent; // Particular to the current arg.
    private Hashtable<Integer,int[]> startends; // The extent of each sequence.
    private Map.Entry<Integer,int[]>[] startendsarray; // This holds a permuted array of the startends.
    private LinkedList<Integer> coalescenceedges; // The edges that should be tested for coalescence.
    private int[] edgepointers; // For each unphased position, point to the edge on which it occurs.
    private int[] conflicts; // For each unphased position, what does it conflict with?
    private int[] dependencies; // For each unphased position, what is it the same as?
    private static Random rand = new Random(); // Random number generator.
    
    // Parameters of the algorithm.
    private static final double HEURISTICP = 0.9; // How frequently the heuristic is used.
    private static boolean VERBOSE = false; // Whether to print out the ARG as it is constructed.
    
    /** Creates a new instance of ArgBuilderForUnphasedData */
    public ArgBuilderForUnphasedData(Environment env, TextOutput output) {
    	this.output = output;
    	this.env = env;
    }
    
    // ******************************************************** CODE TO DRIVE ARG CONSTRUCTION ********************************************************
    
    /**
     * Constructs ARGs. From either phased or unphased data and missing data.
     * If data is phased and there are no missing characters, this implementation
     * is inefficient in memory and time.
     *
     * @param  numargs  The number of ARGs to construct
     * @param  ip       The InputParser containing the sequence and
     * other data from which the ARGs are to inferred.
     */
    public final void buildArgs(final int numargs, final InputParser ip){
        args = new LinkedList[numargs];
        // Initialise the data structures that come immediately from the input parser.
        numsequences = ip.getNumberOfSequences();
        numcases = ip.getNumberOfCases();
        numcontrols = ip.getNumberOfControls();
        nummarkers = ip.getNumberOfMarkers();
        nmm1 = nummarkers-1;
        
        // Initialise and populate the data structures that come from the marker locations.
        markerlocations = ip.getMarkerLocations();
        distancematrix = new double[nummarkers][nummarkers];
        for (int marker2 = nummarkers; --marker2>=1;)
            for (int marker1 = marker2; --marker1>=0;)
                distancematrix[marker1][marker2] = markerlocations[marker2]-markerlocations[marker1];
        if (distancematrix[1][nmm1]>distancematrix[0][nmm1-1]) longestpossible = distancematrix[1][nmm1];
        else longestpossible = distancematrix[0][nmm1-1];
        
        // The data structures that need to be repopulated for every ARG that we construct.
        startends = new Hashtable<Integer,int[]>(numsequences);
        currentsequences = new Hashtable<Integer,int[]>();
        coalescenceedges = new LinkedList<Integer>();
        edgepointers = new int[ip.getNumberOfUnphasedCharacters()+2]; // These are numbered from 2.
        conflicts = new int[edgepointers.length];
        dependencies = new int[edgepointers.length];
        
        // Get some things to help buildilng the final data structures.
        inputsequences = ip.getInputSequences();
        
        // Populate the final data structures and build the ARGs.
        output.writer().println("%ARGINFERENCE");
        output.writer().println("SEQS SNPS MUTS COAS RECS GECS TRCS SECS HEURP");
        long starttime;
        for (int arg = numargs; --arg>=0;){
            //  Initialise the data structures that need to be reset every time.
            numcoalescences = 0;
            numrecombinations = 0;
            numgeneconversions = 0;
            nummutations = 0;
            time = 0;
            nextparent = numsequences;
            currentarg = new LinkedList<ArgStructure>();
            allelecounts = ip.cloneAlleleCounts();
            possiblemutations = ip.clonePossibleMutations();
            
            // Initialise and populate the data structures that are calculated from the sequences and have to be reset every time.
            currentsequences.clear();
            startends.clear();
            coalescenceedges.clear();
            int whichunphasedcharacter = 2; // These are numbered from 2.
            int[] sistersequence = null; // This is a helper data structure.
            for (int seq = 0; seq<numsequences; seq++){ // Must do loop in this order to get brother and sister sequences correct.
                final int[] currentsequence = new int[nummarkers];
                final boolean hassister = (seq%2==1);
                for (int marker = 0; marker<nummarkers; marker++){ // To get correct ordering.
                    switch (inputsequences[seq][marker]){
                        case 0 : {
                            break;
                        } case 1 : {
                            currentsequence[marker] = 1;
                            break;
                        } case 2 : { // Unphased character.
                            if (hassister) { // We can only populate the conflicts once both sister and brother have been read in.
                                conflicts[whichunphasedcharacter] = sistersequence[marker];
                                conflicts[conflicts[whichunphasedcharacter]] = whichunphasedcharacter;
                            }
                            edgepointers[whichunphasedcharacter] = seq;
                            dependencies[whichunphasedcharacter] = whichunphasedcharacter;
                            currentsequence[marker] = whichunphasedcharacter++;
                            break;
                        } case 3 : { // Unphased character or missing data.
                            edgepointers[whichunphasedcharacter] = seq;
                            dependencies[whichunphasedcharacter] = whichunphasedcharacter;
                            conflicts[whichunphasedcharacter] = 0;
                            currentsequence[marker] = whichunphasedcharacter++;
                        }
                    }
                }
                if (seq%2==0) sistersequence = currentsequence;
                currentsequences.put(seq,currentsequence);
                startends.put(seq,new int[]{0,nmm1});
                coalescenceedges.add(seq);
            }
            sistersequence = null;
            // Construct one ARG.
            starttime = System.nanoTime();
            buildArg();
            output.writer().println(numsequences + " " + nummarkers + " " + nummutations + " " + numcoalescences + " " +
                    numrecombinations + " " + numgeneconversions + " " + (numrecombinations+2*numgeneconversions) + " " +
                    ((System.nanoTime()-starttime)/(double)1000000000) + " " + HEURISTICP);
            args[arg] = currentarg;
            if (VERBOSE) output.writer().println();
        }
    }
    
    private final void buildArg(){ // Builds one ARG.
        // 1. Do all possible mutations.
        // 2. Do a coalescence then goto 1. If no coalsecences possible, goto 3.
        // 3. Do a recombination followed by a coalescence. Goto 1.
        Integer[] possiblecoalescence;
        while (true){
            while (true){
                // Check to see whether any mutations are possible.
                if (!possiblemutations.isEmpty()) doMutations();
                // We will now perform a coalescence or a recombination then coalescences.
                permuteStartEnds(); // This randomises the order in which coalescences and recombinations are performed. Done here to reduce ~half the number of permutations performed.
                // Check to see whether any coalescences are possible.
                possiblecoalescence = getACoalescence();
                if (possiblecoalescence==null) break;
                makeCoalescence(possiblecoalescence);
                //output.writer().println("Numrecs " + numrecombinations + " Numcgs " + numgeneconversions + " Numcos " + numcoalescences);
                if (isFinished()) return;
            }
            possiblecoalescence = doARecombination();
            makeCoalescence(possiblecoalescence);
            //output.writer().println("Numrecs " + numrecombinations + " Numcgs " + numgeneconversions + " Numcos " + numcoalescences);
            if (isFinished()) return;
        }
    }
    
    // ******************************************************** METHOD TO DO MUTATIONS ********************************************************
    
    private final void doMutations(){
        // Mutations are of the form {edge, location, allele (that we are mutating to)}
        // Initialise data structures and shuffle the order in which we will do the mutations.
        Collections.shuffle(possiblemutations);
        final int[][] mutationsarray = possiblemutations.toArray(new int[possiblemutations.size()][]);
        possiblemutations.clear();
        int[] sequence;
        int[] startend;
        Integer edgekey;
        // Perform the mutations.
        for (int mut1 = mutationsarray.length; --mut1>=0;){
            nummutations++;
            // Get the edge that will be mutated.
            edgekey = new Integer(mutationsarray[mut1][0]);
            coalescenceedges.remove(edgekey);
            sequence = currentsequences.remove(edgekey);
            startend = startends.remove(edgekey);
            
            // Perform the mutation.
            sequence[mutationsarray[mut1][1]] = mutationsarray[mut1][2];
            for (int marker = startend[0]; marker<=startend[1]; marker++)
                if (sequence[marker]>1)
                    edgepointers[sequence[marker]] = nextparent;
            // Put the mutated edge into the datastructures.
            coalescenceedges.addFirst(nextparent);
            currentsequences.put(nextparent,sequence);
            startends.put(nextparent,startend);
            currentarg.add(new ArgStructure(time++,ArgStructure.Type.Mu,mutationsarray[mut1][0],-1,nextparent,-1,mutationsarray[mut1][1]));
            // Update the other mutations that are also on that edge.
            for (int mut2 = mut1; --mut2>=0;)
                if (mutationsarray[mut2][0]==mutationsarray[mut1][0])
                    mutationsarray[mut2][0] = nextparent;
            // Update the allele counts.
            if (mutationsarray[mut1][2]==0){
                allelecounts[0][mutationsarray[mut1][1]]++;
                allelecounts[1][mutationsarray[mut1][1]]--;
            } else {
                allelecounts[1][mutationsarray[mut1][1]]++;
                allelecounts[0][mutationsarray[mut1][1]]--;
            }
            nextparent++;
            if (VERBOSE) output.writer().println(currentarg.getLast());
        }
    }
    
    // ******************************************************** METHODS TO DO COALESCENCES ********************************************************
    
    private final Integer[] getACoalescence(){ // Returns two edges that may be coalesced.
        int[] sequence1, sequence2, startend1, startend2, overlap;
        Integer seq1edge, seq2edge;
        final boolean[] startendschecked = new boolean[startendsarray.length];
        // Shuffle the current coalescence edges.
        Collections.shuffle(coalescenceedges);
        // Find a possible coalescence.
        final ListIterator<Integer> li = coalescenceedges.listIterator();
        while (li.hasNext()){
            seq1edge = li.next();
            sequence1 = currentsequences.get(seq1edge);
            startend1 = startends.get(seq1edge);
            for (int seq2 = startendschecked.length; --seq2>=0;) {
                if (startendschecked[seq2]==true) continue;
                seq2edge = startendsarray[seq2].getKey();
                if (seq1edge.equals(seq2edge)){
                    startendschecked[seq2]=true;
                    continue;
                }
                if ((overlap = getOverlap(startend1,startendsarray[seq2].getValue()))!=null
                        && possibleToCoalesce(sequence1,currentsequences.get(seq2edge),overlap)){
                    return new Integer[] {seq1edge,seq2edge}; // If it is possible to coalesce, return it.
                }
            }
            li.remove(); // This is no longer a coalescenceedge because it cannot coalesce with anything.
        }
        return null;
    }
    
    private final boolean possibleToCoalesce(final int[] sequence1, final int[] sequence2, final int[] overlap){
        for (int marker = overlap[0]; marker<=overlap[1]; marker++){
            if (sequence1[marker]!=sequence2[marker] && sequence1[marker]<=1 && sequence2[marker]<=1) return false;
            else if (sequence1[marker]>1 && sequence2[marker]>1 && isThisAConflict(sequence1[marker],sequence2[marker])) return false;
        }
        return true;
    }
    
    private final void makeCoalescence(final Integer[] possiblecoalescence){
        numcoalescences++;
        currentarg.addLast(new ArgStructure(time++,ArgStructure.Type.Co,possiblecoalescence[0],possiblecoalescence[1],nextparent,-1,-1));
        if (VERBOSE) output.writer().println(currentarg.getLast());
        updateCoalescenceEdgesAfterCoalescence(possiblecoalescence[0],possiblecoalescence[1],nextparent);
        updateSequencesAfterCoalescence(possiblecoalescence[0],possiblecoalescence[1],nextparent++);
    }
    
    private final void updateCoalescenceEdgesAfterCoalescence(final Integer childedge1, final Integer childedge2, final int parentedge){
        // Removes childedge1 and childedge2 from the coalescence edges and puts parentedge in their place.
        final ListIterator<Integer> li = coalescenceedges.listIterator();
        li.add(parentedge);
        boolean got1 = false;
        boolean got2 = false;
        while (li.hasNext()){
            final Integer coedge = li.next();
            if (!got1 && childedge1.equals(coedge)) {
                li.remove();
                if (got2) return;
                got1 = true;
            } else if (!got2 && childedge2.equals(coedge)){
                li.remove();
                if (got1) return;
                got2 = true;
            }
        }
    }
    
    private final void updateSequencesAfterCoalescence(Integer childedge1, Integer childedge2, final int parentedge){
        // (1) Update the currentsequences.
        // (2) Update the startends.
        // (3) Update the edgepointers.
        // (4) Resolve conflicts and dependencies.
        // (5) Update allelefrequencies and possible mutations.
        
        // We will recycle *1 and discard *2, so swap these around for max performance.
        int[] startend1 = startends.remove(childedge1);
        int[] startend2 = startends.remove(childedge2);
        if (startend2[1]-startend2[0]>startend1[1]-startend1[0]){
            final Integer temp1 = childedge1; // seq1 is the longest, so swap these around to minimise the amount of work we do.
            childedge1 = childedge2;
            childedge2 = temp1;
            final int[] temp2 = startend1;
            startend1 = startend2;
            startend2 = temp2;
        }
        final int[] sequence1 = currentsequences.remove(childedge1);
        final int[] sequence2 = currentsequences.remove(childedge2);
        final int[] overlap = getOverlap(startend1,startend2);
        
        // Do the left of the overlapping region. Do (1), (2) and (3).
        if (startend1[0]<overlap[0]){
            for (int marker = startend1[0]; marker<overlap[0]; marker++){
                if (sequence1[marker]>1) edgepointers[sequence1[marker]] = parentedge; // Do (3).
            }
        } else {
            for (int marker = startend2[0]; marker<overlap[0]; marker++){
                sequence1[marker] = sequence2[marker]; // Do (1).
                if (sequence1[marker]>1) edgepointers[sequence2[marker]] = parentedge; // Do (3).
            }
            startend1[0] = startend2[0]; // Do (2).
        }
        
        // Do the right of the overlapping region. Do (1), (2) and (3).
        if (startend1[1]>overlap[1]){
            for (int marker = overlap[1]+1; marker<=startend1[1]; marker++){
                if (sequence1[marker]>1) edgepointers[sequence1[marker]] = parentedge; // Do (3).
            }
        } else {
            for (int marker = overlap[1]+1; marker<=startend2[1]; marker++){
                sequence1[marker] = sequence2[marker]; // Do (1).
                if (sequence1[marker]>1) edgepointers[sequence2[marker]] = parentedge; // Do (3).
            }
            startend1[1] = startend2[1]; // Do (2).
        }
        
        startends.put(parentedge,startend1); // Finish doing (2).
        
        // Do the overlapping region.
        for (int marker = overlap[0]; marker<=overlap[1]; marker++){
            if (sequence1[marker]>1){
                // Then we need to update the conflicts and dependencies and edgepointers.
                if (sequence2[marker]>1){ // Sequences 1 and 2 are unphased.
                    updateUnphasedUnphased(sequence1[marker],sequence2[marker],marker,parentedge); // Do (3), (4) and (5).
                } else { // Sequence 1 is unphased, Sequence 2 is phased.
                    updatePhasedUnphased(sequence2[marker],sequence1[marker],marker,childedge1,parentedge); // Do (3), (4) and (5).
                    sequence1[marker] = sequence2[marker];
                }
            } else {
                if (sequence2[marker]>1){ // Sequence 1 is phased, Sequence 2 is unphased.
                    updatePhasedUnphased(sequence1[marker],sequence2[marker],marker,childedge2,parentedge); // Do (3), (4) and (5).
                } else { // Sequence 1 is phased, Sequence 2 is phased.
                    updatePhasedPhased(sequence1[marker],marker,parentedge); // Do (4) and (5).
                }
            }
        }
        
        // Put the parent sequence and startend into the datastructures.
        currentsequences.put(parentedge,sequence1); // Finish doing (1).
    }
    
    private final void updateUnphasedUnphased(final int unphasedchar1, final int unphasedchar2, final int marker, final int parentedge){
        // This can be made faster by keeping the dependencies sorted. This section needs to be recoded.
        
        // (1) Make the union of the dependencies.
        int previousdependency, dependency;
        // Update the dependencies, take the union of the dependencies of unphasedchar1 and unphasedchar2, removing unphasedchar2.
        if (dependencies[unphasedchar2]!=unphasedchar2){ // Then there is work to do.
            previousdependency = unphasedchar1;
            dependency = dependencies[previousdependency];
            // Find the end of the unphasedchar1 set of dependencies.
            while (dependency!=unphasedchar1 && dependency!=unphasedchar2){
                previousdependency = dependency;
                dependency = dependencies[dependency];
            }
            if (dependency==unphasedchar1){
                // previousdependency is the end of the unphasedchar1 set of dependencies.
                dependencies[previousdependency] = dependencies[unphasedchar2];
                previousdependency = dependencies[previousdependency]; // Remove unphasedchar2 from the dependency set.
                dependency = dependencies[previousdependency]; // Merge the 2 dependeny sets together.
                while (dependency!=unphasedchar2){
                    previousdependency = dependency;
                    dependency = dependencies[dependency];
                }
                dependencies[previousdependency] = unphasedchar1; // The completes the merging of the dependency sets.
            } else {
                // dependency==unphasedchar2.
                dependencies[previousdependency] = dependencies[unphasedchar2]; // Remove unphasedhcar2 from the dependency set.
            }
        }
        
        // (2) Update the dependencies of the conflicts.
        
        if (conflicts[unphasedchar1]==0){
            // Then uphasedchar1 is a missing character.
            conflicts[unphasedchar1] = conflicts[unphasedchar2];
        }
        
        if (conflicts[unphasedchar1]!=0){
            if (conflicts[unphasedchar2]!=0){
                previousdependency = conflicts[unphasedchar1];
                dependency = dependencies[previousdependency];
                while (dependency!=conflicts[unphasedchar1] && dependency!=conflicts[unphasedchar2]){
                    previousdependency = dependency;
                    dependency = dependencies[dependency];
                }
                if (dependency==conflicts[unphasedchar1]){
                    dependencies[previousdependency] = conflicts[unphasedchar2];
                    previousdependency = conflicts[unphasedchar2];
                    dependency = dependencies[previousdependency];
                    while (dependency!=conflicts[unphasedchar2]){
                        previousdependency = dependency;
                        dependency = dependencies[dependency];
                    }
                    dependencies[previousdependency] = conflicts[unphasedchar1];
                } // else there is no work to do.
            }
            
            // Update the conflicts for the dependents.
            dependency = dependencies[unphasedchar1];
            while (dependency!=unphasedchar1){
                conflicts[dependency] = conflicts[unphasedchar1];
                dependency = dependencies[dependency];
            }
            
            // Update the conflicts for the conflicts.
            dependency = conflicts[unphasedchar1];
            do {
                conflicts[dependency] = unphasedchar1;
                dependency = dependencies[dependency];
            } while (dependency!=conflicts[unphasedchar1]);
        }
        
        // Update the edgepointers.
        edgepointers[unphasedchar1] = parentedge;
        
        // Update the allelecounts.
        allelecounts[2][marker]--;
    }
    
    private final void updatePhasedUnphased(final int phasedchar, final int unphasedchar, final int marker, final int unphasededge, final int parentedge){
        // Update conflicts and dependencies and allele count and possible mutations.
        // (1) Update the currentsequences.
        // (2) Update snpcout.
        // (3) Update the currentsequences not involved in the coalescence.
        // (4) It may be possible to make some new mutations.
        // (5) Update the edgepointers.
        
        final int conflictchar = (phasedchar+1)%2; // This is what the conflict characters will be assigned to.
        int[] modifythissequence;
        
        // Update the conflicts first.
        final int initialconflict = conflicts[unphasedchar]; // This points to an unphased maker that this conflicts with.
        if (initialconflict!=0){ // If initialconflict==0 then this character has no conflicts, it is a missing character rather than an unphased one.
            // Now loop over the dependencies for the conflict, setting them to conflictchar.
            int conflict = initialconflict;
            do {
                // Modify the dependent character.
                modifythissequence = currentsequences.get(edgepointers[conflict]);
                modifythissequence[marker] = conflictchar;
                // Update the allele counts.
                allelecounts[conflictchar][marker]++;
                allelecounts[2][marker]--;
                // Move onto the next dependent character.
                conflict = dependencies[conflict];
            } while (conflict!=initialconflict); // We have finished looping round.
        }
        
        // Now update the dependencies.
        int dependent = unphasedchar;
        do {
            if (edgepointers[dependent]!=unphasededge){ // Check that this isn't the sequence being coalesced, we modify this in the calling method.
                modifythissequence = currentsequences.get(edgepointers[dependent]);
                modifythissequence[marker] = phasedchar;
                allelecounts[phasedchar][marker]++;
                allelecounts[2][marker]--;
            }
            dependent = dependencies[dependent];
        } while (dependent!=unphasedchar);
        
        // Update the allelecounts and see whether any mutations are possible.
        if (--allelecounts[2][marker]==0 && allelecounts[conflictchar][marker]>0){
            if (allelecounts[conflictchar][marker]==1){
                if (allelecounts[phasedchar][marker]==1 && rand.nextBoolean()) possiblemutations.add(new int[]{parentedge,marker,conflictchar});
                else {
                    int[] startend;
                    // Then we are mutating the other allele, so we have to find the edge on which this allele resides.
                    for (Map.Entry<Integer,int[]> edge : currentsequences.entrySet()){
                        startend = startends.get(edge.getKey());
                        if (marker>=startend[0] && marker<=startend[1] && edge.getValue()[marker]==conflictchar)
                            possiblemutations.add(new int[]{edge.getKey(),marker,phasedchar});
                    }
                }
            } else if (allelecounts[phasedchar][marker]==1) possiblemutations.add(new int[]{parentedge,marker,conflictchar});
        }
    }
    
    private final void updatePhasedPhased(final int phasedchar, final int marker, final int parentedge){
        // Decrement the allele count and check whether it is possible to put in a mutation.
        if (--allelecounts[phasedchar][marker]==1 && allelecounts[2][marker]==0 && allelecounts[(phasedchar+1)%2][marker]!=0){
            // Then it is possible to make a mutation.
            possiblemutations.add(new int[]{parentedge,marker,(phasedchar+1)%2});
        }
    }
    
// ******************************************************** METHODS TO DO RECOMBINATIONS ********************************************************
    
    private final Integer[] doARecombination(){
        int[] sharedsegment; // sharedsegments are of the form {edge1, edge2, start, end}
        if (rand.nextDouble()<=HEURISTICP)
            // Use the heuristic to decide which shared segment to base the recombination on.
            sharedsegment = getLongestSharedSegment();
        else sharedsegment = getAnySharedSegment();
        // Does this shared segment correspond to a crossover or a gene conversion?
        if (sharedsegment[2]==0 || sharedsegment[3]==nmm1) return doCrossover(sharedsegment);
        else return doGeneConversion(sharedsegment);
    }
    
    private final int[] getLongestSharedSegment(){
        // Initialise the shared segment.
        double currentlongest = -1.0;
        int start, end;
        int[] longestsharedsegment = null;
        LinkedList<int[]> longestsharedsegments = new LinkedList<int[]>();
        int[] startend1, startend2, overlap, sequence1, sequence2;
        
        // We will use the startends to loop over. These are permuted elsewhere.
        // Loop over the edges, finding the shared segments.
        for (int edg1 = startendsarray.length; --edg1>=1;){
            startend1 = startendsarray[edg1].getValue();
            sequence1 = currentsequences.get(startendsarray[edg1].getKey());
            for (int edg2 = edg1; --edg2>=0;){
                startend2 = startendsarray[edg2].getValue();
                overlap = getOverlap(startend1,startend2);
                if (overlap==null || distancematrix[overlap[0]][overlap[1]]<=currentlongest) continue;
                sequence2 = currentsequences.get(startendsarray[edg2].getKey());
                
                // Find the longest shared segment from these two sequences.
                start = overlap[0];
                do {
                    if (  !(
                            (sequence1[start]!=sequence2[start] && sequence1[start]<=1 && sequence2[start]<=1)
                            || (sequence1[start]>1 && sequence2[start]>1 && isThisAConflict(sequence1[start],sequence2[start]))
                            )){ // Then the sequences may be coalesced for this character.
                        end = start+1;
                        while (end<=overlap[1] &&
                                !(
                                (sequence1[end]!=sequence2[end] && sequence1[end]<=1 && sequence2[end]<=1)
                                || (sequence1[end]>1 && sequence2[end]>1 && isThisAConflict(sequence1[end],sequence2[end]))
                                )){
                            end++;
                        }
                        if (distancematrix[start][--end]>=currentlongest){
                            if (distancematrix[start][end]>currentlongest){
                                longestsharedsegments.clear();
                                currentlongest = distancematrix[start][end];
                            }
                            if (start==overlap[0]) start = 0;
                            if (end==overlap[1]) end = nmm1;
                            longestsharedsegments.add(new int[]{edg1,edg2,start,end}); // We map edg1 and edg2 onto the live edge at the end.
                        }
                        start = end+1;
                    } else start++;
                } while (start<=overlap[1] && distancematrix[start][overlap[1]]>=currentlongest);
                
                // Store the potentially longest sharedsegment, and return this if it is a longest possible shared segment.
                switch (longestsharedsegments.size()){
                    case 0 : break;
                    case 1 : {
                        longestsharedsegment = longestsharedsegments.remove();
                        if (currentlongest==longestpossible){
                            longestsharedsegment[0] = startendsarray[longestsharedsegment[0]].getKey();
                            longestsharedsegment[1] = startendsarray[longestsharedsegment[1]].getKey();
                            return longestsharedsegment;
                        }
                        break;
                    } default : {
                        longestsharedsegment = longestsharedsegments.get(rand.nextInt(longestsharedsegments.size()));
                        longestsharedsegments.clear();
                        if (currentlongest==longestpossible){
                            longestsharedsegment[0] = startendsarray[longestsharedsegment[0]].getKey();
                            longestsharedsegment[1] = startendsarray[longestsharedsegment[1]].getKey();
                            return longestsharedsegment;
                        }
                    }
                }
            }
        }
        if (longestsharedsegment==null){
        	env.log().println("");
        	//System.err.println();
        }
        
        // Return one of the longest shared segments found.
        longestsharedsegment[0] = startendsarray[longestsharedsegment[0]].getKey();
        longestsharedsegment[1] = startendsarray[longestsharedsegment[1]].getKey();
        return longestsharedsegment;
    }
    
    private final int[] getAnySharedSegment(){
        // Initialise the shared segment.
        int start, end;
        LinkedList<int[]> sharedsegments = new LinkedList<int[]>();
        int[] startend1, startend2, overlap, sequence1, sequence2;
        
        // We will use the startends to loop over. These are permuted elsewhere.
        // Loop over the edges, finding the shared segments.
        for (int edg1 = startendsarray.length; --edg1>=1;){
            startend1 = startendsarray[edg1].getValue();
            sequence1 = currentsequences.get(startendsarray[edg1].getKey());
            for (int edg2 = edg1; --edg2>=0;){
                startend2 = startendsarray[edg2].getValue();
                overlap = getOverlap(startend1,startend2);
                if (overlap==null) continue;
                sequence2 = currentsequences.get(startendsarray[edg2].getKey());
                
                // Get the shared segments.
                start = overlap[0];
                do {
                    if (  !(
                            (sequence1[start]!=sequence2[start] && sequence1[start]<=1 && sequence2[start]<=1)
                            || (sequence1[start]>1 && sequence2[start]>1 && isThisAConflict(sequence1[start],sequence2[start]))
                            )){ // Then the sequences may be coalesced for this character.
                        end = start+1;
                        while (end<=overlap[1] &&
                                !(
                                (sequence1[end]!=sequence2[end] && sequence1[end]<=1 && sequence2[end]<=1)
                                || (sequence1[end]>1 && sequence2[end]>1 && isThisAConflict(sequence1[end],sequence2[end]))
                                )){
                            end++;
                        }
                        if (start==overlap[0]) start = 0;
                        if (end>overlap[1]) end = nmm1;
                        else end--;
                        sharedsegments.add(new int[]{edg1,edg2,start,end}); // We map edg1 and edg2 onto the live edge at the end.
                        start = end+1;
                    } else start++;
                } while (start<=overlap[1]);
                
                // If we have found a shared segment, return it.
                switch (sharedsegments.size()){
                    case 0 : break;
                    case 1 : {
                        final int[] sharedsegment = sharedsegments.getFirst();
                        sharedsegment[0] = startendsarray[sharedsegment[0]].getKey();
                        sharedsegment[1] = startendsarray[sharedsegment[1]].getKey();
                        return sharedsegment;
                    } default : {
                        final int[] sharedsegment = sharedsegments.get(rand.nextInt(sharedsegments.size()));
                        sharedsegment[0] = startendsarray[sharedsegment[0]].getKey();
                        sharedsegment[1] = startendsarray[sharedsegment[1]].getKey();
                        return sharedsegment;
                    }
                }
            }
        }
        return null; // We never get here.
    }
    
    private final Integer[] doCrossover(final int[] sharedsegment){
        int reclocation; // Where we are putting the recombination.
        if (sharedsegment[2]==0) reclocation = sharedsegment[3];
        else reclocation = sharedsegment[2]-1;
        final Integer recedge = sharedsegment[0]; // The edge on which the recombination will be put. These edges are in a random order.
        final Integer coedge = sharedsegment[1]; // This edge will not have a recombination put on it.
        
        // Get the recombination edge
        // We don't need to remove recedge from coedges because coedges is empty when we have to perform a recombination.
        final int[] recedgesequence = currentsequences.remove(recedge); // And this will become the right parent.
        final int[] recedgestartend = startends.remove(recedge);
        
        // Form the left parent.
        final int[] leftparent = new int[nummarkers];
        for (int marker = recedgestartend[0]; marker<=reclocation; marker++){
            leftparent[marker] = recedgesequence[marker];
            if (leftparent[marker]>1) edgepointers[leftparent[marker]] = nextparent;
        }
        coalescenceedges.add(nextparent);
        currentsequences.put(nextparent,leftparent);
        startends.put(nextparent,new int[]{recedgestartend[0],reclocation});
        
        currentarg.add(new ArgStructure(time++,ArgStructure.Type.Re,recedge,-1,nextparent,nextparent+1,reclocation));
        if (VERBOSE) output.writer().println(currentarg.getLast());
        
        // Form the rightparent.
        recedgestartend[0] = reclocation+1;
        nextparent++;
        for (int marker = recedgestartend[0]; marker<=recedgestartend[1]; marker++)
            if (recedgesequence[marker]>1) edgepointers[recedgesequence[marker]] = nextparent;
        startends.put(nextparent,recedgestartend);
        currentsequences.put(nextparent,recedgesequence);
        coalescenceedges.add(nextparent++);
        
        numrecombinations++;
        if (reclocation==sharedsegment[3]) return new Integer[]{coedge,nextparent-2};
        else return new Integer[]{coedge,nextparent-1};
    }
    
    private final Integer[] doGeneConversion(final int[] sharedsegment){
        // Get the sequence that we will put the gene conversion on.
        final int[] recedgesequence = currentsequences.remove(sharedsegment[0]); // The edges in the shared segment are already in a random order.
        final int[] recedgestartend = startends.remove(sharedsegment[0]);
        // We don't need to remove recedge from coedges because coedges is empty when we have to perform a recombination.
        
        // Create the sequence to the left of the first break point.
        final int[] leftsequence = new int[nummarkers];
        final int[] leftstartend = new int[]{recedgestartend[0],sharedsegment[2]-1};
        for (int marker = leftstartend[0]; marker<=leftstartend[1]; marker++){
            leftsequence[marker] = recedgesequence[marker];
            if (leftsequence[marker]>1) edgepointers[leftsequence[marker]] = nextparent;
        }
        // Put this left sequence into the data structures.
        currentarg.add(new ArgStructure(time++,ArgStructure.Type.Re,sharedsegment[0],-1,nextparent,nextparent+1,leftstartend[1]));
        if (VERBOSE) output.writer().println(currentarg.getLast());
        currentsequences.put(nextparent,leftsequence);
        startends.put(nextparent,leftstartend);
        coalescenceedges.add(nextparent);
        nextparent+=2; // To update from the child of the second recombination.
        
        // Create the sequence between the first and second breakpoint.
        final int[] middlesequence = new int[nummarkers];
        final int[] middlestartend = new int[]{sharedsegment[2],sharedsegment[3]};
        for (int marker = middlestartend[0]; marker<=middlestartend[1]; marker++){
            middlesequence[marker] = recedgesequence[marker];
            if (middlesequence[marker]>1) edgepointers[middlesequence[marker]] = nextparent;
        }
        // Put this middle sequence into the data structures.
        currentarg.add(new ArgStructure(time++,ArgStructure.Type.Re,nextparent-1,-1,nextparent,nextparent+1,sharedsegment[3]));
        if (VERBOSE) output.writer().println(currentarg.getLast());
        currentsequences.put(nextparent,middlesequence);
        startends.put(nextparent,middlestartend);
        coalescenceedges.add(nextparent++);
        
        // Recycle the recedge to be the edge to the right of the gene conversion.
        recedgestartend[0] = sharedsegment[3]+1;
        for (int marker = recedgestartend[0]; marker<=recedgestartend[1]; marker++){
            if (recedgesequence[marker]>1) edgepointers[recedgesequence[marker]] = nextparent;
        }
        // Put this right sequence into the data structures.
        currentsequences.put(nextparent,recedgesequence);
        startends.put(nextparent,recedgestartend);
        coalescenceedges.add(nextparent++);
        
        numgeneconversions++;
        return new Integer[]{sharedsegment[1],nextparent-2};
    }
    
// ******************************************************** UTILITIES ********************************************************
    
    private final void permuteStartEnds(){
        final LinkedList<Map.Entry<Integer,int[]>> startendstemp = new LinkedList(startends.entrySet());
        Collections.shuffle(startendstemp);
        startendsarray = startendstemp.toArray(new Map.Entry[startendstemp.size()]);
    }
    
    private static final int[] getOverlap(final int[] startend1, final int[] startend2){
        if (startend1[1]<startend2[0] || startend2[1]<startend1[0]) return null;
        else {
            final int[] overlap = new int[2];
            if (startend1[0]<startend2[0]) overlap[0] = startend2[0];
            else overlap[0] = startend1[0];
            if (startend1[1]>startend2[1]) overlap[1] = startend2[1];
            else overlap[1] = startend1[1];
            return overlap;
        }
    }
    
    private final boolean isThisAConflict(final int unphasedmarker1, final int unphasedmarker2){
        // Optimised version that requires the conflicts to be linked in increasing order.
        final int initialconflict = conflicts[unphasedmarker1];
        int conflict = initialconflict;
        int dependency;
        while (true){
            dependency = dependencies[conflict];
            if (dependency==unphasedmarker2) return true; // unphasedmarker2 is one of the conflicts.
            if (dependency==initialconflict) return false;
            conflict = dependency;
        }
    }
    
    private final boolean isFinished(){
        // To check whether we are finished after performing a coalescence.
        for (int marker = nummarkers; --marker>=0;) if (allelecounts[0][marker]+allelecounts[1][marker]+allelecounts[2][marker]>1) return false;
        return true;
    }
    
    
// ******************************************************** INTERFACE ********************************************************
    
    /**
     * Returns the inferred ARGs.
     *
     * @return  The ARGs.
     */
    public final LinkedList<ArgStructure>[] getArgs(){
        return args;
    }
    
    /**
     * Print ARGs/
     *
     * Prints the ARGs to the terminal.
     */
    public final void printArgs(){
        output.writer().println("%ARGS");
        output.writer().println("TIME OPERATION CHILD1 {CHILD2} PARENT1 {PARENT2} {LOCATION}");
        for (int arg = 0; arg<args.length; arg++){
            output.writer().println("ARG " + arg);
            for (ArgStructure struct : args[arg]){
                output.writer().println(struct);
            }
        }
    }
    
    /**
     * Print marginal trees.
     *
     * Prints all marginal trees to the terminal.
     */
    public final void printTrees(){
        output.writer().println("%TREES");
        output.writer().println("TIME CHILD1 CHILD2 PARENT");
        for (int arg = 0; arg<args.length; arg++)
            for (int marker = 0; marker<nummarkers; marker++)
                printTree(marker,arg);
    }
    
    /**
     * Print one marginal tree.
     *
     * Prints one marginal tree to the terminal.
     */
    public final void printTree(final int marker, final int whicharg){
        output.writer().println("TREE: ARG " + whicharg + " MARKER " + marker);
        final Hashtable<Integer,Short> argtreemap = new Hashtable<Integer,Short>();
        final short[][] marginaltree = new short[3][numsequences-1];
        final LinkedList<ArgStructure> arg = args[whicharg]; // Iterate over the ARGs.
        for (int i = numsequences; --i>=0;) argtreemap.put(i,(short)i);
        short parent = (short)numsequences;
        int node = 0;
        boolean firstone = true;
        for (ArgStructure a : arg){
            switch (a.t){
                case Mu:{
                    if (argtreemap.containsKey(a.child1)) argtreemap.put(a.parent1, argtreemap.remove(a.child1));
                    break;
                } case Co:{
                    if (argtreemap.containsKey(a.child1) && argtreemap.containsKey(a.child2)){
                        marginaltree[0][node] = argtreemap.remove(a.child1);
                        marginaltree[1][node] = argtreemap.remove(a.child2);
                        marginaltree[2][node] = parent;
                        output.writer().println(marginaltree[0][node] + " " + marginaltree[1][node] + " " + marginaltree[2][node]);
                        node++;
                        argtreemap.put(a.parent1, parent++);
                    } else if (argtreemap.containsKey(a.child1)){argtreemap.put(a.parent1, argtreemap.remove(a.child1));
                    } else if (argtreemap.containsKey(a.child2)){argtreemap.put(a.parent1, argtreemap.remove(a.child2));}
                    break;
                } case Re:{
                    if (argtreemap.containsKey(a.child1)){
                        if (marker<=a.location) argtreemap.put(a.parent1, argtreemap.remove(a.child1));
                        else argtreemap.put(a.parent2, argtreemap.remove(a.child1));
                    }
                }
            }
        }
        argtreemap.clear();
    }    
}

