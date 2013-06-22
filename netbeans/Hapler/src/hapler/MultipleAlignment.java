/*
 *  Copyright 2011 Shawn Thomas O'Neil
 *
 *  This file is part of Hapler.
 *
 *  Hapler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Hapler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Hapler.  If not, see <http://www.gnu.org/licenses/>.
 */

package hapler;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.jgrapht.graph.SimpleGraph;

/**
 *
 * @author soneil
 * Note: There are a number of invariants which we need to ensure in this class
 * 1) Every time a sequence is added or removed, we may need to update the multiple alignment length
 *     (we don't recompute it every time it's asked for because that would be slow)
 * 2) Every time a sequence is added, we have to ensure that the snps covering it are added to it's list
 *      we also have to make sure the sequence is added to all the SNPs it covers
 * 3) Every time a SNP is added, we have to ensure that it is added to sequences covering it
 *      we also have to make sure the sequences covering it are added to the SNPs list
 * 4) Every time a sequence is removed, we have to remove it from all the snps that cover it
 *      we may as well remove the snps from the seq as well (but we never remove them from the master list of SNPs)
 *      This keeps an invariant that for every sequence in the mult. align, it knows about all the snps,
 *      and for every snp in the multiple align, it knows about every sequence (and only those seqs) in the alignm.'
 *      (So would leaving them I guess, but we'll just assume that seqs outside a mult align. don't know about SNPs)
 * IGNORE 5) every time a SNP is removed, we have to remove it from all the sequences which cover it
 *
 * DONE: What if removing a sequence makes a SNP no longer a SNP by some rule? So what. Leave it.
 * DONE: What if removing a SNP makes a sequence be contained within another? You can't delete SNPS.
 *
 * Well, how about we don't allow SNPS to be removed. They were proper snps at one point, just ignoring some
 *   reads doesn't mean they're not variant positions anymore.
 */
public class MultipleAlignment {
	
	private String name;
	private String givenConsensus;
	private ArrayList<Sequence> sequences;
	private HashSet<SNP> snps;
	private ArrayList<HaplotypeBlock> haploBlocks;
	private Haplotype universalHaplotype;
	private int length;  // number of columns in the multiple alignment
	private int numUnmaskedSequences;
	private int numSNPs;
	private int maxSeqNameLength;
	private HashSet<Character> symbolDictionary;




	public MultipleAlignment() {
		name = "";
		givenConsensus = "";
		sequences = new ArrayList<Sequence>();
		snps = new HashSet<SNP>();
		haploBlocks = new ArrayList<HaplotypeBlock>();
		length = -1;
		numUnmaskedSequences = 0;
		numSNPs = 0;
		maxSeqNameLength = 0;
		symbolDictionary = new HashSet<Character>();
	}


	public String SNPsString() {
		StringBuilder sb = new StringBuilder();

		for(int i = 0; i < length; i++) {
			sb.append('_');
		}

		for(SNP snp: snps) {
		//for(int i = 0; i < snps.size(); i++) {
			//SNP snp = snps.get(i);
			//if(snp.isMasked()) {
			//	sb.replace(snp.getPosition(), snp.getPosition()+1, "+");
			//}
			//else {
				sb.replace(snp.getPosition(), snp.getPosition()+1, "*");
			//}
		}

		return sb.toString();
	}

	public String SNPsAsString() {
		StringBuilder sb = new StringBuilder();
		int snpsCount = snps.size();
		if(snps.size() != 0) {
			int count = 1;
			for(SNP snpi : snps)
			//for(int i = 0; i < snps.size()-1; i++) {
				//SNP snpi = snps.get(i);
				if(count != snpsCount) sb.append(snpi.getPosition() + ",");
				else sb.append(snpi.getPosition());
				count = count + 1;
			}
		return sb.toString();
	}

	public String toString() {
		Collections.sort(sequences, new SequenceStartPosLengthComparator());
		StringBuffer sb = new StringBuffer();
		sb.append("Num Sequences: " + sequences.size() + System.getProperty("line.separator"));
		sb.append("Num SNPS: " + snps.size() + System.getProperty("line.separator"));
		sb.append("Num Unmasked Sequences: " + numUnmaskedSequences + System.getProperty("line.separator"));
		sb.append("Num Unmasked SNPS: " + numSNPs + System.getProperty("line.separator"));

		

		char[] snpCharArray = new char[length+maxSeqNameLength+5];
		for(int i = 0; i < snpCharArray.length; i++) {
			snpCharArray[i] = ' ';
		}
		for(SNP snp : snps) {
		//for(int i = 0; i < snps.size(); i++) {
			//SNP snp = snps.get(i);
			//if(!snp.isMasked()){
			//	snpCharArray[snp.getPosition()+maxSeqNameLength+5] = '*';
			//}
			//else {
				snpCharArray[snp.getPosition()+maxSeqNameLength+5] = '+';
			//}
		}
		//for(char c : snpCharArray) {
		for(int i = 0; i < snpCharArray.length; i++) {
			char c = snpCharArray[i];
			sb.append(c);
		}
		sb.append(System.getProperty("line.separator"));

		// Append the unmasked sequences
		//for(Sequence seq : sequences) {
		for(int i = 0; i < sequences.size(); i++) {
			Sequence seq = sequences.get(i);
			if(!seq.isMasked() || true) {
				sb.append(String.format("%1$" + (maxSeqNameLength +1) + "s    ", seq.getName()));
				if(!seq.isMasked()) { sb.append( seq.sequenceToString(false) + System.getProperty("line.separator"));}
				else {sb.append( seq.sequenceToStringReverseCase(false) + System.getProperty("line.separator"));}
				// print the seqs this one is masking below it, reversing the case?
				ArrayList<Sequence> masking = seq.getMasking();

			}
		}

		return sb.toString();
	}

	/**
	 * Returns an array of chars representing the multiple alignment of this position
	 * _of unmasked sequences_
	 * no particular order is gauranteed. ~ represents a read not covering the position.
	 * Yes, it is currently possible to get a column of the alignment that is way off the end of the
	 * alignment, which will consist of sequences.length() '~'s.
	 * @param pos
	 * @return
	 */
	public char[] unmaskedColumn(int pos) {
		char[] column = new char[numUnmaskedSequences];
		int i = 0;
		//for(Sequence seq : sequences) {
		for(int j = 0; j < sequences.size(); j++) {
			Sequence seq = sequences.get(j);
			if(!seq.isMasked()) {
				column[i] = seq.baseAtAlignmentPosition(pos);
				i = i + 1;
			}
		}
		return column;
	}

	public void maskSequence(Sequence seq) {
		// This is an expensive check...
		//if(sequences.contains(seq)) {
			if(seq.setMasked(true)) numUnmaskedSequences = numUnmaskedSequences - 1;
		//}
	}

	public void unmaskSequence(Sequence seq) {
		// This is an expensive check...
		//if(sequences.contains(seq)) {
			if(seq.setMasked(false)) numUnmaskedSequences = numUnmaskedSequences + 1;
		//}
	}


	/**
	 * returns the majority vote consensus of the overall, you know, thingy
	 * amongst ALL masked and non-masked sequences
	 * @return
	 */
	public String majorityVoteConsensus() {
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < this.length; i++) {
			sb.append(this.majorityVoteConsensusAtPosition(i));
		}

		return sb.toString();
	}


	/**
	 * returns the majority vote, amongst ALL masked and non-masked sequences in the
	 * alignment.
	 * @param pos
	 * @return
	 */
	public char majorityVoteConsensusAtPosition(int pos) {
		HashMap<Character, Integer> baseCountHash = new HashMap<Character, Integer>();
		//for(Sequence seq: haplotype) {
		for(int j = 0; j < sequences.size(); j++) {
			Sequence seq = sequences.get(j);
			//if(i == 0) System.out.println("Considering sequence " + seq.getName() + " masked is: " + seq.isMasked());
			//if(seq.isMasked()) System.out.println("this shouldnt be masked...");
			char baseAti = seq.baseAtAlignmentPosition(pos);
			if(!baseCountHash.containsKey(baseAti)) {
				baseCountHash.put(baseAti, 1);
			}
			else {
				baseCountHash.put(baseAti, baseCountHash.get(baseAti)+1);
			}
		}
		int maxCount = -1;
		char maxBase = '~';
		for(Character base: baseCountHash.keySet()) {
			int count = baseCountHash.get(base);
			if(count > maxCount && base != '~') {
				maxCount = count;
				maxBase = base;
			}
		}
		return maxBase;
	}



	/**
	 * Returns an array of strings representing the multiple alignment of positions
	 * startPos to endPos _of unmasked sequences_.
	 * No particular order is gauranteed. ~ represents a read not
	 * covering a position.
	 * It is currently possible to get an alignment off the end of the alignment;
	 * these position will consist of '~'.
	 * @param startPos
	 * @param endPos
	 * @return
	 */
	public String[] unmaskedColumns(int startPos, int endPos) {
		String[] columns = new String[numUnmaskedSequences];
		int i = 0;
		//for(Sequence seq: sequences) {
		for(int j = 0; j < sequences.size(); j++) {
			Sequence seq = sequences.get(j);
			if(!seq.isMasked()) {
				columns[i] = seq.basesFromAlignmentPositionTo(startPos, endPos);
				i = i + 1;
			}
		}
		return columns;
	}

	/*
	 * Represents a column in the alignment (including masked sequences) as a hash from allele-> count, not including ~ positions
	 */
	public HashMap<Character, Integer> variantsHash(int pos) {
		HashMap<Character, Integer> variantCounts = new HashMap<Character, Integer>();
		for(int j = 0; j < sequences.size(); j++) {
			Sequence seq = sequences.get(j);
			char charAtPos = seq.baseAtAlignmentPosition(pos);
			if(charAtPos != '~') {
				if(!variantCounts.containsKey(charAtPos)) {
					variantCounts.put(charAtPos, 1);
				}
				else {
					int count = variantCounts.get(charAtPos);
					variantCounts.put(charAtPos,count + 1);
				}
			}
		}
		return variantCounts;
	}

	/**
	 * Determines whether there is any actual variation at a position, amongst unmasked and masked sequences.
	 * IE, given all bases covering position, are the non-~ characters all matching?
	 * @return
	 */
	public boolean containsVariationAtPosition(int pos) {
		String firstBaseSeen = "~";

		for(int i = 0; i < sequences.size(); i++) {
			Sequence seq = sequences.get(i);
			String seqAtPos = seq.basesFromAlignmentPositionTo(pos, pos);
			if(!seqAtPos.equals("~")) {
				if(firstBaseSeen.equals("~")) {
					firstBaseSeen = seqAtPos;
				}
				else {
					if(!firstBaseSeen.equals(seqAtPos)) return true;
				}
			}
		}

		return false;
	}

	
	/**
	 * WARNING: This method allows you to add a sequence twice. That wouldn't be smart.
	 * We have to make sure to update the multiple alignment length, if necessary
	 * And we have to make sure to add all the SNPs which the new read cover to that read
	 * And we have to make sure to add this read to all the SNPs which it covers
	 * @param seq
	 */
	public void addSequence(Sequence seq) {
		//if(!sequences.contains(seq)) {
			char[] charSeq = seq.sequenceToString(true).toCharArray();
			for(int i = 0; i < charSeq.length; i++) {
				if(charSeq[i] != '~') {
					symbolDictionary.add(charSeq[i]);
				}
			}
			sequences.add(seq);
			seq.setAlignment(this);
			if(maxSeqNameLength < seq.getName().length()) {
				maxSeqNameLength = seq.getName().length();
			}
			if(length < seq.endPosition() + 1) {
				length = seq.endPosition() + 1;
			}
			for(SNP snp : snps) {
			//for(int j = 0; j < snps.size(); j++) {
				//SNP snp = snps.get(j);
				if(snp.getPosition() <= seq.endPosition() && snp.getPosition() >= seq.getStartPosition()) {
					try {
						snp.addCoveredSequence(seq);
						seq.addCoveredSNP(snp);
					}
					catch(Exception e){
						System.err.println("#############################");
						System.err.println("Error: " + e.getMessage());
						System.err.println("#############################");
						System.err.println();
						e.printStackTrace();
					}
				}
			}
			numUnmaskedSequences = numUnmaskedSequences + 1;
		//}
	}


	/**
	 * Compares each unmasked sequence against each unmasked sequence. If sequence a is covered by sequence
	 * b, and they do not conflict at any unmasked SNP positions, then a will be masked. This ensures
	 * that any conflict graph drawn on the sequences is a co-comparability graph.
	 *
	 * Further, during this process each sequence i which is masked by j is added to j's
	 * list of masked sequences. This is useful for determing where sequences "went".
	 *
	 * This process is done according to start position, and then length, this maximally groups
	 * masked sequences together into what masks them. 
	 */
	public void maskUnmaskedRedundantSeqsByUnmasked() throws Exception {
		// If the alignment is _just_ two sequences of the same length, they should not both be masked.
		HashSet<Sequence> toMask = new HashSet<Sequence>();

		ArrayList<Sequence> unMaskedSeqs = new ArrayList<Sequence>();
		//for(Sequence seq : sequences) {
		for(int j = 0; j < sequences.size(); j++) {
			Sequence seq = sequences.get(j);
			if(!seq.isMasked()) {
				unMaskedSeqs.add(seq);
			}
		}
		// O(n^2*)
		Collections.sort(unMaskedSeqs, new SequenceFirstSNPLengthComparator());
		int numUnmasked = unMaskedSeqs.size();
		for(int i = 0; i < numUnmasked-1; i++) {
			for(int j = i+1; j < numUnmasked; j++) {
				// TODO: seq i should mask j if seq j's SNPs are a subset of i's (and they don't conflict), --- but do they need to overlap?
				// What do we do about reads which don't cover _any_ snps?
				////if(unMaskedSeqs.get(i).covers(unMaskedSeqs.get(j)) && !unMaskedSeqs.get(i).conflictsWithAtSNP(unMaskedSeqs.get(j)) && !toMask.contains(unMaskedSeqs.get(j))) {
				//if(unMaskedSeqs.get(j).agreesWithAtCoveredSNPs(unMaskedSeqs.get(i)) && unMaskedSeqs.get(j).overlapsAtCoveredSNPs(unMaskedSeqs.get(i)) && !toMask.contains(unMaskedSeqs.get(j))) {
				if(!toMask.contains(unMaskedSeqs.get(j)) && unMaskedSeqs.get(i).shouldMask(unMaskedSeqs.get(j))) {
					toMask.add(unMaskedSeqs.get(j));
					unMaskedSeqs.get(i).addMaskedSequence(unMaskedSeqs.get(j));
					//System.out.println("Seq " + unMaskedSeqs.get(i).getName() + " masks " + unMaskedSeqs.get(j).getName());
				}
			}
		}

		maskSequences(toMask);
	}

	/**
	 * Masks all currently unmasked Sequences that do not cover any SNP position, and returns the
	 * list of those sequences that were masked.
	 * @return
	 */
	public HashSet<Sequence> maskUnmaskedNoSNPSequences(){
		HashSet<Sequence> masking = new HashSet<Sequence>();
		for(int i = 0; i < sequences.size(); i++){
			Sequence seq = sequences.get(i);
			if(!seq.isMasked() && seq.getSNPsCovered().size() == 0){
				masking.add(seq);
			}
		}
		maskSequences(masking);
		return masking;
	}


	/**
	 * Batch masking of sequences.
	 * @param seqs
	 */
	public void maskSequences(HashSet<Sequence> seqs) {
		for(Sequence seq : seqs) {
		//for(int j = 0; j < seqs.size(); j++) {
			//Sequence seq = seqs.get(j);
			maskSequence(seq);
		}
	}

	/**
	 * Returns a represtation of the connnected components of the conflict graph of
	 * unmasked sequences. O(n^2*maxSNPcountPerSeq*setUnionTime)
	 * @return
	 */
	public ArrayList<ArrayList<Sequence>> unmaskedConflictGraphConnectedComponents() {
		ArrayList<ArrayList<Sequence>> ccs = new ArrayList<ArrayList<Sequence>>();
		ArrayList<DisjointSetNode> djNodes = new ArrayList<DisjointSetNode>();
		//for(Sequence seq: sequences) {
		// O(n)
		for(int j = 0; j < sequences.size(); j++) {
			Sequence seq = sequences.get(j);
			if(!seq.isMasked()) {
				djNodes.add(new DisjointSetNode(seq));
			}
		}

		// union the connected components, O(n^2*maxSNPcountPerSeq*setUnionTime)
		//for(DisjointSetNode inode : djNodes) {
		for(int i = 0; i < djNodes.size()-1; i++) {
			DisjointSetNode inode = djNodes.get(i);
			//for(DisjointSetNode jnode: djNodes) {
			for(int j = i+1; j < djNodes.size(); j++) {
				DisjointSetNode jnode = djNodes.get(j);
				Sequence iseq = (Sequence)inode.getData();
				Sequence jseq = (Sequence)jnode.getData();
				if(iseq != jseq && iseq.conflictsWithAtSNP(jseq)) {
					inode.union(jnode);
				}
			}
		}

		// Create a map from root nodes to lists of sequences
		//O(n)
		HashMap<DisjointSetNode, ArrayList<Sequence>> rootsToLists = new HashMap<DisjointSetNode, ArrayList<Sequence>>();
		//for(DisjointSetNode inode: djNodes) {
		for(int i = 0; i < djNodes.size(); i++) {
			DisjointSetNode inode = djNodes.get(i);
			DisjointSetNode root = inode.find();
			if(!rootsToLists.containsKey(root)){
				ArrayList<Sequence> newList = new ArrayList<Sequence>();
				Sequence theSeq = (Sequence)inode.getData();
				newList.add(theSeq);
				rootsToLists.put(root, newList);
			}
			else
			{
				Sequence theSeq = (Sequence)inode.getData();
				ArrayList<Sequence> theList = rootsToLists.get(root);
				theList.add(theSeq);
			}
		}

		// Now translate the list of lists to the data structure we are returning
		//ArrayList<Sequence>[] rootsToListValues = (ArrayList<Sequence>[])rootsToLists.values().toArray();
		for(ArrayList<Sequence> theList :rootsToLists.values()) {
		//for(int i = 0; i < rootsToListValues.length; i++) {
			//ArrayList<Sequence> theList = rootsToListValues[i];
			ccs.add(theList);
		}
		
		return ccs;
	}


	/**
	 * Returns a graph, where each node is an unmasked sequence. Any two nodes whose
	 * sequences conflict at an unmasked snp have an edge between them.
	 * @return
	 */
	public SimpleGraph<GraphNode, DefaultEdge> conflictGraph() {
		SimpleGraph<GraphNode, DefaultEdge> conflictGraph = new SimpleGraph<GraphNode,DefaultEdge>(DefaultEdge.class);
		ArrayList<GraphNode> graphNodes = new ArrayList<GraphNode>();


		//for(Sequence seq : sequences) {
		for(int i = 0; i < sequences.size(); i++) {
			Sequence seq = sequences.get(i);
			if(!seq.isMasked()) {
				GraphNode newNode = new GraphNode(seq);
				graphNodes.add(newNode);
				conflictGraph.addVertex(newNode);
			}
		}


		for(int i = 0; i < graphNodes.size() - 1; i++) {
			for(int j = i + 1; j < graphNodes.size(); j++) {
				GraphNode inode = graphNodes.get(i);
				GraphNode jnode = graphNodes.get(j);
				Sequence iseq = (Sequence)inode.getData();
				Sequence jseq = (Sequence)jnode.getData();

				if(iseq != jseq && iseq.conflictsWithAtSNP(jseq)) {
					conflictGraph.addEdge(inode, jnode);
				}
			}
		}

		return conflictGraph;
	}


	/**
	 * Given a list of sequences, creates a directed, transitively reduced,
	 * "positive association" graph. In this graph, node i connects to node j
	 * if sequence i is before j, i overlaps j, i and j do not conflict at any SNP,
	 * and i and j share at least one snp allele. 
	 * @param sequences
	 * @return
	 */
	public SimpleDirectedGraph<GraphNode, DefaultEdge> positiveTRGraphFromSequences(ArrayList<Sequence> sequences) {
		System.out.println("# COMPUTING POSITIVE TR GRAPH");
		Collections.sort(sequences, new SequenceStartPosLengthComparator());

		SimpleDirectedGraph<GraphNode, DefaultEdge> newGraph = new SimpleDirectedGraph<GraphNode, DefaultEdge>(DefaultEdge.class);


		// Let's keep a list of ancestors for each sequence so that we can make sure our
		// graph is transitively reduced. This will use O(n^2) memory, but should be faast.
		HashMap<Sequence,HashSet<Sequence>> ancestorSet = new HashMap<Sequence,HashSet<Sequence>>(sequences.size());
		for(int j = 0; j < sequences.size(); j++) {
			ancestorSet.put(sequences.get(j), new HashSet<Sequence>(sequences.size()));
		}

		// We will also need a mapping from sequences to graphnodes containing them
		HashMap<Sequence, GraphNode> seqsToNodes = new HashMap<Sequence, GraphNode>();
		for(int j = 0; j < sequences.size(); j++) {
			Sequence seqj = sequences.get(j);
			GraphNode newNode = new GraphNode(seqj);
			newGraph.addVertex(newNode);
			seqsToNodes.put(seqj, newNode);
		}

		// Here's the sorted set (sorted by end position) for the current coverage list
		LinkedList<Sequence> coveredQueue = new LinkedList();
		LinkedList<Sequence> tempQueue = new LinkedList();

		for(int j = 0; j < sequences.size(); j++) {
			Sequence seqj = sequences.get(j);
			//System.out.println("#DOT \"" + seqj.getName() + "\" [ pos = \"0," + -1*seqj.getStartPosition() + "\" ];");
			System.err.println("#DOT \"" + seqj.getName()+"XXX" + "\" ");
			//if(j < sequences.size() - 1) {
				//System.out.println("#DOT \"" + seqj.getName() + "\" -> \"" + sequences.get(j+1).getName() + "\" [style=invis,len="+ (sequences.get(j+1).getStartPosition() - seqj.getStartPosition()) + "]");
			//}
			tempQueue.addFirst(seqj);
			//System.out.println("Adding seq starting at position: " + seqj.getStartPosition() + " Ending at: "+ seqj.getEndPosition());
			while(!coveredQueue.isEmpty()) {
				Sequence popped = coveredQueue.removeLast();
				if(!(popped.endPosition() < seqj.getStartPosition())) {
					//System.out.println("   -Popped and adding seq starting at position: " + popped.getStartPosition() + " Ending at: "+ popped.getEndPosition());
					tempQueue.addFirst(popped);
					//if(popped.agreesWithAtCoveredSNPs(seqj)) { // just fucking around ...
					if(popped.sharesSNPWith(seqj) && popped.agreesWithAtCoveredSNPs(seqj)) {
						if(!ancestorSet.get(seqj).contains(popped)) {
							ancestorSet.get(seqj).add(popped);
							ancestorSet.get(seqj).addAll(ancestorSet.get(popped));
							newGraph.addEdge(seqsToNodes.get(popped), seqsToNodes.get(seqj));
							System.err.println("#DOT \"" + popped.getName()+"XXX" + "\" -> \"" + seqj.getName()+"XXX" + "\"");
						}
					}
					//else {
					//	System.out.println("found a disagreement!");
					//}
				}
			}
			LinkedList ttq = coveredQueue;
			coveredQueue = tempQueue;
			tempQueue = ttq;
		}



		return newGraph;
		}



	/**
	 * Given a list of sequences, creates a directed compatibility graph.
	 * In this graph, node i is compatible with j (there is an edge from i to j)
	 * if i starts before j, and they don't conflict at an unmasked SNP.
	 * @param sequences
	 * @return
	 */
	public SimpleDirectedGraph<GraphNode, DefaultEdge> compatabilityGraphFromSequences(ArrayList<Sequence> sequences) {
		Collections.sort(sequences, new SequenceStartPosLengthComparator());
		SimpleDirectedGraph<GraphNode, DefaultEdge> newGraph = new SimpleDirectedGraph<GraphNode, DefaultEdge>(DefaultEdge.class);
		ArrayList<GraphNode> nodes = new ArrayList<GraphNode>();
		for(int i = 0; i < sequences.size(); i++) {
			Sequence seq = sequences.get(i);
			GraphNode newNode = new GraphNode(seq);
			newGraph.addVertex(newNode);
			nodes.add(newNode);
		}
		for(int i = 0; i < nodes.size(); i++) {
			GraphNode nodei = nodes.get(i);
			Sequence seqi = (Sequence)nodei.getData();
			for(int j = 0; j < nodes.size(); j++ ) {
				GraphNode nodej = nodes.get(j);
				Sequence seqj = (Sequence)nodej.getData();
				if(i != j) {
					if(seqi.startsBefore(seqj) && !seqi.conflictsWithAtSNP(seqj)) {
						newGraph.addEdge(nodei, nodej);
						System.err.println("#DOT \"" + seqi.getName() + "\" -> \"" + seqj.getName() + "\"");

					}
				}
			}
		}
		return newGraph;
	}






	
	/**
	 * Do nothing if a snp already exists at this position
	 * We make sure that every read covering this snp is added to the SNPs list
	 * We also make sure this SNP is added to every read covering it
	 * Runs in Theta(numReads)
	 * @param theSNP
	 */
	public void addSNP(SNP theSNP) throws Exception {
		if(!snps.contains(theSNP)) {
			snps.add(theSNP);
			//for(Sequence seq : sequences) {
			for(int i = 0; i < sequences.size(); i++) {
				Sequence seq = sequences.get(i);
				if(seq.containsPosition(theSNP.getPosition())) {
					//try {
						theSNP.addCoveredSequence(seq);
						seq.addCoveredSNP(theSNP);
					//}
					/*catch(Exception e) {
						System.err.println("#############################");
						System.err.println("Error: " + e.getMessage());
						System.err.println("#############################");
						System.err.println();
						e.printStackTrace();
					}*/
				}
			}
			numSNPs = numSNPs + 1;
		}

	}

	/**
	 * Adds an ArrayList<SNP> one at a time.
	 * Runs in time Theta(numSNPs * numReads)
	 * @param snpList
	 */
	public void addSNPs(ArrayList<SNP> snpList) throws Exception {
		//for(SNP snp : snpList) {
		for(int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			addSNP(snp);
		}
	}


	public String getGivenConsensus() {
		return givenConsensus;
	}

	public void setGivenConsensus(String givenConsensus) {
		this.givenConsensus = givenConsensus;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}


	public int getNumUnmaskedSequences() {
		return numUnmaskedSequences;
	}

	public ArrayList<Sequence> getSequences() {
		return sequences;
	}

	public int numSNPs() {
		return numSNPs;
	}

	public HashSet<SNP> getSNPs() {
		/*ArrayList<SNP> toRet = new ArrayList<SNP>();
		for(int i = 0; i < snps.size(); i++) {
			SNP theSNP = snps.get(i);
			//if(!theSNP.isMasked()) {
				toRet.add(theSNP);
			//}
		}
		return toRet;*/
		return snps;
	}

	public int getLength() {
		return length;
	}

	public int getMaxSeqNameLength() {
		return maxSeqNameLength;
	}


	public ArrayList<HaplotypeBlock> getHaploBlocks() {
		return haploBlocks;
	}

	public void setHaploBlocks(ArrayList<HaplotypeBlock> haploBlocks) {
		this.haploBlocks = haploBlocks;
	}

	public void addHaploBlock(HaplotypeBlock haploBlock) {
		haploBlocks.add(haploBlock);
	}

	public Haplotype getUniversalHaplotype() {
		return universalHaplotype;
	}

	public void setUniversalHaplotype(Haplotype universalHaplotype) {
		this.universalHaplotype = universalHaplotype;
	}

	public HashSet<Character> getSymbolDictionary() {
		return symbolDictionary;
	}
}
