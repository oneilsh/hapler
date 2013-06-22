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
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

/**
 * Greedy alg refs: http://books.google.com/books?id=RrtXSKMAmWgC&pg=PA246&lpg=PA246&dq=%22
 * perfect+ordering%22+%22comparability+graph%22&source=bl&ots=yGo1F84zIT&sig=S12WLzpTSFM
 * 7BaWJnY3uGD1FR_c&hl=en&ei=fxDgTLOlCtGgngfEws20Dw&sa=X&oi=book_result&ct=result&resnum=
 * 1&ved=0CBMQ6AEwAA#v=onepage&q=%22perfect%20ordering%22%20%22comparability%20graph%22&f=false
 * http://en.wikipedia.org/wiki/Linear_extension
 * AW HELL. These are for coloring, not clique cover!
 */

/**
 *
 * @author soneil
 */
public class HaplotyperCompat {

	protected MultipleAlignment alignment;
	//private ArrayList<HaplotypeBlock> haplotypeBlocks;
	protected boolean useEpsilon;
	protected int numReps, maxReps;
	protected SimpleWeightedGraph<Sequence, DefaultWeightedEdge> togetherCounts;



	public HaplotyperCompat(MultipleAlignment theAlignment, boolean theUseEpsilon, int theNumReps, int theMaxReps) throws Exception {
		alignment = theAlignment;
		useEpsilon = theUseEpsilon;
		numReps = theNumReps;
		// togetherCounts is a graph of all the sequences in this block, edges between them count how many times sequences showed up in the same haplotype.
		togetherCounts = new SimpleWeightedGraph<Sequence, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		maxReps = theMaxReps;
	}



	/**
	 * Given a multiple alignment of sequences and identified snp positions, first
	 * masks redundant sequences, then for each haplotype block (identified connected
	 * components in the conflict graph), identifies the sequences belonging to each
	 * haplotype. This is done by first computing a TRANSITIVE REDUCTION of the
	 * AGREEMENT GRAPH, then using the path cover algorithm in repetition.
	 * @throws Exception
	 */
	public void executeTR() throws Exception {
		// See notes for execute, below
		
		}


	/**
	 * Given a multiple alignment of sequences and identified snp positions,
	 * first masks redundant sequences, then for each haplotype block (identified
	 * connected componenets in the conflict graph), identifies the sequences belonging
	 * to each haplotype.
	 */
	public void execute() throws Exception {
		// Seqs which contain a SNP _must_ conflict with something. The parser should catch this, if not, we'll throw an exception below.
		// Seqs which don't contain a SNP must be univeral, and will go into the universal haplotype.
		// These sequences are masked, but they aren't maked _by_ anything.
		//TODO: what about user-defined SNPs that don't actually cover any variant positions? The parser should probably catch this...
		HashSet<Sequence> universalSeqs = alignment.maskUnmaskedNoSNPSequences();
		Haplotype universalHap = this.createUniversalHaplotype(universalSeqs);
		HaplotypeBlock universalBlock = this.createUniversalBlock(universalHap);
		alignment.setUniversalHaplotype(universalHap);
		alignment.addHaploBlock(universalBlock);

		// Mask the rest of the sequences that need to be masked
		System.err.println(alignment.getName() + ": Masking redundant sequences...");

		alignment.maskUnmaskedRedundantSeqsByUnmasked();

		System.err.println(alignment.getName() + ": Masking of redundant seqs completed. Computing haplotype blocks.");

		ArrayList<ArrayList<Sequence>> connectedComponents = alignment.unmaskedConflictGraphConnectedComponents();

		System.err.println(alignment.getName() + ": Blocks computed. Processing alignment " + alignment.getName() + ", " + connectedComponents.size() + " haplotype blocks...");

		// After we make sure that all the additional small edges are used,
		// Each CC defines a haplotype block
		Integer hapBlockNum = 1;
		Integer minColors = -1;

		System.err.println("#DOT strict digraph G{");

		for (int i = 0; i < connectedComponents.size(); i++) {
			ArrayList<Sequence> ccSeqList = connectedComponents.get(i);
			System.err.println(alignment.getName() + ": Processing Haplotype Block " + hapBlockNum + ", number non-redundant sequences: " + ccSeqList.size());
			//if(true) throw new Exception("Debug exception...");
			// Moved the sort into the graph computation.
			//Collections.sort(ccSeqList, new SequenceStartPosLengthComparator());


			for (int j = 0; j < ccSeqList.size(); j++) {
				Sequence seq = ccSeqList.get(j);
				togetherCounts.addVertex(seq);
			}

			// We do the bipartite matching a number of times, shuffling the input after each one
			// Create the bipartite representation
			WeightedBipartiteGraph bpGraph = new WeightedBipartiteGraph(alignment.compatabilityGraphFromSequences(ccSeqList), useEpsilon);

			int repsCompleted = 0; int repsUnchanged = 0; int numCCs = 0;
			if(numReps == -1) { // If the numReps is -1, we want "auto," ie., run until we get 10% repetitions with no change in the number of haplotypes, or 10 reps, whichever is larger
										// but don't do more than maxReps...
				while(repsCompleted <= maxReps && (repsCompleted < 20 || repsUnchanged < 0.5*repsCompleted)) {

					ArrayList<ArrayList<Sequence>> haplotypes = bipartiteMatchAndGetHaplotypes(bpGraph, useEpsilon);
					minColors = haplotypes.size();
					repsCompleted = repsCompleted + 1;
					updateTogetherCounts(haplotypes);
					int newNumCCs = findTogethersByCountThreshold(ccSeqList, repsCompleted).size();
					if(numCCs != newNumCCs) {
						repsUnchanged = 0;
						numCCs = newNumCCs;
					}
					else {
						repsUnchanged = repsUnchanged + 1;
					}
					bpGraph.shuffleOrder();
				}
			}
			else {
				for (int rep = 0; rep < numReps; rep++) {
					ArrayList<ArrayList<Sequence>> haplotypes = bipartiteMatchAndGetHaplotypes(bpGraph, useEpsilon);
					minColors = haplotypes.size();
					updateTogetherCounts(haplotypes);
					repsCompleted = repsCompleted + 1;
					if(repsCompleted >= maxReps) break;
					bpGraph.shuffleOrder();
				}
			}


			HashMap<DisjointSetNode, Set<Sequence>> rootsToNodeLists = findTogethersByCountThreshold(ccSeqList, repsCompleted);

			HaplotypeBlock newBlock = createNewHapBlockFromRootsToNodeLists(rootsToNodeLists, hapBlockNum.toString());
			newBlock.setRepsToFinish(repsCompleted);
			newBlock.setMinColors(minColors);
			alignment.addHaploBlock(newBlock);
			hapBlockNum = hapBlockNum + 1;
		}
		System.err.println("#DOT }");

	}





	/**
	 * Given a list of sequences that are "universal" (don't conflict with anything),
	 * returns a placeholder haplotype in a placeholder block (both named "U") for them to
	 * go in.
	 * @param universalSeqs
	 * @return
	 */

	protected Haplotype createUniversalHaplotype(HashSet<Sequence> universalSeqs) {
		// Make a block and haplotype for sequences that don't conflict with anything.
		Haplotype universalHaplotype = new Haplotype(alignment);

		universalHaplotype.setName("U");
		for(Sequence seq : universalSeqs) {
		//for(int i = 0; i < universalSeqs.size(); i++){
		//	universalHaplotype.addSequence(universalSeqs.get(i));
			universalHaplotype.addSequence(seq);
		}
		return universalHaplotype;
	}

	protected HaplotypeBlock createUniversalBlock(Haplotype universalHaplotype) {
		HaplotypeBlock universalBlock = new HaplotypeBlock(alignment);
		universalBlock.addHaplotype(universalHaplotype);
		universalHaplotype.setHaploBlock(universalBlock);
		universalBlock.setName("U");
		return universalBlock;
	}




	/**
	 * Given a graph with nodes containing sequences (edges containing counts of times those sequences appeared together),
	 * a list of those sequences, and a threshold, returns a hash structure where those sequences which all occured together more
	 * than threshold time are contained together in a list.
	 * O(n^2)
	 * @param togetherCounts
	 * @param ccSeqList
	 * @param threshold
	 * @return
	 */
	protected HashMap<DisjointSetNode, Set<Sequence>> findTogethersByCountThreshold(ArrayList<Sequence> ccSeqList, int threshold) {
		// Now we need to group up sequences in the togetherCounts graph that were connected numReps times
		// Making them haplotypes in the new block.

		// First we do union find to get the connected componenets of sequences that are always together
		ArrayList<DisjointSetNode> djNodes = new ArrayList<DisjointSetNode>();
		for (int j = 0; j < ccSeqList.size(); j++) {
			Sequence seqj = ccSeqList.get(j);
			djNodes.add(new DisjointSetNode(seqj));
		}
		for (int j = 0; j < djNodes.size() - 1; j++) {
			Sequence seqj = (Sequence)djNodes.get(j).getData();
			for (int k = j + 1; k < djNodes.size(); k++) {
				Sequence seqk = (Sequence)djNodes.get(k).getData();
				// TODO: memory hotspot, getEdge invokes a lot of iterators
				int edgeWeight = (int) togetherCounts.getEdgeWeight(togetherCounts.getEdge(seqk, seqj));
				//System.out.println("Comparing " + edgeWeight + " and " + numReps);
				// TODO: memory hotspot: containsEdge creates a lot of iterators
				if (togetherCounts.containsEdge(seqj, seqk) && edgeWeight >= threshold) {
					//System.out.println("Unioning!");
					djNodes.get(j).union(djNodes.get(k));
				}
			}
		}

		// Then we create lists of them that go together
		HashMap<DisjointSetNode, Set<Sequence>> rootsToNodeLists = new HashMap<DisjointSetNode, Set<Sequence>>();
		for (int j = 0; j < djNodes.size(); j++) {
			DisjointSetNode djNodej = djNodes.get(j);
			if (!rootsToNodeLists.containsKey(djNodej.find())) {
				Sequence thisSeq = (Sequence) djNodej.getData();
				Set<Sequence> newList = new HashSet<Sequence>();
				newList.add(thisSeq);
				rootsToNodeLists.put(djNodej.find(), newList);
			} else {
				Sequence thisSeq = (Sequence) djNodej.getData();
				Set<Sequence> theList = rootsToNodeLists.get(djNodej.find());
				theList.add(thisSeq);
				rootsToNodeLists.put(djNodej.find(), theList); // This probably isn't necessary, everything being pointers and all
			}
		}
		return rootsToNodeLists;
	}


	/**
	 * Given a hash from DisjointSetNode to sets of sequence, each set of sequences representing a haplotype, and all lists
	 * representing a haplotype block, builds this block and returns it.
	 * @param rootsToNodeLists
	 * @return
	 */
	protected HaplotypeBlock createNewHapBlockFromRootsToNodeLists(HashMap<DisjointSetNode, Set<Sequence>> rootsToNodeLists, String blockName) {
		// Create haplotypes from the groupings
		HaplotypeBlock newBlock = new HaplotypeBlock(alignment);

		Integer hapNum = 0;
		for (DisjointSetNode root : rootsToNodeLists.keySet()) {
			ArrayList<Sequence> hapList = extractSequencesFromSequenceSet(rootsToNodeLists.get(root));
			hapNum = hapNum + 1;
			Haplotype newHap = new Haplotype(alignment);
			newHap.setName(hapNum.toString());
			for (int k = 0; k < hapList.size(); k++) {
				Sequence seq = hapList.get(k);
				//System.out.println("Seq " + seq.getName() + " masked is " + seq.isMasked() );
				//why am I adding some masked sequences hur? I'm not: but I AM adding them to the universal haplotype
				newHap.addSequence(seq);
			}
			newBlock.addHaplotype(newHap);
			newHap.setHaploBlock(newBlock);
		}
		newBlock.setName(blockName);
		return newBlock;
	}

	/**
	 * Updates the togetherCounts graph structure based on list of haplotypes (haplotype=list of sequences)
	 * @param haplotypes
	 * @return
	 */
	protected void updateTogetherCounts(ArrayList<ArrayList<Sequence>> haplotypes) {
		for (int j = 0; j < haplotypes.size(); j++) {
			ArrayList<Sequence> haplotypeSeqs = haplotypes.get(j);

			// Fore each sequence pair in the haplotype created, increase their count int he togetherCounts graph
			for (int k = 0; k < haplotypeSeqs.size() - 1; k++) {
				Sequence seqk = haplotypeSeqs.get(k);
				for (int l = k + 1; l < haplotypeSeqs.size(); l++) {
					Sequence seql = haplotypeSeqs.get(l);
					// TODO: containsEdge is a memory hotspot, creates a lot of iterators. also getEdge
					if (togetherCounts.containsEdge(seqk, seql)) {
						double oldWeight = togetherCounts.getEdgeWeight(togetherCounts.getEdge(seqk, seql));
						double newWeight = oldWeight + 1;
						togetherCounts.setEdgeWeight(togetherCounts.getEdge(seqk, seql), newWeight);
					} else {
						//System.out.println("Trying to add and edge from " + seqk.getName() + " to " + seql.getName());
						togetherCounts.addEdge(seqk, seql);
						togetherCounts.setEdgeWeight(togetherCounts.getEdge(seqk, seql), 1.0);
					}
				}
			}

		}
	}

	protected ArrayList<ArrayList<Sequence>> bipartiteMatchAndGetHaplotypes(WeightedBipartiteGraph bpGraph, boolean useEpsilon) {


		// Match it, getting back our container graph
		//WeightedBipartiteGraph matchedBipartiteGraph = new HungarianAlgorithm().maxWeightedBipartiteMatching(bpGraph);
		WeightedBipartiteGraph matchedBipartiteGraph = bpGraph.maxWeightedBipartiteMatching();

		// After we've max matched it, we need to add edges from sequences on the left to their matches on the right
		SimpleWeightedGraph<GraphNode, DefaultWeightedEdge> matchedGraph = matchedBipartiteGraph.getGraph();
		ArrayList<GraphNode> leftNodes = matchedBipartiteGraph.getLeftNodes();
		ArrayList<GraphNode> rightNodes = matchedBipartiteGraph.getRightNodes();
		for (int j = 0; j < leftNodes.size(); j++) {
			matchedGraph.addEdge(leftNodes.get(j), rightNodes.get(j));
		}

		// After that, our connected components represent haplotypes.
		List<Set<GraphNode>> haplotypeCCs = new ConnectivityInspector(matchedGraph).connectedSets();


		ArrayList<ArrayList<Sequence>> haplotypes = new ArrayList<ArrayList<Sequence>>();
		for (int j = 0; j < haplotypeCCs.size(); j++) {
			Set<GraphNode> haplotypeCC = haplotypeCCs.get(j);
			ArrayList<Sequence> haplotypeSeqs = extractSequencesFromGraphNodeSet(haplotypeCC);
			haplotypes.add(haplotypeSeqs);
		}
		return haplotypes;
	}

	public String printBlockLex(HaplotypeBlock theBlock) {
		StringBuffer sb = new StringBuffer();

		ArrayList<Double> lexSize = theBlock.lexicographicSize();
		for (Double value : lexSize) {
			sb.append(value + " : ");
		}

		return sb.toString();
	}

	public ArrayList<HaplotypeBlock> getHaplotypeBlocks() {
		return alignment.getHaploBlocks();
	}

	/**
	 * Returns the summary strings returned by each HaplotypeBlock prepended with
	 * an index for which haplotype block they are in.
	 * @return
	 */
	public ArrayList<StringBuilder> summaryStrings() {
		ArrayList<StringBuilder> summaryStrings = new ArrayList<StringBuilder>();

		//for(HaplotypeBlock haploBlock: haplotypeBlocks) {
		ArrayList<HaplotypeBlock> haplotypeBlocks = alignment.getHaploBlocks();
		for (int j = 0; j < haplotypeBlocks.size(); j++) {
			HaplotypeBlock haploBlock = haplotypeBlocks.get(j);
			ArrayList<StringBuilder> haps = haploBlock.summaryStrings();
			//for(StringBuilder hap: haps) {
			for (int k = 0; k < haps.size(); k++) {
				StringBuilder hap = haps.get(k);
				StringBuilder thisBuilder = new StringBuilder();
				thisBuilder.append(haploBlock.getName());
				thisBuilder.append("\t");
				thisBuilder.append(hap.toString());
				summaryStrings.add(thisBuilder);
			}
		}

		return summaryStrings;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		ArrayList<HaplotypeBlock> haplotypeBlocks = alignment.getHaploBlocks();
		for(HaplotypeBlock block: haplotypeBlocks) {
			sb.append(block.toString());
			sb.append(System.getProperty("line.separator"));
		}
		return sb.toString();
	}

	public String haplotypeBlocksToString() {
		ArrayList<HaplotypeBlock> haplotypeBlocks = alignment.getHaploBlocks();

		StringBuffer sb = new StringBuffer();
		String blockSep = stringOfChars(alignment.getLength() + alignment.getMaxSeqNameLength() + 5, '#');

		//for(HaplotypeBlock block: haplotypeBlocks) {
		for (int i = 0; i < haplotypeBlocks.size(); i++) {
			HaplotypeBlock block = haplotypeBlocks.get(i);
			sb.append(blockSep);
			sb.append(System.getProperty("line.separator"));
			sb.append(block.toString());
			sb.append(System.getProperty("line.separator"));
		}

		return sb.toString();
	}

	protected String stringOfChars(int length, char theChar) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < length; i++) {
			sb.append(theChar);
		}
		return sb.toString();
	}

	/**
	 * Extracts unique sequences from a set of sequences (so if SeqA and SeqB both point to SEQ1, SEQ1 is
	 * only added to the returned list once)
	 * @param nodeSet
	 * @return
	 */
	protected ArrayList<Sequence> extractSequencesFromSequenceSet(Set<Sequence> nodeSet) {
		HashSet<Sequence> seqSet = new HashSet<Sequence>();
		ArrayList<Sequence> seqList = new ArrayList<Sequence>();

		for (Sequence theSeq : nodeSet) {
			//Sequence theSeq = (Sequence) node.getData();
			if (!seqSet.contains(theSeq)) {
				seqSet.add(theSeq);
			}
		}

		for (Sequence seq : seqSet) {
			seqList.add(seq);
		}

		return seqList;
	}

	/**
	 * Extracts unique sequences from a set of sequences (so if SeqA and SeqB both point to SEQ1, SEQ1 is
	 * only added to the returned list once)
	 * @param nodeSet
	 * @return
	 */
	protected ArrayList<Sequence> extractSequencesFromGraphNodeSet(Set<GraphNode> nodeSet) {
		HashSet<Sequence> seqSet = new HashSet<Sequence>();
		ArrayList<Sequence> seqList = new ArrayList<Sequence>();

		for (GraphNode node : nodeSet) {
			Sequence theSeq = (Sequence) node.getData();
			if (!seqSet.contains(theSeq)) {
				seqSet.add(theSeq);
			}
		}

		for (Sequence seq : seqSet) {
			seqList.add(seq);
		}

		return seqList;
	}
}




