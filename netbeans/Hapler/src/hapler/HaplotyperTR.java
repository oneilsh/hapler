/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hapler;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Is inheriting from HaplotyperCompat a good way to do this? I dunno, see:
 * http://stackoverflow.com/questions/49002/prefer-composition-over-inheritance
 * In the future maybe I will define an interface. Blech.
 * @author soneil
 */
public class HaplotyperTR extends HaplotyperCompat {

	public HaplotyperTR(MultipleAlignment theAlignment, boolean theUseEpsilon, int theNumReps, int theMaxReps) throws Exception {
		super(theAlignment, theUseEpsilon, theNumReps, theMaxReps);
	}

	
	/**
	 * 
	 * @throws Exception
	 */
	@Override
	public void execute() throws Exception {
		// Seqs which contain a SNP _must_ conflict with something. The parser should catch this, if not, we'll throw an exception below.
		// Seqs which don't contain a SNP must be univeral, and will go into the universal haplotype.
		// These sequences are masked, but they aren't maked _by_ anything.
		//TODO: what about user-defined SNPs that don't actually cover any variant positions? The parser should probably catch this...

		// just fucking around...
		HashSet<Sequence> universalSeqs = alignment.maskUnmaskedNoSNPSequences();
		Haplotype universalHap = this.createUniversalHaplotype(universalSeqs);
		HaplotypeBlock universalBlock = this.createUniversalBlock(universalHap);
		alignment.setUniversalHaplotype(universalHap);
		alignment.addHaploBlock(universalBlock);

		// Mask the rest of the sequences that need to be masked
		System.err.println(alignment.getName() + ": Masking redundant sequences...");
		alignment.maskUnmaskedRedundantSeqsByUnmasked(); // just fucking around...

		System.err.println(alignment.getName() + ": Masking of redundant seqs completed. Computing haplotype blocks (connected components of the conflict graph)...");
		// just fucking around...
		//ArrayList<ArrayList<Sequence>> connectedComponents = new ArrayList<ArrayList<Sequence>>();
		//connectedComponents.add(alignment.getSequences());
		ArrayList<ArrayList<Sequence>> connectedComponents = alignment.unmaskedConflictGraphConnectedComponents();

		System.err.println(alignment.getName() + ": Blocks computed. Processing " + connectedComponents.size() + " haplotpe blocks...");

		Integer hapBlockNum = 1;
		Integer minColors = -1;

		System.err.println("#DOT strict digraph G{");
		System.err.println("#DOT nodesep = 0.02;");
		for(int i = 0; i < connectedComponents.size(); i++) {
			ArrayList<Sequence> ccSeqList = connectedComponents.get(i);
			System.err.println(alignment.getName() + ": Processing Haplotype Block " + hapBlockNum + ", number non-redundant sequences: " + ccSeqList.size());

			// Maintain a graph counting how many times each sequence was put together into the same haplotype
			for(int j = 0; j < ccSeqList.size(); j++) {
				Sequence seqj = ccSeqList.get(j);
				togetherCounts.addVertex(seqj);
			}


			// We do the bipartite matching a number of times, shuffling the input after each one
			// Create the bipartite representation
			WeightedBipartiteGraph bpGraph = new WeightedBipartiteGraph(alignment.positiveTRGraphFromSequences(ccSeqList), useEpsilon);

			int repsCompleted = 0; int repsUnchanged = 0; int numCCs = 0;
			if(numReps == -1) { // If the numReps is -1, we want "auto," ie., run until we get 10% repetitions with no change in the number of haplotypes, or 10 reps, whichever is larger
										// but don't do more than maxReps...
				while(repsCompleted <= maxReps && (repsCompleted < 20 || repsUnchanged < 0.5*repsCompleted)) {

					ArrayList<ArrayList<Sequence>> haplotypes = super.bipartiteMatchAndGetHaplotypes(bpGraph, useEpsilon);
					minColors = haplotypes.size();
					repsCompleted = repsCompleted + 1;
					super.updateTogetherCounts(haplotypes);
					int newNumCCs = super.findTogethersByCountThreshold(ccSeqList, repsCompleted).size();
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
					ArrayList<ArrayList<Sequence>> haplotypes = super.bipartiteMatchAndGetHaplotypes(bpGraph, useEpsilon);
					minColors = haplotypes.size();
					super.updateTogetherCounts(haplotypes);
					repsCompleted = repsCompleted + 1;
					if(repsCompleted >= maxReps) break;
					bpGraph.shuffleOrder();
				}
			}


			HashMap<DisjointSetNode, Set<Sequence>> rootsToNodeLists = super.findTogethersByCountThreshold(ccSeqList, repsCompleted);

			HaplotypeBlock newBlock = super.createNewHapBlockFromRootsToNodeLists(rootsToNodeLists, hapBlockNum.toString());
			newBlock.setRepsToFinish(repsCompleted);
			newBlock.setMinColors(minColors);
			alignment.addHaploBlock(newBlock);
			hapBlockNum = hapBlockNum + 1;

		}
		System.err.println("#DOT }");
		//System.exit(0);
	}
}
