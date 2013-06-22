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

/**
 * TODO: all the methods in this class are pretty stateful. I should use fewer class variables...
 * or at least make most of the methods not depend on the state so much
 * @author soneil
 */
public class Reconstructor {

		
	//private ArrayList<HaplotypeBlock> haploBlocks;
	//private ArrayList<SNP> snpList;


	// I should probably just make this static final...
	public Reconstructor() {

		//haploBlocks = new ArrayList<HaplotypeBlock>();
		//snpList = new ArrayList<SNP>();
	}




	/**
	 * Given a set of haplotypes, reconstructs an "optimal" path through them at SNP positions to create
	 * a consensusHumanReadable. Returns an object counting crossover points, haps used, etc.
	 * Optimizes on:
	 *  1) number of crossover points
	 *  2) number of crossovers into low-coverage haps
	 *  3) number of unique haplotypes used
	 *  4) total coverage of SNP alleles chosen
	 * @param haploBlocks
	 * @param alignment
	 * @return
	 */
	public ReconstructedConsensus computeReconstructedConsensus(MultipleAlignment alignment) {
		HashSet<SNP> theSnpList = alignment.getSNPs();

		if(theSnpList.size() == 0) {
			return noSNPsReconstruction(alignment);
		}
		else {
			ArrayList<HaplotypeBlock> haploBlocks = alignment.getHaploBlocks();

			ReconstructedConsensus recon = new ReconstructedConsensus();

			HashMap<SNP, ArrayList<SnpDAGNode>> snpToSnpDAGNodes = setupNodesAllHaps(haploBlocks, theSnpList);
			computeBackPointers(snpToSnpDAGNodes, theSnpList);
			ArrayList<SnpDAGNode> minSnpCompatibilityList = this.minSNPCompatibilityList(snpToSnpDAGNodes, theSnpList);
			String consensus = computeConsensus(minSnpCompatibilityList, alignment);
			ArrayList<SnpDAGNode> crossOverPoints = this.computeCrossOverPoints(minSnpCompatibilityList);

			recon.setConsensus(consensus);
			recon.setCrossOverPoints(crossOverPoints);
			recon.setMinSnpCompatibilityList(minSnpCompatibilityList);

			HashSet<Haplotype> uniqueHapsUsed = this.computeUniqueHapsUsed(minSnpCompatibilityList, crossOverPoints);
			//HashSet<Haplotype> uniqueLowCoverageHapsUsed = this.computeUniqueLowCoverageHapsUsed(minSnpCompatibilityList, crossOverPoints);

			recon.setUniqueHapsUsed(uniqueHapsUsed);
			recon.setSNPAlleleCoverage(computeSupportedCoverage(minSnpCompatibilityList));
			//recon.setUniqueLowCoverageHapsUsed(uniqueLowCoverageHapsUsed);

			return recon;
		}
	}



	/**
	 * Given a consensusHumanReadable sequence, computes an "optimal" path through the haplotypes that reconstructs it at the
	 * given snp positions (see computeReconstructedConsensus). 
	 * If no haplotype agrees at a non-gap snp position, this incurs a new crossover point into a low-coverage (1.0) haplotype.
	 * If punishGapsInConsensus is true, then this punishment occurs even if the consensusHumanReadable position is unknown (a ~ character)
	 * @param haploBlocks
	 * @param alignment
	 * @param consensusHumanReadable
	 * @return
	 * @throws Exception
	 */
	public ReconstructedConsensus evaluateConsensus(MultipleAlignment alignment, Sequence consensus, boolean punishGapsInConsensus) throws Exception {
		HashSet<SNP> theSnpList = alignment.getSNPs();

		if(theSnpList.size() == 0) {
			return noSNPsReconstruction(alignment);
		}
		else {
			ArrayList<HaplotypeBlock> haploBlocks = alignment.getHaploBlocks();

			ReconstructedConsensus recon = new ReconstructedConsensus();
			HashSet<SNP> coveredSnps = new HashSet<SNP>();

			for(SNP snpi : theSnpList) {
			//for(int i = 0; i < theSnpList.size(); i++) {
				//SNP snpi = theSnpList.get(i);
				if(consensus.containsPosition(snpi.getPosition()) || punishGapsInConsensus) {
					coveredSnps.add(snpi);
				}
			}

			HashMap<SNP, ArrayList<SnpDAGNode>> snpToSnpDAGNodes = setupNodesHapsComptabileWithSeq(haploBlocks, coveredSnps, consensus);
			computeBackPointers(snpToSnpDAGNodes, coveredSnps);
			ArrayList<SnpDAGNode> minSnpCompatibilityList = this.minSNPCompatibilityList(snpToSnpDAGNodes, coveredSnps);
			ArrayList<SnpDAGNode> crossOverPoints = this.computeCrossOverPoints(minSnpCompatibilityList);

			recon.setConsensus(consensus.sequenceToString(true));
			recon.setCrossOverPoints(crossOverPoints);
			recon.setMinSnpCompatibilityList(minSnpCompatibilityList);

			HashSet<Haplotype> uniqueHapsUsed = this.computeUniqueHapsUsed(minSnpCompatibilityList, crossOverPoints);
			//HashSet<Haplotype> uniqueLowCoverageHapsUsed = this.computeUniqueLowCoverageHapsUsed(minSnpCompatibilityList, crossOverPoints);

			recon.setUniqueHapsUsed(uniqueHapsUsed);
			recon.setSNPAlleleCoverage(computeSupportedCoverage(minSnpCompatibilityList));
			//recon.setUniqueLowCoverageHapsUsed(uniqueLowCoverageHapsUsed);

			return recon;
		}
	}

	/**
	 * wraper for using evaluateConsensus but given a string consensus rather than a sequence consensus
	 * @param alignment
	 * @param consensus
	 * @param punishGapsInConsensus
	 * @return
	 * @throws Exception
	 */
	public ReconstructedConsensus evaluateConsensus(MultipleAlignment alignment, String consensus, boolean punishGapsInConsensus) throws Exception {
			ArrayList<String> stringArray = new ArrayList<String>();
			stringArray.add(consensus);
			Sequence newSeq = new Sequence();
			newSeq.addPieces(stringArray, 0);

			return this.evaluateConsensus(alignment, newSeq, punishGapsInConsensus);
	}

	/**
	 * Given a list of SnpDAGNodes which describes the haplotypes to be used for each snp in
	 * the reconstruction, compute the total haplotype supporting coverage at those snp loci
	 * @param minSnpCompatibilityList
	 * @return
	 */
	private int computeSupportedCoverage(ArrayList<SnpDAGNode> minSnpCompatibilityList) {
		int sum = 0;

		for(int i = 0; i < minSnpCompatibilityList.size(); i++) {
			SnpDAGNode nodei = minSnpCompatibilityList.get(i);
			sum = sum + nodei.getCompatibleHap().coverageAtPosition(nodei.getSnp().getPosition());
		}

		return sum;
	}


	/**
	 * If there are no SNPs in the alignment, return a reconstructed consensus which is just the majority
	 * vote of the alignment and specifies the universal haplotype as the only one to use
	 * @param alignment
	 * @return
	 */
	private ReconstructedConsensus noSNPsReconstruction(MultipleAlignment alignment) {
		ReconstructedConsensus recon = new ReconstructedConsensus();
		recon.setConsensus(alignment.majorityVoteConsensus());
		recon.setSNPAlleleCoverage(0);

		ArrayList<SnpDAGNode> universalDAGNodeList = new ArrayList<SnpDAGNode>();
		universalDAGNodeList.add(new SnpDAGNode(new SNP(0, alignment),alignment.getUniversalHaplotype()));

		recon.setCrossOverPoints(universalDAGNodeList);
		recon.setMinSnpCompatibilityList(universalDAGNodeList);

		HashSet<Haplotype> hapsUsed = new HashSet<Haplotype>();
		hapsUsed.add(alignment.getUniversalHaplotype());
		recon.setUniqueHapsUsed(hapsUsed);
		//recon.setUniqueLowCoverageHapsUsed(uniqueLowCoverageHapsUsed);
		return recon;
	}




		/**
	 * Given a list of haploBlocks, extracts the haplotypes from them, and for each
	 * snp creates a series of SnpDAGNodes which are put into the main data structure, but only using
	 * haplotypes that are consistent with the given sequence at the SNP. If the given sequence at a SNP position
	 * '~,' then all haps are compatible with it.
	 * This way we can evaluate a consensusHumanReadable from a set of reads which does't completely cover the ground truth
	 * (for example, near the ends where coverage is 0), even if the ground truth has snps there.
	 * IF NO HAP IS CONSISTENT AT THE GIVEN SNP THEN A "NEW" HAP IS CREATED JUST FOR THAT SNP
	 * WHICH MUST BE CROSSED OVER INTO AND BACK OUT OF. This also happens if no haplotype _covers_ a given snp
	 * (which can happen if we are evaluating a consensusHumanReadable with "missing" information relative to the haploblocks.
	 * @param haploBlocks
	 * @param snpList
	 */
	private HashMap<SNP, ArrayList<SnpDAGNode>> setupNodesHapsComptabileWithSeq(ArrayList<HaplotypeBlock> haploBlocks, HashSet<SNP> theSnpList, Sequence theSeq) throws Exception {
		// get snps and haplotypes
		ArrayList<Haplotype> haplotypes = new ArrayList<Haplotype>();
		HashMap<SNP, ArrayList<SnpDAGNode>> snpToSnpDAGNodes = new HashMap<SNP, ArrayList<SnpDAGNode>>();


		HashSet<SNP> snpList = theSnpList;
		for(int i = 0; i < haploBlocks.size(); i++) {
			HaplotypeBlock blocki = haploBlocks.get(i);
			ArrayList<Haplotype> blockiHaps = blocki.getHaplotypes();
			for(int j = 0; j < blockiHaps.size(); j++) {
				haplotypes.add(blockiHaps.get(j));
			}
		}

		for(SNP snpi : snpList) {
		//for(int i = 0; i < snpList.size(); i++) {
			//SNP snpi = snpList.get(i);
			int snpiPos = snpi.getPosition();
			if(!snpToSnpDAGNodes.containsKey(snpi)) snpToSnpDAGNodes.put(snpi, new ArrayList<SnpDAGNode>());
			for(int j = 0; j < haplotypes.size(); j++) {
				Haplotype hapj = haplotypes.get(j);
				char hapjCons = hapj.consensusOverallAtPosition(snpiPos);
				char seqCons = theSeq.baseAtAlignmentPosition(snpiPos);
				if(hapj.coversPosition(snpiPos)) {
					/*if(seqCons == '~') {
						throw new Exception("Oops, you want to check to see whether hap " + hapj.fullName() + " is compatible with sequence " + theSeq.getName() + " at position " + snpi.getPosition() + ", but there is a ~ character there.");
					}*/
					if(hapjCons == seqCons || seqCons == '~') { // In this version, we only add haps compatible with the seq, and of course covering the snp position
						SnpDAGNode newNode = new SnpDAGNode(snpi, hapj);
						snpToSnpDAGNodes.get(snpi).add(newNode);
					}
				}
			}
			if(snpToSnpDAGNodes.get(snpi).size() == 0) { // if that's the case, then we create a brand new haplotype of weight 1.0 (containing just the seq and a brand new name), so that the weight will be shitty. hehehe
				Haplotype newHap = new Haplotype(haplotypes.get(0).getHaploBlock().getAlignment());
				newHap.setHaploBlock(haplotypes.get(0).getHaploBlock());
				newHap.setName("UnsupportedAt" + theSeq.getName() + snpi.getPosition());
				newHap.addSequence(theSeq);
				SnpDAGNode newNode = new SnpDAGNode(snpi,newHap);
				snpToSnpDAGNodes.get(snpi).add(newNode);
			}
		}

		return snpToSnpDAGNodes;
	}




	/**
	 * Given a list of haploBlocks, extracts the haplotypes from them, and for each
	 * snp creates a series of SnpDAGNodes which are put into the main data structure
	 * @param haploBlocks
	 * @param snpList
	 */
	private HashMap<SNP, ArrayList<SnpDAGNode>> setupNodesAllHaps(ArrayList<HaplotypeBlock> haploBlocks, HashSet<SNP> theSnpList) {
		// get snps and haplotypes
		ArrayList<Haplotype> haplotypes = new ArrayList<Haplotype>();
		HashMap<SNP, ArrayList<SnpDAGNode>> snpToSnpDAGNodes = new HashMap<SNP, ArrayList<SnpDAGNode>>();


		HashSet<SNP> snpList = theSnpList;
		for(int i = 0; i < haploBlocks.size(); i++) {
			HaplotypeBlock blocki = haploBlocks.get(i);
			ArrayList<Haplotype> blockiHaps = blocki.getHaplotypes();
			for(int j = 0; j < blockiHaps.size(); j++) {
				haplotypes.add(blockiHaps.get(j));
			}
		}

		for(SNP snpi : snpList) {
		//for(int i = 0; i < snpList.size(); i++) {
			//SNP snpi = snpList.get(i);
			int snpiPos = snpi.getPosition();
			if(!snpToSnpDAGNodes.containsKey(snpi)) snpToSnpDAGNodes.put(snpi, new ArrayList<SnpDAGNode>());
			for(int j = 0; j < haplotypes.size(); j++) {
				Haplotype hapj = haplotypes.get(j);
				if(hapj.coversPosition(snpiPos)) {
					SnpDAGNode newNode = new SnpDAGNode(snpi, hapj);
					snpToSnpDAGNodes.get(snpi).add(newNode);
				}
			}
		}

		return snpToSnpDAGNodes;
	}

	/**
	 * Ok, now we can do the dynamic programming to find the min-chimeric path through the snps
	 */
	private void computeBackPointers(HashMap<SNP, ArrayList<SnpDAGNode>> snpToSnpDAGNodes, HashSet<SNP> snpList) {
		ArrayList<SNP> snpListArray = new ArrayList<SNP>();
		for(SNP snpi : snpList) {
			snpListArray.add(snpi);
		}
		Collections.sort(snpListArray, new SNPPosComparator());
		// base of the dynamic program
		ArrayList<SnpDAGNode> firstNodes = snpToSnpDAGNodes.get(snpListArray.get(0));
		for(int i = 0; i < firstNodes.size(); i++) {
			SnpDAGNode nodei = firstNodes.get(i);
			nodei.setMinChimerisms(0);
			if(nodei.getCompatibleHap().averageCoverage() < 2) nodei.setMinLowCount(1);
			else nodei.setMinLowCount(0);
		}
		// we have dependecies between adjacent snps, we update information in a feed-forward fashion
		for(int i = 0; i < snpListArray.size()-1; i++) {
			SNP thisSnp = snpListArray.get(i);
			SNP nextSnp = snpListArray.get(i+1);
			ArrayList<SnpDAGNode> thisNodeList = snpToSnpDAGNodes.get(thisSnp);
			ArrayList<SnpDAGNode> nextNodeList = snpToSnpDAGNodes.get(nextSnp);
			for(int j = 0; j < thisNodeList.size(); j++) {
				SnpDAGNode thisNodej = thisNodeList.get(j);
				for(int k = 0; k < nextNodeList.size(); k++) {
					SnpDAGNode nextNodek = nextNodeList.get(k);
					if(thisNodej.getCompatibleHap() == nextNodek.getCompatibleHap()) {
						// Here's the candidate information, which we'll compare to
						int candidateChimerisms = thisNodej.getMinChimerisms() + 0;
						int candidateLowCount = thisNodej.getMinLowCount() + 0; // if we're staying the same hap, we won't add more lowCount even if it is low since we already added it
						int candidateMinUsageSize = thisNodej.getMinUsage().size();
							if(!thisNodej.getMinUsage().contains(nextNodek.getCompatibleHap())) candidateMinUsageSize = candidateMinUsageSize + 1; //I don't think this will ever fire since they're equal...
						int candidateMaxCoverage = thisNodej.getMaxCoverage() + nextNodek.getCompatibleHap().coverageAtPosition(nextSnp.getPosition());
						// Now we compare it...
						if(isCandidateBetterThanBest(candidateChimerisms, nextNodek.getMinChimerisms(), candidateLowCount, nextNodek.getMinLowCount(), candidateMinUsageSize, nextNodek.getMinUsage().size(), candidateMaxCoverage, nextNodek.getMaxCoverage()) >= 0) {
							nextNodek.setBackPointer(thisNodej);
							nextNodek.setMinChimerisms(candidateChimerisms);
							nextNodek.setMinLowCount(candidateLowCount);
							nextNodek.setMinUsage(thisNodej.getMinUsage()); //we don't need to add the next hap to the usage, since it's already in the present's list, since they are equal
							nextNodek.setMaxCoverage(candidateMaxCoverage);
						}
					}
					else {
						// Here's the candidate information, which we'll compare to
						int candidateChimerisms = thisNodej.getMinChimerisms() + 1;
						int candidateLowCount = thisNodej.getMinLowCount();
							if(nextNodek.getCompatibleHap().averageCoverage() < 2) candidateLowCount = candidateLowCount + 1;
						int candidateMinUsageSize = thisNodej.getMinUsage().size();
							if(!thisNodej.getMinUsage().contains(nextNodek.getCompatibleHap())) candidateMinUsageSize = candidateMinUsageSize + 1;
						// Now we compare it...
						int candidateMaxCoverage = thisNodej.getMaxCoverage() + nextNodek.getCompatibleHap().coverageAtPosition(nextSnp.getPosition());
						if(isCandidateBetterThanBest(candidateChimerisms, nextNodek.getMinChimerisms(), candidateLowCount, nextNodek.getMinLowCount(), candidateMinUsageSize, nextNodek.getMinUsage().size(), candidateMaxCoverage, nextNodek.getMaxCoverage()) >= 0) {
							nextNodek.setBackPointer(thisNodej);
							nextNodek.setMinChimerisms(candidateChimerisms);
							nextNodek.setMinLowCount(candidateLowCount);
							Haplotype candidateNewHap = nextNodek.getCompatibleHap();
							HashSet<Haplotype> oldUsage = thisNodej.getMinUsage();
							oldUsage.add(candidateNewHap);
							nextNodek.setMinUsage(oldUsage);
							nextNodek.setMaxCoverage(candidateMaxCoverage);
						}
					}
				}
			}
		}
	}

	/**
	 * Traces the list of backpointers back to get the consensusHumanReadable
	 * @return
	 */
	private ArrayList<SnpDAGNode> minSNPCompatibilityList(HashMap<SNP, ArrayList<SnpDAGNode>> snpToSnpDAGNodes, HashSet<SNP> snpList) {
		//SNP lastSnp = snpList.get(snpList.size()-1);
		SNP lastSnp = null;
		int lastSNPposition = -1;
		for(SNP snpi : snpList) {
			if(snpi.getPosition() > lastSNPposition) {
				lastSnp = snpi;
				lastSNPposition = snpi.getPosition();
			}
		}
		ArrayList<SnpDAGNode> lastNodes = snpToSnpDAGNodes.get(lastSnp);
		SnpDAGNode bestNode = lastNodes.get(0);
		int bestChimerisms = bestNode.getMinChimerisms();
		int bestLowCount = bestNode.getMinLowCount();
		int bestMinUsageSize = bestNode.getMinUsage().size();
		int bestMaxCoverage = bestNode.getMaxCoverage();
		for(int i = 0; i < lastNodes.size(); i++) {
			SnpDAGNode candidateNode = lastNodes.get(i);
			int candidateChimerisms = candidateNode.getMinChimerisms();
			int candidateLowCount = candidateNode.getMinLowCount();
			int candidateMinUsageSize = candidateNode.getMinUsage().size();
			int candidateMaxCoverage = candidateNode.getMaxCoverage();
			if(isCandidateBetterThanBest(candidateChimerisms, bestChimerisms, candidateLowCount, bestLowCount, candidateMinUsageSize, bestMinUsageSize, candidateMaxCoverage, bestMaxCoverage) >= 0) {
				bestNode = candidateNode;
				bestChimerisms = candidateChimerisms;
				bestLowCount = candidateLowCount;
				bestMinUsageSize = candidateMinUsageSize;
				bestMaxCoverage = candidateMaxCoverage;
			}
		}
		return bestNode.backPointersList();
	}



	private void printDebugConsensus(ArrayList<SnpDAGNode> nodeList) {
		for(int i = 0; i < nodeList.size(); i++) {
			System.out.println("Snp position: " + nodeList.get(i).getSnp().getPosition() + " hap: " + nodeList.get(i).getCompatibleHap().getName());
		}
	}


	/**
	 * Returns the consensusHumanReadable; which is the consensusHumanReadable of the haplotypes determined by the
	 * minSNPCompatibilityList with the consensusHumanReadable of the overall alignment filling in in between
	 * @return
	 */
	private String computeConsensus(ArrayList<SnpDAGNode> minSnpCompatibilityList, MultipleAlignment alignment) {
		StringBuilder sb = new StringBuilder();
		int alignmentLength = alignment.getLength();

		for(int i = 0; i < alignmentLength; i++) {
			sb.append(alignment.majorityVoteConsensusAtPosition(i));
		}

		Collections.sort(minSnpCompatibilityList, new SnpDAGNodePosComparator());

		for(int i = 0; i < minSnpCompatibilityList.size(); i++) {
			SnpDAGNode nodei = minSnpCompatibilityList.get(i);
			int posi = nodei.getSnp().getPosition();
			char consChar = nodei.getCompatibleHap().consensusOfHapAtPosition(posi);
			sb.setCharAt(posi, consChar);
		}

		return sb.toString();
	}

	/*
	 * Computes the crossOverPoints (snp positions in the minSNPCompatibilityList that indicate a
	 * change in haplotype, including the first haplotype)
	 */
	private ArrayList<SnpDAGNode> computeCrossOverPoints(ArrayList<SnpDAGNode> minSnpCompatibilityList) {
		ArrayList<SnpDAGNode> coPoints = new ArrayList<SnpDAGNode>();
		Haplotype currentHap = minSnpCompatibilityList.get(0).getCompatibleHap();
		coPoints.add(minSnpCompatibilityList.get(0));
		for(int i = 1; i < minSnpCompatibilityList.size(); i++) {
			SnpDAGNode nodei = minSnpCompatibilityList.get(i);
			Haplotype nextHap = nodei.getCompatibleHap();
			if(nextHap != currentHap) coPoints.add(nodei);
			currentHap = nextHap;
		}
		return coPoints;
	}

	private HashSet<Haplotype> computeUniqueHapsUsed(ArrayList<SnpDAGNode> minSnpCompatibilityList, ArrayList<SnpDAGNode> crossOverPoints) {
		SnpDAGNode firstNode = minSnpCompatibilityList.get(0);
		Haplotype firstHap = firstNode.getCompatibleHap();
		HashSet<Haplotype> uniqueHapsUsed = new HashSet<Haplotype>();
		uniqueHapsUsed.add(firstHap);

		for(int i = 0; i < crossOverPoints.size(); i++) {
			SnpDAGNode nodei = crossOverPoints.get(i);
			Haplotype hapi = nodei.getCompatibleHap();
			if(!uniqueHapsUsed.contains(hapi)) {
				uniqueHapsUsed.add(hapi);
			}
		}
		return uniqueHapsUsed;
	}

	private HashSet<Haplotype> computeUniqueLowCoverageHapsUsed(ArrayList<SnpDAGNode> minSnpCompatibilityList, ArrayList<SnpDAGNode> crossOverPoints) {
		SnpDAGNode firstNode = minSnpCompatibilityList.get(0);
		HashSet<Haplotype> uniqueLowCoverageHapsUsed = new HashSet<Haplotype>();
		Haplotype firstHap = firstNode.getCompatibleHap();
		if(firstHap.averageCoverage() < 2) {
			uniqueLowCoverageHapsUsed.add(firstHap);
		}
		for(int i = 0; i < crossOverPoints.size(); i++) {
			SnpDAGNode nodei = crossOverPoints.get(i);
			Haplotype hapi = nodei.getCompatibleHap();

			if(hapi.averageCoverage() < 2 && ! uniqueLowCoverageHapsUsed.contains(hapi)) {
				uniqueLowCoverageHapsUsed.add(hapi);
			}
		}
		return uniqueLowCoverageHapsUsed;
	}


	/**
	 * REturns 1 if the cadidate set of params is better than the best so far, 0 if they are equal, -1 if they are worse
	 * @param candidateChimerisms
	 * @param currentMinChimerisms
	 * @param candidateLowCount
	 * @param currentMinLowCount
	 * @param candidateUsageSize
	 * @param currentMinUsageSize
	 * @param candidateMaxCoverage
	 * @param currentMaxCoverage
	 * @return
	 */
	public int isCandidateBetterThanBest(int candidateChimerisms, int currentMinChimerisms, int candidateLowCount, int currentMinLowCount, int candidateUsageSize, int currentMinUsageSize, int candidateMaxCoverage, int currentMaxCoverage) {
		if(candidateChimerisms < currentMinChimerisms) {
			return 1;
		}
		else if(candidateChimerisms == currentMinChimerisms && candidateMaxCoverage > currentMaxCoverage) {
			return 1;
		}
		else if(candidateChimerisms == currentMinChimerisms && candidateMaxCoverage == currentMaxCoverage) {
			return 0;
		}
		else {
			return -1;
		}
		/*else if(candidateChimerisms == currentMinChimerisms && candidateLowCount < currentMinLowCount) {
			return 1;
		}
		else if(candidateChimerisms == currentMinChimerisms && candidateLowCount == currentMinLowCount && candidateUsageSize < currentMinUsageSize) {
			return 1;
		}
		else if(candidateChimerisms == currentMinChimerisms && candidateLowCount == currentMinLowCount && candidateUsageSize == currentMinUsageSize && candidateMaxCoverage > currentMaxCoverage) {
			return 1;
		}
		else if(candidateChimerisms == currentMinChimerisms && candidateLowCount == currentMinLowCount && candidateUsageSize == currentMinUsageSize && candidateMaxCoverage == currentMaxCoverage) {
			return 0;
		}
		else {
			return -1;
		}*/
	}

	
}
