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
 * A haplotype is a set of consistent sequences. No sequence may be present more than once.
 * @author soneil
 */
public class Haplotype {
	private HashSet<Sequence> haplotype;
	private MultipleAlignment alignment;
	private HaplotypeBlock haploBlock;
	private HashMap<Integer,SNP> snpsCovered;
	private String name;



	public Haplotype(MultipleAlignment theAlignment) {
		haplotype = new HashSet<Sequence>();
		alignment = theAlignment;
		snpsCovered = new HashMap<Integer,SNP>();
	}

	public void addSequence(Sequence theSequence) {
		if(!haplotype.contains(theSequence)) {
			haplotype.add(theSequence);
			HashSet<SNP> possibleNewSnps = theSequence.getSNPsCovered();
			for(SNP possibleNewSnp : possibleNewSnps) {
			//for(int i = 0; i < possibleNewSnps.size(); i++) {
			//	SNP possibleNewSnp = possibleNewSnps.get(i);
				if(!snpsCovered.containsKey(possibleNewSnp.getPosition())) {
					snpsCovered.put(possibleNewSnp.getPosition(), possibleNewSnp);
				}
			}
		}
	}



	/**
	 * Given a list of sequecces, returns a list of them that is compatible with this
	 * haplotype (ie, don't conflict at covered SNP positions in the given)
	 * @param hapList
	 * @return
	 */
	public HashMap<Sequence,Integer> mismatchCounts(ArrayList<Sequence> seqList) {
		HashMap<Sequence,Integer> seqToMisMatchCount = new HashMap<Sequence, Integer>();
		for(int i = 0; i < seqList.size(); i++) {
			Sequence seqi = seqList.get(i);
			seqToMisMatchCount.put(seqi, sequenceMismatchCount(seqi));
		}

		return seqToMisMatchCount;
	}




	/**
	 * returns true if the start position of this haplotype is BEFORE the start position of the other
	 * @param other
	 * @return
	 */
	public boolean startsBefore(Haplotype other) {
		if(this.startPos() < other.startPos()) return true;
		else return false;
	}

	/**
	 * returns true if the end position of this haplotype is BEFORE the end position of the other
	 * @param other
	 * @return
	 */
	public boolean endsBefore(Haplotype other) {
		if(this.endPos() < other.endPos()) return true;
		else return false;
	}

	/**
	 * returns true if this sequence overlaps the other, as determined by start and end positions
	 * of both
	 * @param other
	 * @return
	 */
	public boolean overlaps(Haplotype other) {
		boolean toRet = true;
		if(this.startPos() > other.endPos() || other.startPos() > this.endPos()) {
			toRet = false;
		}
		return toRet;
	}


	/**
	 * Returns a list of subhaplotypes: a subhaplotype is defined as a set of sequences that are part of the
	 * same haplotype AND are part of the same connected component in the overlap graph. Thus, if the sequences
	 * are gapless, each subhaplotype will be contiguous.
	 * @return
	 */
	/*public ArrayList<Haplotype> getSubHaplotypes() {
		ArrayList<Haplotype> subHaps = new ArrayList<Haplotype>();

		
		ArrayList<DisjointSetNode> djNodes = new ArrayList<DisjointSetNode>();
		for(int j = 0; j < haplotype.size(); j++) {
			Sequence seq = haplotype.get(j);
			if(!seq.isMasked()) {
				djNodes.add(new DisjointSetNode(seq));
			}
		}

		// union the connected components
		for(int i = 0; i < djNodes.size()-1; i++) {
			DisjointSetNode inode = djNodes.get(i);
			for(int j = i+1; j < djNodes.size(); j++) {
				DisjointSetNode jnode = djNodes.get(j);
				Sequence iseq = (Sequence)inode.getData();
				Sequence jseq = (Sequence)jnode.getData();
				if(iseq != jseq && iseq.overlaps(jseq)) {
					inode.union(jnode);
				}
			}
		}

		// each list represents a subhaplotype
		HashMap<DisjointSetNode, ArrayList<Sequence>> rootsToSeqLists = new HashMap<DisjointSetNode, ArrayList<Sequence>>();
		for(int i = 0; i < djNodes.size(); i++) {
			Sequence seqi = (Sequence)djNodes.get(i).getData();
			DisjointSetNode root = djNodes.get(i).find();
			if(!rootsToSeqLists.containsKey(root)) {
				rootsToSeqLists.put(root, new ArrayList<Sequence>());
			}
			rootsToSeqLists.get(root).add(seqi);
		}

		int subHapNumber = 1;
		for(DisjointSetNode djNode : rootsToSeqLists.keySet()) {
			ArrayList<Sequence> subHapSeqList = rootsToSeqLists.get(djNode);
			Haplotype subHap = new Haplotype(alignment);
			for(int i = 0; i < subHapSeqList.size(); i++) {
				subHap.addSequence(subHapSeqList.get(i));
			}
			subHap.setName(this.getName());
			subHaps.add(subHap);
			subHapNumber = subHapNumber+1;
		}

		return subHaps;
	}*/


	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(this.getName());
		sb.append(System.getProperty("line.separator"));


		String consensus = consensusHumanReadable();
		int consLength = consensus.length();
		int maxSeqNameLength = alignment.getMaxSeqNameLength();



		// Display all SNP positions covered.
		ArrayList<SNP> snps = this.consistentSNPsCovered();

		char[] snpCharArray = new char[consLength+maxSeqNameLength+5];
		for(int i = 0; i < snpCharArray.length; i++) {
			snpCharArray[i] = ' ';
		}
		//for(SNP snp : snps) {
		for(int i = 0; i < snps.size(); i++) {
			SNP snp = snps.get(i);
			snpCharArray[snp.getPosition()+maxSeqNameLength+5] = '*';
		}
		for(int i = 0; i < snpCharArray.length; i++) {
			char c = snpCharArray[i];
			sb.append(c);
		}
		sb.append(System.getProperty("line.separator"));



		// Display all internally variant positions
		snps = this.allVariantPositionsAsSNPs();
		snpCharArray = new char[consLength+maxSeqNameLength+5];
		for(int i = 0; i < snpCharArray.length; i++) {
			snpCharArray[i] = ' ';
		}
		//for(SNP snp : snps) {
		for(int i = 0; i < snps.size(); i++) {
			SNP snp = snps.get(i);
			snpCharArray[snp.getPosition()+maxSeqNameLength+5] = '+';
		}
		for(int i = 0; i < snpCharArray.length; i++) {
			char c = snpCharArray[i];
			sb.append(c);
		}
		sb.append(System.getProperty("line.separator"));



		sb.append(String.format("%1$" + (maxSeqNameLength +1) + "s    ", "CONS"));
		sb.append(consensus);
		sb.append(System.getProperty("line.separator"));


		ArrayList<Sequence> hapSeqs = new ArrayList<Sequence>(haplotype);
		Collections.sort(hapSeqs, new SequenceStartPosLengthComparator());
		for(Sequence seq : hapSeqs) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seq = haplotype.get(i);
			if(!seq.isMasked()) {
				sb.append(String.format("%1$" + (maxSeqNameLength +1) + "s    ", seq.getName()));
				sb.append(seq.sequenceToString(false) + System.getProperty("line.separator"));
				// print the seqs this one is masking below it, reversing the case?
				ArrayList<Sequence> masking = seq.getMasking();
				//for(Sequence masked : seq.getMasking()) {
				for(int j = 0; j < masking.size(); j++) {
					Sequence masked = masking.get(j);
					sb.append("+");
					sb.append(String.format("%1$" + (maxSeqNameLength) + "s    ", masked.getName()));
					sb.append(masked.sequenceToStringReverseCase(false) + System.getProperty("line.separator"));
				}
				//if(seq.getMasking().size() != 0) sb.append(System.getProperty("line.separator"));
			}
		}
		return sb.toString();
	}

	/**
	 * Returns the number of non-~ characters in this haplotype
	 * @return
	 */
	public int numBasesCovered() {
		int count = 0;
		for(int i = 0; i < alignment.getLength(); i++) {
			if(this.coversPosition(i)) {
				count = count + 1;
			}
		}
		return count;
	}


	/**
	 *
	 * @return
	 */
	public StringBuilder summaryString() {
		StringBuilder sb = new StringBuilder();
		sb.append(numNonRedundantSequences() + "\t");
		sb.append(numRedundantSequences() + "\t");
		sb.append(numSNPsCovered() + "\t");
		//sb.append(numInconsistentSNPs() + "\t");
		sb.append(alignment.numSNPs() - numSNPsCovered() + "\t");
		sb.append(startPos() + "\t");
		sb.append(endPos() + "\t");
		sb.append(length() + "\t");
		sb.append(numPieces() + "\t");
		sb.append(String.format("%.2f", averageCoverage()) + "\t");
		sb.append(this.coverageOfAllSNPs() + "\t");
		//sb.append(averageCoverage() + "\t");
		sb.append(getConsensusAtAllSNPs() + "\t");
		sb.append(consensusHumanReadable() + "\t");
		for(Sequence seq: haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seq = haplotype.get(i);
			sb.append(seq.getName());
			sb.append("(" + seq.getStartPosition() + ")");
			sb.append(',');
			ArrayList<Sequence> masking = seq.getMasking();
			//for(Sequence maskedSeq: seq.getMasking()) {
			for(int j = 0; j < masking.size(); j++) {
				Sequence maskedSeq = masking.get(j);
				sb.append('+');
				sb.append(maskedSeq.getName());
				sb.append("(" + maskedSeq.getStartPosition() + ")");
				sb.append(',');
			}
		}
		sb.deleteCharAt(sb.length()-1);
		return sb;
	}


	/**
	 * Returns the sequences in this haplotype as a list of sequence names, prefixed by a
	 * plus if they are masked, with start position in parens
	 * Returns a '-' if there are no sequences in this haplotype
	 * @return
	 */
	public String sequenceNamesAsString() {
		if(haplotype.size() == 0) return "-";
		StringBuilder sb = new StringBuilder();
		for(Sequence seq : haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seq = haplotype.get(i);
			sb.append(seq.getName());
			sb.append("(" + seq.getStartPosition() + ")");
			sb.append(',');
			ArrayList<Sequence> masking = seq.getMasking();
			//for(Sequence maskedSeq: seq.getMasking()) {
			for(int j = 0; j < masking.size(); j++) {
				Sequence maskedSeq = masking.get(j);
				sb.append('+');
				sb.append(maskedSeq.getName());
				sb.append("(" + maskedSeq.getStartPosition() + ")");
				sb.append(',');
			}
		}
		sb.deleteCharAt(sb.length()-1);
		return sb.toString();
	}

	/**
	 * Returns the number of "pieces" the haplotype is in: the number of connected
	 * components in the graph composed of sequences as nodes and edges where two
	 * sequences cover a non-masked SNP.
	 *
	 * Runs in O(n^2*numSNPs) (assuming the graph library can get connected componens in linear time,
	 * which it damn well should).
	 * @return
	 */
	public int numPieces() {
		
		ArrayList<DisjointSetNode> djNodes = new ArrayList<DisjointSetNode>();
		for(Sequence seq: haplotype) {
		//for(int j = 0; j < haplotype.size(); j++) {
			//Sequence seq = haplotype.get(j);
			if(!seq.isMasked()) {
				djNodes.add(new DisjointSetNode(seq));
			}
		}

		// union the connected components
		//for(DisjointSetNode inode : djNodes) {
		for(int i = 0; i < djNodes.size()-1; i++) {
			DisjointSetNode inode = djNodes.get(i);
			//for(DisjointSetNode jnode: djNodes) {
			for(int j = i+1; j < djNodes.size(); j++) {
				DisjointSetNode jnode = djNodes.get(j);
				Sequence iseq = (Sequence)inode.getData();
				Sequence jseq = (Sequence)jnode.getData();
				if(iseq != jseq && iseq.sharesSNPWith(jseq)) {
					inode.union(jnode);
				}
			}
		}

		HashSet<DisjointSetNode> roots = new HashSet<DisjointSetNode>();
		for(int i = 0; i < djNodes.size(); i++) {
			DisjointSetNode root = djNodes.get(i).find();
			if(!roots.contains(root)) {
				roots.add(root);
			}
		}

		return roots.size();
	}

	/**
	 * Returns the haplotype as an array of contiguous sequences (no '~'s),
	 * complete with the SNPs covering each.
	 *
	 * This MAY return a list of size 0, if there is no sequence in this haplotype
	 * (usually just for the the universal haplotype, which may be empty)
	 * @return
	 * @throws Exception
	 */
	/*public ArrayList<Sequence> getContiguousSeqs() throws Exception {
		Collections.sort(snpsCovered, new SNPPosComparator()); int snpIndex = 0;
		ArrayList<SNP> currentSNPList = new ArrayList<SNP>();

		ArrayList<Sequence> contiguousSeqs = new ArrayList<Sequence>();
		int seqNum = 0;
		StringBuilder fullSeq = new StringBuilder(this.consensusHumanReadable());
		int startPosition = 0; // the consensusHumanReadable given starts at 0
		StringBuilder pieceBuffer = null;
		int pieceIndex = -1;
		for(int i = 0; i < fullSeq.length(); i++) {
			//System.out.println("Num SNPS: " + snpsCovered.size() + " snpIndex: " + snpIndex + " number seqs: " + haplotype.size());
			if(snpsCovered.size() > 0) {
				if(snpsCovered.get(snpIndex).getPosition() == i) {
					currentSNPList.add(snpsCovered.get(snpIndex));
					if(snpIndex < snpsCovered.size()-1) snpIndex = snpIndex + 1;
				}
			}
			if(fullSeq.charAt(i) == '~') {
				if(pieceBuffer != null) { // we had a seq that we are ending
					// add the piece with the appropriate index and start a new one
					//seq.addPiece(pieceBuffer.toString(), pieceIndex+startPosition);
					Sequence newSeq = new Sequence();
					newSeq.addPiece(pieceBuffer.toString(), pieceIndex+startPosition);
					newSeq.setName(name + "_" + seqNum);
					newSeq.addCoveredSNPs(currentSNPList);
					currentSNPList = new ArrayList<SNP>();

					seqNum = seqNum + 1;
					contiguousSeqs.add(newSeq);
					pieceIndex = -1;
					pieceBuffer = null;
				}
				else {
					// do nothing, we're just dropping ~ chars
				}
			}
			else {
				if(pieceBuffer != null) {
					// we're continuing the pieceBuffer
					pieceBuffer.append(fullSeq.charAt(i));
				}
				else {
					// here we're starting a new piece, so we need to record the index
					// and start a new pieceBuffer, not forgetting to add the current char
					pieceBuffer = new StringBuilder();
					pieceIndex = i;
					pieceBuffer.append(fullSeq.charAt(i));
				}
			}
		}
		if(fullSeq.charAt(fullSeq.length()-1) != '~') {
			// if we didnt' end on a ~, we need to add the final crated pieceBuffer
			//seq.addPiece(pieceBuffer.toString(), pieceIndex+startPosition);
			Sequence newSeq = new Sequence();
			newSeq.addPiece(pieceBuffer.toString(), pieceIndex+startPosition);
			newSeq.setName(name + "_" + seqNum);
			newSeq.addCoveredSNPs(currentSNPList);
			currentSNPList = new ArrayList<SNP>();
			
			contiguousSeqs.add(newSeq);
			seqNum = seqNum + 1;
		}
		return contiguousSeqs;
	}*/


	public String consensusTrimmed() {
		int startPos = this.startPos();
		int endPos = this.endPos();
		return this.consensusInRange(startPos, endPos);
	}

	/**
	 * Returns the consensusHumanReadable of the haplotype over the entire multiple alignment.
	 * For columns which the (masked and unmasked) sequences of this haplotype cover,
	 * but are not SNP positions, this is the consensusHumanReadable of the multiple alignment as a whole.
	 * For columsn which are SNPs that are covered, this is the consensusHumanReadable of just this haplotype.
	 *
	 * O(m+n) where m= length of reads n= number of reads
	 * @return
	 */
	public String consensusHumanReadable() {
		//System.err.println("Getting consensusHumanReadable...");
		StringBuffer sb = new StringBuffer();
		char[] consensusArray = new char[alignment.getLength()];
		for(int pos = 0; pos < alignment.getLength(); pos++) {
			if(this.coversPosition(pos)) {
				consensusArray[pos] = alignment.majorityVoteConsensusAtPosition(pos);
				//consensusArray[pos] = this.consensusOfHapAtPosition(pos);
			}
			else {
				consensusArray[pos] = '~';
			}
		}

		for(SNP theSnp : snpsCovered.values()) {
			//theSnp = snpsCovered.get(j);
			int pos = theSnp.getPosition();
			consensusArray[pos] = this.consensusOfHapAtPosition(pos);
		}

		// This takes any internally variant positions and masks them in the consensusHumanReadable. Turning it off cuz it
		// doesn't seem to help.
		/*ArrayList<SNP> withinVariants = allVariantPositionsAsSNPs();
		for(int i = 0; i < withinVariants.size(); i++) {
			int pos = withinVariants.get(i).getPosition();
			consensusArray[pos] = '~';
		}*/

		for(int k = 0; k < consensusArray.length; k++) {
			sb.append(consensusArray[k]);
		}

		return sb.toString();
	}

	/**
	 * Given a range, returns the first covered SNP that occurs in that range.
	 * May return null if no SNPs are covered in the given range.
	 * @param start
	 * @param end
	 * @return
	 */
	public SNP firstCoveredSNPInRange(int start, int end) {
		SNP firstCoveredSNP = null;
		ArrayList<SNP> justSnps = new ArrayList<SNP>();
		justSnps.addAll(snpsCovered.values());
		Collections.sort(justSnps, new SNPPosComparator());

		for(int i = 0; i < justSnps.size(); i++) {
			SNP snpi = justSnps.get(i);
			if(start <= snpi.getPosition() && snpi.getPosition() <= end) {
				firstCoveredSNP = snpi;
				break;
			}
		}

		return firstCoveredSNP;
	}


	/**
	 * Returns the consensusHumanReadable, but only from positions start to end
	 * @param start
	 * @param end
	 * @return
	 */
	public String consensusInRange(int start, int end) {
		//System.err.println("Getting consensusHumanReadable...");
		StringBuffer sb = new StringBuffer();
		char[] consensusArray = new char[end-start+1];
		for(int pos = start; pos <= end; pos++) {
			if(this.coversPosition(pos)) {
				consensusArray[pos-start] = alignment.majorityVoteConsensusAtPosition(pos);
				//consensusArray[pos] = this.consensusOfHapAtPosition(pos);
			}
			else {
				consensusArray[pos-start] = '~';
			}
		}

		for(SNP theSnp : snpsCovered.values()) {
			//SNP theSnp = snpsCovered.get(j);
			int pos = theSnp.getPosition();
			if(start <= pos && pos <= end) {
				consensusArray[pos-start] = this.consensusOfHapAtPosition(pos);
			}
		}

		// This takes any internally variant positions and masks them in the consensusHumanReadable. Turning it off cuz it
		// doesn't seem to help.
		/*ArrayList<SNP> withinVariants = allVariantPositionsAsSNPs();
		for(int i = 0; i < withinVariants.size(); i++) {
			int pos = withinVariants.get(i).getPosition();
			consensusArray[pos] = '~';
		}*/

		for(int k = 0; k < consensusArray.length; k++) {
			sb.append(consensusArray[k]);
		}

		return sb.toString();
	}

	/**
	 * returns true if some sequence (masked or unmasked) covers the position, false otherwise
	 * @param pos
	 * @return
	 */
	public boolean coversPosition(int pos) {
		for(Sequence seq : haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seq = haplotype.get(i);
			if(seq.containsPosition(pos)) return true;
			ArrayList<Sequence> masking = seq.getMasking();
			for(int j = 0; j < masking.size(); j++) {
				Sequence masked = masking.get(j);
				if(masked.containsPosition(pos)) return true;
			}
		}

		return false;
	}

	/*
	 * Returns the number of sequences in this haplotype covering the position, 
	 * including masked sequences
	 */
	public int coverageAtPosition(int pos) {
		int toRet = 0;
		for(Sequence seqi : haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seqi = haplotype.get(i);
			if(seqi.containsPosition(pos)) toRet = toRet + 1;
			ArrayList<Sequence> masked = seqi.getMasking();
			for(int j = 0; j < masked.size(); j++) {
				Sequence maskedj = masked.get(j);
				if(maskedj.containsPosition(pos)) toRet = toRet + 1;
			}
		}

		return toRet;
	}

	/*
	 * Returns the number of bases covering all snp positions in this haplotype
	 */
	public int coverageOfAllSNPs() {
		int covers = 0;
		for(SNP theSnp : snpsCovered.values()) {
			covers = covers + coverageAtPosition(theSnp.getPosition());
		}
		return covers;
	}


	/**
	 * Given a position, returns the consensusHumanReadable of this hap if it is at a SNP position, otherwise
	 * returns the consensusHumanReadable of the overall alignment
	 * @param pos
	 * @return
	 */
	public char consensusOverallAtPosition(int pos) {
		if(snpsCovered.containsKey(pos)) {
			return this.consensusOfHapAtPosition(pos);
		}
		else {
			return alignment.majorityVoteConsensusAtPosition(pos);
		}
	}


	/**
	 * returns the majority vote, amongst masked and non-masked sequences in THIS
	 * haplotype. Note that this will be true even if you ask about a position not
	 * covered by a SNP (and thus is likely to contain a sequencing error)
	 * @param pos
	 * @return
	 */
	public char consensusOfHapAtPosition(int pos) {
		HashMap<Character, Integer> baseCountHash = new HashMap<Character, Integer>();
		for(Sequence seq: haplotype) {
		//for(int j = 0; j < haplotype.size(); j++) {
			//Sequence seq = haplotype.get(j);
			//if(i == 0) System.out.println("Considering sequence " + seq.getName() + " masked is: " + seq.isMasked());
			//if(seq.isMasked()) System.out.println("this shouldnt be masked...");
			char baseAti = seq.baseAtAlignmentPosition(pos);
			if(!baseCountHash.containsKey(baseAti)) {
				baseCountHash.put(baseAti, 1);
			}
			else {
				baseCountHash.put(baseAti, baseCountHash.get(baseAti)+1);
			}
			ArrayList<Sequence> masking = seq.getMasking();
			for(int k = 0; k < masking.size(); k++) {
				Sequence seqk = masking.get(k);
				//if(i == 0) System.out.println("Considering masked sequence " + seqk.getName() + " which is masked by " + seq.getName());
				char baseAtik = seqk.baseAtAlignmentPosition(pos);
				if(!baseCountHash.containsKey(baseAtik)) {
					baseCountHash.put(baseAtik, 1);
				}
				else {
					baseCountHash.put(baseAtik, baseCountHash.get(baseAtik)+1);
				}
				//System.out.println("considering a masked sequence... "+ " base at pos " + i + " is " + baseAti + " current count is " + baseCountHash.get(baseAti));

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

	public String getConsensusAtAllSNPs() {
		StringBuffer sb = new StringBuffer();
		HashSet<SNP> allSnps = alignment.getSNPs();
		ArrayList<SNP> allSnpsList = new ArrayList<SNP>(allSnps);
		Collections.sort(allSnpsList, new SNPPosComparator());

		if(allSnps.size() == 0) {
			sb.append("!");
		}
		
		//for(SNP theSNP : allSnps) {
		for(int k = 0; k < allSnpsList.size(); k++) {
			SNP theSNP = allSnpsList.get(k);
			int i = theSNP.getPosition();
			HashMap<Character, Integer> baseCountHash = new HashMap<Character, Integer>();
			for(Sequence seq: haplotype) {
			//for(int j = 0; j < haplotype.size(); j++) {
				//Sequence seq = haplotype.get(j);
				char baseAti = seq.baseAtAlignmentPosition(i);
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
			sb.append(maxBase);
		}

		return sb.toString();
	}

	/**
	 * Returns the number of unmasked SNPs that are covered which are showing inconsistencies.
	 * @return
	 */
	public int numInconsistentSNPs() {
		int count = 0;
		for(SNP theSnp : snpsCovered.values()) {
			//SNP theSnp = snpsCovered.get(i);
			// for each unmasked snp...
			//if(!theSnp.isMasked()) {
				int position = theSnp.getPosition();

				char firstSeen = '~';
				for(Sequence seq : haplotype) {
				//for(int j = 0; j < haplotype.size(); j++) {
					//Sequence seq = haplotype.get(j);
					char c = seq.baseAtAlignmentPosition(position);

					if(c != '~' && firstSeen == '~') {
						firstSeen = c;
					}
					else if(c != firstSeen && c != '~') {
						count = count + 1;
						break;
					}
				}
			//}
		}
		return count;
	}


	/**
	 * Returns an arraylist of all the consistent snps that are covered.
	 * @return
	 */
	public ArrayList<SNP> consistentSNPsCovered() {
		ArrayList<SNP> consistentSnps = new ArrayList<SNP>();
		for(SNP theSnp : snpsCovered.values()) {
			//SNP theSnp = snpsCovered.get(i);
			// for each unmasked snp...
			//if(!theSnp.isMasked()) {
				int position = theSnp.getPosition();
				char firstSeen = '~';
				boolean consistent = true;
				for(Sequence seq : haplotype) {
				//for(int j = 0; j < haplotype.size(); j++) {
					//Sequence seq = haplotype.get(j);
					char c = seq.baseAtAlignmentPosition(position);

					if(c != '~' && firstSeen == '~') {
						firstSeen = c;
					}
					else if(c != firstSeen && c != '~') {
						consistent = false;
					}
				}
				if(consistent) {
					consistentSnps.add(theSnp);
				}
			//}
		}
		return consistentSnps;
	}



/**
	 * Returns an arraylist of SNPs, where each SNP describes the position
	 * of some variant position _within_ this haplotype, regardless of whether
	 * it coincides with a SNP or not, including masked sequences. Useful for debugging.
	 * Runs in Theta(alignment.length() * alignment.getNumReads)
	 * @return
	 */
	public ArrayList<SNP> allVariantPositionsAsSNPs() {
		ArrayList<SNP> snpList = new ArrayList<SNP>();
		int alignmentLength = alignment.getLength();

		for(int pos = 0; pos < alignmentLength; pos++) {
			// For each column, if there is more than one type of character that is
			// not a '~' this column is a SNP.
			char[] col = this.hapColumnAtPosition(pos);

			boolean isSNP = false;
			char firstSeen = '~';
			//for(char c : col) {
			for(int i = 0; i < col.length; i++) {
				char c = col[i];
				//System.out.println("firstSeen is " + firstSeen + " and current char is " + c);
				//if(c != '~') System.out.println("currently at position " + pos + ", i see a " + c);
				if(c != '~' && firstSeen == '~') {
					firstSeen = c;
				}
				else if(c != firstSeen && c != '~') {
					//System.out.println("position " + pos + " is a snp");
					isSNP = true;
					//System.out.println("#########first seen was "+ firstSeen +" and now I see a " + c + ", position is " + pos);
				}
			}
			if(isSNP) {
				SNP newSNP = new SNP(pos, alignment);
				snpList.add(newSNP);
			}
		}

		return snpList;
	}



	/**
	 * Returns a column representing the alleles at the given position of all sequences in this
	 * haplotype, including masked sequences. Useful for debugging.
	 * @param pos
	 * @return
	 */
	public char[] hapColumnAtPosition(int pos) {
		int numSeqs = this.numNonRedundantSequences() + this.numRedundantSequences();
		char[] column = new char[numSeqs];
		//System.out.print("hap " + this.getName() + " has " + numSeqs + " total Seqs. Alleles at " + pos + ": ");
		int i = 0;
		for(Sequence seq : haplotype) {
		//for(int j = 0; j < haplotype.size(); j++) {
			//Sequence seq = haplotype.get(j);
			if(!seq.isMasked()) {
				column[i] = seq.baseAtAlignmentPosition(pos);
				//System.out.print(i+" : "+column[i] + "   ");
				i = i + 1;
				ArrayList<Sequence> masking = seq.getMasking();
				for(int k = 0; k < masking.size(); k++) {
					Sequence masked = masking.get(k);
					column[i] =  masked.baseAtAlignmentPosition(pos);
					//System.out.print(i+" : "+column[i] + "   ");
					i = i + 1;
				}
			}
		}
		//System.out.println();
		return column;
	}




	/**
	 * Returns true if theSeq is consistent with all sequences in this haplotype at the given SNP positions
	 * (regardless of the maskedness of the snps in the snp list.)
	 * @param theSeq
	 * @param snpList
	 * @return
	 */
	public boolean sequenceConsistentAtSNPs(Sequence theSeq, ArrayList<SNP> snpList) {
		for(Sequence hapSeq : haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence hapSeq = haplotype.get(i);
			if(hapSeq.conflictsWithAtAnySNP(theSeq, snpList)) {
				return false;
			}
		}
		return true;
	}


	/**
	 * Returns the number of mismatches from a given sequence to the consensusHumanReadable of this haplotype
	 * ( the REPORTED consensusHumanReadable = consensusHumanReadable of this
	 * haplotype at SNP positions, of the alignment overall at other positions)
	 * if the given sequence is longer than the multiple alignment, we don't consider that a "mismatch"
	 * @param theSeq
	 * @param pos
	 * @return
	 */
	public int sequenceMismatchCount(Sequence theSeq) {
		// the consensusHumanReadable _might_ be shorter than the compared sequence if we are comparing to
		// a ground truth sequence. So we gotta make sure we are only checking the area of the multiple alignment
		String consensus = this.consensusHumanReadable();
		/*System.out.println("**************");
		System.out.println(consensusHumanReadable);
		System.out.println(theSeq.sequenceToString());
		System.out.println("Start is: " + start + " end is: " + end);
		System.out.println("**************");*/
		int count = 0;
		for(int pos = 0; pos < alignment.getLength(); pos++) {
			char consensusChar = consensus.charAt(pos);
			char seqChar = theSeq.baseAtAlignmentPosition(pos);
			if(seqChar != '~' && consensusChar != '~' && seqChar != consensusChar) count = count + 1;
		}
		return count;
	}

	/**
	 * Returns the number of non-masked sequences in this haplotype
	 * @return
	 */
	public int numNonRedundantSequences() {
		return haplotype.size();
	}

	/**
	 * Returns number of sequences that are being masked by sequences in this haplotype
	 * TODO: how do I make sure that multiple sequences don't mask the same sequence?
	 * @return
	 */
	public int numRedundantSequences() {
		int numRedundant = 0;
		for(Sequence seq: haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seq = haplotype.get(i);
			numRedundant = numRedundant + seq.getMasking().size();
		}
		return numRedundant;
	}

	/**
	 * Returns the number of non-masked SNPs sequences in this haplotype are covering.
	 * Works in O(n*q), where n is the number of sequences and q is the number of SNPs covered by
	 * all sequences (in the worst case when all sequenes cover the same SNPS).
	 * @return
	 */
	public int numSNPsCovered() {
		return snpsCovered.size();
	}

	/**
	 * Returns the starting position within the multiple alignment of the first sequence in the haplotype
	 * Runs in O(n), n = number of sequences in the alignment
	 * Returns -1 if there are no sequences in this haplotype.
	 * @return
	 */
	public int startPos() {
		if(haplotype.size() == 0) {
			return -1;
		}

		int startPos = alignment.getLength();

		for(Sequence seq: haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seq = haplotype.get(i);
			if(seq.getStartPosition() < startPos) {
				startPos = seq.getStartPosition();
			}
			ArrayList<Sequence> masking = seq.getMasking();
			for(int j = 0; j < masking.size(); j++) {
				Sequence masked = masking.get(j);
				if(masked.getStartPosition() < startPos) {
					startPos = masked.getStartPosition();
				}
			}
		}

		return startPos;
	}


	/**
	 * Returns the ending position within the multiple alignment of the last sequence in the haplotype
	 * Runs in O(n), n = number of seqs
	 * Returns -1 if there are no sequences in this haplotype
	 * @return
	 */
	public int endPos() {
		int endPos = -1;

		for(Sequence seq: haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence seq = haplotype.get(i);
			if(seq.endPosition() > endPos) {
				endPos = seq.endPosition();
			}
			ArrayList<Sequence> masking = seq.getMasking();
			for(int j = 0; j < masking.size(); j++) {
				Sequence masked = masking.get(j);
				if(masked.endPosition() > endPos) {
					endPos = masked.endPosition();
				}
			}
		}

		return endPos;
	}

	/**
	 * Returns an arrayList of _all_ sequences in this haplotype;
	 * ie, the sequences defining the haplotype and those that
	 * are covered.
	 * @return
	 */
	public ArrayList<Sequence> allSequences() {
		ArrayList<Sequence> allSeqs = new ArrayList<Sequence>();
		for(Sequence hapSeq : haplotype) {
		//for(int i = 0; i < haplotype.size(); i++) {
			//Sequence hapSeq = haplotype.get(i);
			ArrayList<Sequence> masking = hapSeq.getMasking();
			allSeqs.add(hapSeq);
			allSeqs.addAll(masking);
		}

		return allSeqs;
	}

	/**
	 * Returns the average coverage of the haplotype, from first non-~ character
	 * to last non-~ character, including masked sequences. Note that the since
	 * a haplotype might have ~'s a consensusHumanReadable in the middle, these columns have coverage
	 * 0
	 * @return
	 */
	public double averageCoverage() {
		int count = 0;
		ArrayList<Sequence> allSeqs = allSequences();
		if(allSeqs.size() == 0) return 0.0;
		for(int i = 0; i < allSeqs.size(); i++) {
			count = count + allSeqs.get(i).numBases();
		}

		return (double)count/(double)length();
	}

	/**
	 * Returns the length of the consensusHumanReadable, from first non-~ base to last.
	 * If there are no sequences in this haplotype, returns 0
	 * @return
	 */
	public int length() {
		if(haplotype.size() == 0) {
			return 0;
		}
		return endPos() - startPos() + 1;
	}

	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	/**
	 * returns the full name of this haplotype, as ALIGNMENTNAME_BLOCKNAME_HAPNAME
	 * @return
	 */
	public String fullName() {
		StringBuilder sb = new StringBuilder();
		HaplotypeBlock myBlock = getHaploBlock();
		MultipleAlignment myAlignment = myBlock.getAlignment();
		String myAlName = myAlignment.getName();
		sb.append(myAlName);
		sb.append("_");
		sb.append(getHaploBlock().getName());
		sb.append("_");
		sb.append(getName());

		return sb.toString();
	}

	public HaplotypeBlock getHaploBlock() {
		return haploBlock;
	}

	public void setHaploBlock(HaplotypeBlock haploBlock) {
		this.haploBlock = haploBlock;
	}



}
