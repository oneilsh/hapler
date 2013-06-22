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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

/**
 *
 * @author soneil
 */
public class HaplotypeBlock {

	private ArrayList<Haplotype> haplotypes;
	private MultipleAlignment alignment;

	public MultipleAlignment getAlignment() {
		return alignment;
	}
	private String name;
	private int minColors;
	private int repsToFinish;





	public HaplotypeBlock(MultipleAlignment theAlignment) {
		haplotypes = new ArrayList<Haplotype>();
		alignment = theAlignment;
		minColors = 0;
		repsToFinish = 0;
	}

	public ArrayList<Haplotype> getHaplotypes() {
		return haplotypes;
	}

	public void setHaplotypes(ArrayList<Haplotype> haplotypes) {
		this.haplotypes = haplotypes;
	}

	public void addHaplotype(Haplotype haplotype) {
		haplotypes.add(haplotype);
	}

	public String toString() {
		int alignmentLength = alignment.getLength();
		int maxSeqNameLength = alignment.getMaxSeqNameLength();

		String hapSep = stringOfChars(alignmentLength + maxSeqNameLength + 5, ':');

		StringBuffer sb = new StringBuffer();
		//for(Haplotype haplotype: haplotypes) {
		for(int i = 0; i < haplotypes.size(); i++) {
			Haplotype haplotype = haplotypes.get(i);
			sb.append(hapSep);
			sb.append(System.getProperty("line.separator"));
			sb.append(System.getProperty("line.separator"));
			sb.append(haplotype.toString());
			sb.append(System.getProperty("line.separator"));
		}

		return sb.toString();
	}



	/**
	 * Returns an ArrayList of strings describing each haplotype, with a unique number
	 * before each indicating it's number within this haplotype block.
	 * TODO: This is slow, and and sucky. 
	 * @return
	 */
	public ArrayList<StringBuilder> summaryStrings() {
		ArrayList<StringBuilder> summaryStrings = new ArrayList<StringBuilder>();

		for(int i = 0; i < haplotypes.size(); i++) {
			StringBuilder thisBuilder = new StringBuilder();
			Haplotype hap = haplotypes.get(i);
			thisBuilder.append(this.getMinColors() + "\t");
			thisBuilder.append(this.getRepsToFinish() + "\t");
			thisBuilder.append(hap.getName());
			thisBuilder.append("\t");
			thisBuilder.append(hap.summaryString().toString());
			summaryStrings.add(thisBuilder);
		}

		return summaryStrings;
	}


	/**
	 * Returns all the contiguous haplotype sequences in all haplotypes in this block.
	 * Complete with SNPs covering each.
	 * This MAY return a list with 0 elements (in the case of empty haplotype blocks, which
	 * should only be the universal block.)
	 * @return
	 */
	/*public ArrayList<Sequence> getContiguousSeqs() throws Exception {
		ArrayList<Sequence> toRet = new ArrayList<Sequence>();
		for(int i = 0; i < haplotypes.size(); i++) {
			ArrayList<Sequence> theseSeqs = haplotypes.get(i).getContiguousSeqs();
			for(int j = 0; j < theseSeqs.size(); j++) {
				toRet.add(theseSeqs.get(j));
			}
		}

		return toRet;
	}*/

	/**
	 * Returns all the subhaplotypes defined by each haplotype in this block.
	 * This MAY return a list with 0 elements (in the case of and empty haplotype block,
	 * which should only be the universal block).
	 * @return
	 */
	/*public ArrayList<Haplotype> getSubHaplotypes() {
		ArrayList<Haplotype> toRet = new ArrayList<Haplotype>();
		for(int i = 0; i < haplotypes.size(); i++) {
			ArrayList<Haplotype> theseSubHaps = haplotypes.get(i).getSubHaplotypes();
			for(int j = 0; j < theseSubHaps.size(); j++) {
				Haplotype subHap = theseSubHaps.get(j);
				subHap.setName(alignment.getName() + "_" + this.getName() + "_" + subHap.getName());
				toRet.add(subHap);
			}
		}
		return toRet;
	}*/


	private String stringOfChars(int length, char theChar) {
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < length; i++) {
			sb.append(theChar);
		}
		return sb.toString();
	}

	/**
	 * Returns the sum of numPieces of each haplotype
	 * @return
	 */
	public int numPieces() {
		int numPieces = 0;
		for(int i = 0; i < haplotypes.size(); i++) {
			Haplotype haplotype = haplotypes.get(i);
			numPieces = numPieces + haplotype.numPieces();
		}
		return numPieces;
	}

	/**
	 * Returns the total number of times haplotypes are inconsistent at unmasked SNPs.
	 * For example, suppose HAP1 is inconsistent at SNPS A, B, and D,
	 * And HAP2 is inconsistent at SNPS B, D, and G, then the total number of inconsistent SNPS is 6.
	 * @return
	 */
	public int numInconsistentSNPs() {
		int count = 0;
		for(int i = 0; i < haplotypes.size(); i++) {
			Haplotype hap = haplotypes.get(i);
			count = count + hap.numInconsistentSNPs();
		}
		return count;
	}

	/**
	 * Returns an array of the number of non-redundant sequences in each haplotype,
	 * ordered in increasing order by the number of non-redundant sequences.
	 * EG: [2,2,3,3,4] indicates 5 haplotypes, one with 4 non-redundant seqs, two with 3, and two with 2
	 * @return
	 */
	/*public int[] lexicographicSize() {
		int[] sizesArray = new int[haplotypes.size()];
		for(int i = 0; i < haplotypes.size(); i++) {
			Haplotype haplotype = haplotypes.get(i);
			sizesArray[i] = haplotype.numNonRedundantSequences();
		}
		Arrays.sort(sizesArray);
		return sizesArray;
	}*/


	/*
	 * Returns an array of "support" for each haplotype, sorted.
	 * Each sequence in the haplotype block may "support" (not conflict with) multiple haplotypes;
	 * if a sequence supports n haplotypes, it contributes 1/n to each. Note that support is only computed over
	 * SNPs that are consistent within each haplotype.
	 *
	 * For example, if a haplotype is inconsistent at a SNP B for HAP1, but consistent at HAP2, then
	 * B is ignored when computing support for HAP1, but considered for HAP2.
	 */
	public ArrayList<Double> lexicographicSize() {
		//Double[] supports = haplotypeSupports().values().toArray(new Integer[0]);
		ArrayList<Double> supports = new ArrayList<Double>();

		Collection<Double> doubleCollection = haplotypeSupports().values();
		for(Double value: doubleCollection) {
			supports.add(value);
		}
		Collections.sort(supports);
		
		return supports;
	}

	/**
	 * Returns a hashmap giving the numerical support of each haplotype.
	 * EG: suppose seq1, seq2, seq3, and seq4 are sequences in this haplotype block (which may be masked/covered)
	 * and hapA and hapB are haplotypes. if seq1 supports A and seq2 supports B and seq3 supports A and B and seq4
	 * supports B, then A will have a support of 1.5 and B will have support of 2.5.
	 * @return
	 */
	public HashMap<Haplotype,Double> haplotypeSupports() {
		ArrayList<Sequence> allSeqs = this.allSequences();
		HashMap<Sequence,ArrayList<Haplotype>> seqsToHapsSupporting = seqsToHapsSupporting(allSeqs);

		HashMap<Haplotype,Double> haplotypeSupports = new HashMap<Haplotype,Double>();

		for(int i = 0; i < haplotypes.size(); i++) {
			Haplotype hap = haplotypes.get(i);
			haplotypeSupports.put(hap, 0.0);
		}

		for(int i = 0; i < allSeqs.size(); i++) {
			Sequence seq = allSeqs.get(i);
			ArrayList<Haplotype> supportingList = seqsToHapsSupporting.get(seq);
			double numSupporting = supportingList.size();
			for(int j = 0; j < supportingList.size(); j++) {
				Haplotype supportedHap = supportingList.get(j);
				haplotypeSupports.put(supportedHap, haplotypeSupports.get(supportedHap) + 1.0/numSupporting);
			}


		}

		return haplotypeSupports;
	}

	/**
	 * Given a list of sequences, returns a hashmap which maps each sequence to the haplotypes that it is supporting.
	 * Note that a sequence supports a haplotype if it doesn't conflict at any unmasked, non-conflicting snps.
	 * @param theSeqs
	 * @return
	 */
	private HashMap<Sequence,ArrayList<Haplotype>> seqsToHapsSupporting(ArrayList<Sequence> theSeqs) {

		HashMap<Sequence,ArrayList<Haplotype>> seqsToHapsSupporting = new HashMap<Sequence,ArrayList<Haplotype>>();

		for(int i = 0; i < theSeqs.size(); i++) {
			Sequence seq = theSeqs.get(i);
			for(int j = 0; j < haplotypes.size(); j++) {
				Haplotype hap = haplotypes.get(j);
				ArrayList<SNP> consistentSnps = hap.consistentSNPsCovered();
				if(hap.sequenceConsistentAtSNPs(seq, consistentSnps)) {
					if(!seqsToHapsSupporting.containsKey(seq)) {
						ArrayList<Haplotype> newList = new ArrayList<Haplotype>();
						newList.add(hap);
						seqsToHapsSupporting.put(seq, newList);
					}
					else {
						ArrayList<Haplotype> theList = seqsToHapsSupporting.get(seq);
						theList.add(hap);
					}
				}
			}
		}
		return seqsToHapsSupporting;
	}

	/**
	 * Returns an arraylist of all sequences (including covered/masked ones) from all
	 * haplotypes in this haplotype block.
	 * @return
	 */
	public ArrayList<Sequence> allSequences() {
		ArrayList<Sequence> allSeqs = new ArrayList<Sequence>();
		for(int i = 0; i < haplotypes.size(); i++) {
			Haplotype hap = haplotypes.get(i);
			allSeqs.addAll(hap.allSequences());
		}

		return allSeqs;
	}


	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getMinColors() {
		return minColors;
	}

	public void setMinColors(int minColors) {
		this.minColors = minColors;
	}


	public int getRepsToFinish() {
		return repsToFinish;
	}

	public void setRepsToFinish(int repsToFinish) {
		this.repsToFinish = repsToFinish;
	}



}
