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
import java.util.HashMap;

/**
 * This class returns a list of SNPs based on a multiple alignment, returning any
 * variant position, but not if the only positions are '-' and some base, and those
 * sequences having the base are part of a sequence of 3 or more of that base.
 * If the multiple alignment's length is less than or equal to 7, will return an emtpy list of SNPs no matter what.
 * Not necessarily very fast (probably two to five times slower than it could be),
 * Still, Theta(numReads * alignmetnLength)
 *
 * Eg:
 *                     1 1 1
 * 0 1 2 3 4 5 6 7 8 9 0 1 2
 * A A A T T T G G C C T T T
 * A - A T T T G G C C T T T
 * - A A T - T G - C C - T T
 * A A - - T T G G C C T T -
 *
 * In the above, the only "real" SNP position is position 7, all others look like homopolymer runs
 * @author soneil
 */
public class Basic454SNPListCreator extends AbstractSNPListCreator {




	public Basic454SNPListCreator() {
		// Empty constructor
	}



	public ArrayList<SNP> computeSNPList(MultipleAlignment alignment) {
		ArrayList<SNP> snpList = new ArrayList<SNP>();
		if(alignment.getLength() <= 7) {
			return snpList;
		}
		//Is the first column a snp?
		String[] first3Cols = alignment.unmaskedColumns(0, 2);
		if(baseInAlignmentChunk454Variant(first3Cols, 0)) snpList.add(new SNP(0, alignment));
		//Is the second column a snp?
		String[] first4Cols = alignment.unmaskedColumns(0, 3);
		if(baseInAlignmentChunk454Variant(first4Cols, 1)) snpList.add(new SNP(1, alignment));
		//Is the last Column a snp?
		String[] last3Cols = alignment.unmaskedColumns(alignment.getLength()-3, alignment.getLength()-1);
		if(baseInAlignmentChunk454Variant(last3Cols, 2)) snpList.add(new SNP(alignment.getLength()-1, alignment));
		//Is the second to last Column a snp?
		String[] last4Cols = alignment.unmaskedColumns(alignment.getLength()-4, alignment.getLength()-1);
		if(baseInAlignmentChunk454Variant(last4Cols, 2)) snpList.add(new SNP(alignment.getLength()-2, alignment));

		//are any of the columns in the middle snps?
		for(int i = 2; i <= alignment.getLength()-3; i++) {
			String[] fiveCols = alignment.unmaskedColumns(i-2, i+2);
			if(baseInAlignmentChunk454Variant(fiveCols,2)) snpList.add(new SNP(i, alignment));
		}


		return snpList;
	}


	private boolean baseInAlignmentChunk454Variant(String[] seqs, int pos) {

		// This keys in this hash correspond to the possible bases at position pos
		// in the multiple alignment. If every instance of that base in the
		// strings at that position are part of homopolymer runs, the associated
		// boolean will be TRUE. else, it will be false.
		//
		// Thus, if there is only one non-'~' non-'-' key in the hash, and it's
		// boolean is true, it's a 454 error and not a snp (since all instances of this base
		// existing are part of homopolymer runs)
		// TODO: BUG: man, this is fucking BROKEN. A hashmap is so wrong, we want to be counting alleles and such,
		// which a hashmap won't let us do
		// TODO: IS this still broken? I don't think so.
		HashMap<Character,Integer> baseCounts = new HashMap<Character,Integer>();
		boolean allBaseTypesHomoRuns = true;


		//for(String str : seqs) {
		for(int i = 0; i < seqs.length; i++) {
			String str = seqs[i];
			//System.out.println("str is " + str + " pos is " + pos + " i is " + i);
			char base = str.charAt(pos);
			boolean homoRunHere = baseInStringPartOfHomoPolymerRun(str, pos);
			if(base != '-' && base != '~' && !homoRunHere) allBaseTypesHomoRuns = false;
			if(!baseCounts.containsKey(base)) baseCounts.put(base, 1);
			else baseCounts.put(base, baseCounts.get(base)+1);
			//else if(!homoRunHere) baseCounts.put(base, false);
		}

		int numGaps = 0;
		int numNonGaps = 0;
		int numBaseTypes = 0;
		//Character[] baseToHomoRunsKeySet = (Character[])baseCounts.keySet().toArray();
		for(char base : baseCounts.keySet()) {
		//for(int i = 0; i < baseToHomoRunsKeySet.length; i++) {
			//char base = baseToHomoRunsKeySet[i];
			//boolean allHomoRuns = baseCounts.get(base);
			if(base != '-' && base != '~') numBaseTypes = numBaseTypes + 1;
			if(base == '-') numGaps = baseCounts.get(base);
			else if(base != '~') numNonGaps = numNonGaps + baseCounts.get(base);
			// only say something isn't a homopolymer run if it's an actual base
			//if(!allHomoRuns && base != '~' && base != '-') allBaseTypesHomoRuns = false;
		}

		// 
		// This will be called a snp when it shouldn't be:
		// TTT-AAAA
		// TTTTAAAA
		// TTT-AAAA
		// TTTAAAAA
		// TTT-AAAA
		// FIXED By adding !allBaseTypesHomoRuns to first check

		if(numBaseTypes > 1 && numNonGaps > numGaps && !allBaseTypesHomoRuns) return true;
		// if numBaseTypes is 1 (or godforbid 0), it's a SNP if any of them DON'T look like a homopolymer run
		if(!allBaseTypesHomoRuns && numBaseTypes == 1 && baseCounts.containsKey('-') && numNonGaps > numGaps) return true;
		return false;
	}

	/**
	 * Given a portion of a multiple alignment, tells us whether position pos
	 * in the chunk of multiple alignment is variant (ie, contains two characters
	 * other than ~).
	 * 0: not variant at all
	 * 1: gap-variant
	 * 2: two or more bases seen
	 * @param columns
	 * @param pos
	 * @return
	 */
	/*private int baseInAlignmentChunkIsVariant(String[] columns, int pos) {
		HashMap<Character,Integer> baseCountsHash = new HashMap<Character,Integer>();
		for(String str : columns) {
			char base = str.charAt(pos);
			if(!baseCountsHash.containsKey(base)) baseCountsHash.put(base, 1);
			else baseCountsHash.put(base, baseCountsHash.get(base) + 1);
		}
		if(baseCountsHash.size() == 0) return -1; // shouldn't happen
		if(baseCountsHash.containsKey('~')) {
			if(baseCountsHash.containsKey('-')) {
				// ~ and - and possibly other stuff
				if(baseCountsHash.size() == 2) return 0; // just gaps, shouldn't happen unless our alignment is stupid
				else if(baseCountsHash.size() == 3) return 1; // gaps and one other base, gap-variant
				else return 2; // at least two other things, really variant
			}
			else {
				// ~ but no gaps, possibly other stuff
				if(baseCountsHash.size() == 1) return 0; // only ~'s, shouldn't happen unless our alignment is stupid
				else if (baseCountsHash.size() == 2) return 0; // not variant at all, only one type of base
				else return 2; // at least two other things, really variant
			}
		}
		else {
			if(baseCountsHash.containsKey('-')) {
				// - but no ~'s
				if(baseCountsHash.size() == 1) return 0; // only -'s shouldnt happen unless our alignment is stupid
				else if(baseCountsHash.size() == 2) return 1; // only one other base, gap-variant
				else return 2; // at least two other bases, really variant
			}
			else {
				// no -s or ~s, everything we saw was a base
				if(baseCountsHash.size() == 1) return 0; // all the same base
				else return 2; // at least two other things, really variant
			}
		}

	}*/


	/**
	 * Given a string and a position in the string, tells us whether the base at that
	 * position is part of a homopolymer run (ie, is part of a run of 3 or more of that base)
	 * Returns false if the position is ~ or -.
	 * @param str
	 * @param pos
	 * @return
	 */
	private boolean baseInStringPartOfHomoPolymerRun(String str, int pos) {
		char centralBase = str.charAt(pos);
		// return false if the central base is a '-' or a '~' : these can't be homopolymer runs
		if(centralBase == '-' || centralBase == '~') return false;
		// is it a run starting at pos?
		if(pos <= str.length() - 3 && str.charAt(pos + 1) == centralBase && str.charAt(pos+2) == centralBase) return true;
		// is it a run ending at pos?
		if(pos >= 2 && str.charAt(pos-1) == centralBase && str.charAt(pos-2) == centralBase) return true;
		// is it a run with the central base in the middle?
		if(pos <= str.length() - 2 && pos >= 1 && str.charAt(pos-1) == centralBase && str.charAt(pos+1) == centralBase) return true;
		return false;
	}

	
}
