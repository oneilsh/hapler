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

/**
 *
 * @author soneil
 */
public class HaplotypeBlockGreedyComparator {


	/**
	 * Returns which haplotype block is "better"
	 * First, checks for consistency: if A has fewer inconsistencies than B, A is better
	 * Next, checks for number of pieces: if A has fewer overall pieces than B, A is better
	 * Next, checks for number of haplotypes (although this SHOULD be unnecessary in the contexts in which I plan on using it)
	 * Finally, checks for greedyness: The lexicographically larger is preferred.
	 * All these being equal, the HaplotypeBlocks are equal.
	 * @param seqi
	 * @param seqj
	 * @return
	 */
	public int compare(HaplotypeBlock hapa, HaplotypeBlock hapb) {
		int aConflicts = hapa.numInconsistentSNPs();
		int bConflicts = hapb.numInconsistentSNPs();
		if(aConflicts < bConflicts) {
			return -1;
		}
		else if(aConflicts > bConflicts) {
			return 1;
		}
		else {
			int aPieces = hapa.numPieces();
			int bPieces = hapb.numPieces();
			if(aPieces < bPieces) {
				return -1;
			}
			else if(bPieces < aPieces) {
				return 1;
			}
			else {
				ArrayList<Double> aLexicographicSize = hapa.lexicographicSize();
				ArrayList<Double> bLexicographicSize = hapb.lexicographicSize();
				if(aLexicographicSize.size() < bLexicographicSize.size()) {
					return -1;
				}
				else if(aLexicographicSize.size() > bLexicographicSize.size()) {
					return 1;
				}
				else {
					for(int i = aLexicographicSize.size() - 1; i >= 0; i--) {
						Double aVal = aLexicographicSize.get(i);
						Double bVal = bLexicographicSize.get(i);
						if(aVal > bVal) {
							return -1;
						}
						else if(aVal < bVal) {
							return 1;
						}
					}
					return 0;
				}
			}
		}
	}
}
