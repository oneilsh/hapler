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

import java.util.Comparator;

/**
 *
 * @author soneil
 */
public class SequenceFirstSNPLengthComparator implements Comparator<Sequence> {	/**

 * Returns which sequence contains the first SNP. If they share the same first SNP, returns the longer one.
	 * If they share the same first SNP and are the same length, they are "equal"
	 * If one of them contains no SNPS, it is first. If neither contain no SNPs, then the longer one is first.
	 * If neither contain snps and they are the same length... they are equal
	 * @param seqi
	 * @param seqj
	 * @return
	 */
	public int compare(Sequence seqi, Sequence seqj) {
		SNP firstSNPi = seqi.firstSNP();
		SNP firstSNPj = seqj.firstSNP();
		// if either are null
		if(firstSNPi == null || firstSNPj == null) {
			if(firstSNPi == null && firstSNPj == null) { // both null
				if(seqi.length() > seqj.length()) { // decide based on length
					return -1;
				}
				else if(seqi.length() < seqj.length()) {
					return 1;
				}
				else {
					return 0;
				}
			}
			else if(firstSNPi == null) { // just i is null, so it comes first
				return -1;
			}
			else { // just j is null, so it comes first
				return 1;
			}
		}

		// ok, neither of them are null
		if(firstSNPi.getPosition() < firstSNPj.getPosition()) {
			return -1;
		}
		else if(firstSNPi.getPosition() > firstSNPj.getPosition()) {
			return 1;
		}
		else {
			if(seqi.length() > seqj.length()) {
				return -1;
			}
			else if(seqi.length() < seqj.length()) {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
}
