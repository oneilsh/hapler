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
public class BasicAlignmentSNPListCreator extends AbstractSNPListCreator {


	public BasicAlignmentSNPListCreator() {
		// Empty constructor
	}



	/**
	 * Returns an ArrayList of SNPs, with SNPs corresponding to any column of the
	 * multiple alignment that is variant (has more than one letter other than ~)
	 * Runs in Theta(alignment.getLength() * alignment.getNumReads)
	 * @return
	 */
	public ArrayList<SNP> computeSNPList(MultipleAlignment alignment) {
		ArrayList<SNP> snpList = new ArrayList<SNP>();
		int alignmentLength = alignment.getLength();

		for(int pos = 0; pos < alignmentLength; pos++) {
			// For each column, if there is more than one type of character that is
			// not a '~' this column is a SNP.
			char[] col = alignment.unmaskedColumn(pos);

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





}
