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

/**
 *
 * @author soneil
 */
public class BasicStrictAlignmentSNPListCreator extends AbstractSNPListCreator {


	public BasicStrictAlignmentSNPListCreator() {
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


			HashMap<Character, Integer> alleleCounts = new HashMap<Character, Integer>();
			//for(char c : col) {
			for(int i = 0; i < col.length; i++) {
				char c = col[i];
				//System.out.println("firstSeen is " + firstSeen + " and current char is " + c);
				if(!alleleCounts.containsKey(c)) {
					alleleCounts.put(c, 1);
				}
				else {
					alleleCounts.put(c, alleleCounts.get(c)+1);
				}
			}

			alleleCounts.put('~', 0);

			ArrayList<Integer> mapValues = new ArrayList(alleleCounts.values());
			Collections.sort(mapValues);
			Collections.reverse(mapValues);

			int totalCoverage = 0;
			for(int k = 0; k < mapValues.size(); k++) {
				totalCoverage = totalCoverage + mapValues.get(k);
			}

			/*if(mapValues.size() >= 2) {
				if(mapValues.get(1) >= 0.25*totalCoverage && totalCoverage >= 6) {
					snpList.add(new SNP(pos));
				}
			}*/			
			if(mapValues.size() >= 2) {
				//if(mapValues.get(1) >= 0.25*totalCoverage) {
				//	System.err.println(mapValues.get(1) + " >= 0.25 * " + totalCoverage + " = " + 0.25*totalCoverage);
				//}
				if(mapValues.get(1) >= 2 || (totalCoverage <= 10 && mapValues.get(1) >= 1)) {
					snpList.add(new SNP(pos, alignment));
				}
			}

			/*for(Character c: alleleCounts.keySet()) {
				Integer count = alleleCounts.get(c);
				//System.out.print(c + ":" + count + "  ");
			}
			for(int k = 0; k < mapValues.size(); k++) {
				//System.out.print(mapValues.get(k) + ":");
			}
			if(mapValues.size() > 1 && mapValues.get(1) >= 2){
				//System.out.println("************* POSITION "+ pos);
			}
			else {
				//System.out.println();
			}*/
		}

		return snpList;
	}





}
