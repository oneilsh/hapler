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
import java.util.HashSet;

/**
 *
 * @author soneil
 */
public class SNP {

	private int position;
	private HashSet<Sequence> sequencesCovered;
	private MultipleAlignment alignment;
	private HashMap<Character, Integer> variantsHash;
	//private boolean masked;




	/*public SNP() {
		sequencesCovered = new ArrayList<Sequence>();
		position = -1;
		//masked = false;
		// Empty Constructor
	}*/
	
	public SNP(int thePosition, MultipleAlignment theAlignment) {
		alignment = theAlignment;
		sequencesCovered = new HashSet<Sequence>();
		position = thePosition;
		variantsHash = theAlignment.variantsHash(thePosition);
		//masked = false;
	}

	/**
	 * Returns true if there are more than three alleles in variantsHash
	 * @return
	 */
	public boolean isBiAllelic() {
		/*for(Character key : variantsHash.keySet()) {
			System.err.print(" " + key + " : " + variantsHash.get(key));
		}
		System.err.println();*/
		if(variantsHash.size() == 2) {
			return true;
		}
		return false;
	}

	/**
	 * Only adds the sequence if it isn't already there.
	 * Throw a hissy fit if the sequence doesn't cover this snp
	 * @param theSequence
	 * @throws Exception
	 */
	public void addCoveredSequence(Sequence theSequence) throws Exception {
		if(!sequencesCovered.contains(theSequence)) {
			if(position <= theSequence.endPosition() && position >= theSequence.getStartPosition()) {
				sequencesCovered.add(theSequence);
			}
			else {
				throw new Exception("Sorry, you are trying to say that the SNP at position" + position
						  + " occurs in the sequence "
						  + theSequence.getName()
						  + " which ranges from "
						  + theSequence.getStartPosition()
						  + " to "
						  + theSequence.endPosition());
			}
		}
	}

	/**
	 * Only tries to remove it if it's already in there
	 * @param theSequence
	 */
	public void removeCoveredSequence(Sequence theSequence) {
		if(sequencesCovered.contains(theSequence)) {
			sequencesCovered.remove(theSequence);
		}
	}

	/**
	 * Returns true if this SNP covers (ie is covered by) the sequence
	 * @param theSequence
	 * @return
	 */
	public boolean coversSequence(Sequence theSequence) {
		return sequencesCovered.contains(theSequence);
	}


	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	/*public boolean isMasked() {
		return masked;
	}

	public void setMasked(boolean masked) {
		this.masked = masked;
	}*/
}
