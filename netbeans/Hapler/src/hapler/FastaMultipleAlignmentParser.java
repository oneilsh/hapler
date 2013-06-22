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
 * Parses a fasta file, assuming that each sequence is a "read" starting at position 0.
 * Useful for evaluating real "ground truth" haplotypes...
 * Does not set an id for the alignment...
 * @author soneil
 */
public class FastaMultipleAlignmentParser  extends AbstractMultipleAlignmentParser {

	public FastaMultipleAlignmentParser( ) {
	}

	public ArrayList<MultipleAlignment> openFile(String fileName, String allowGaps) throws Exception {
		ArrayList<MultipleAlignment> alignmentList = new ArrayList<MultipleAlignment>();
		
		FastaParser parser = new FastaParser();
		HashMap<String, Sequence> idsToSeqs = parser.openFile(fileName);

		MultipleAlignment alignment = new MultipleAlignment();


		for(String id : idsToSeqs.keySet()) {
			Sequence seq = idsToSeqs.get(id);
			alignment.addSequence(seq);
		}

		alignmentList.add(alignment);

		return alignmentList;
	}


}
