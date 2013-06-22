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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Parses a fasta file, reading only ids, returning a hashmap from id to sequence.
 * Currently, gaps (~) are left in the sequences and not treated specially.
 * @author soneil
 */
public class FastaParser {

	public FastaParser() {
	}


	public HashMap<String, Sequence> openFile(String fileName) throws Exception {
		HashMap<String,ArrayList<String>> idsToSeqArrays = new HashMap<String,ArrayList<String>>();
		HashMap<String, Sequence> idsToSeqs = new HashMap<String, Sequence>();
		try{
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));


			String currentId = null;

			String strLine = null;
			while ((strLine = br.readLine()) != null) {

				if(strLine.matches("^>.*")) {
					String[] lineArray = strLine.split("[ \t\n\r]+");
					String id = lineArray[0].replaceFirst("^>", "");

					currentId = id;
					if(!idsToSeqArrays.containsKey(id)) {
						idsToSeqArrays.put(id, new ArrayList<String>());
					}
					else {
						throw new Exception("Ooops, we've already seen the id " + id + " in the fasta file " + fileName + ". Erroring out.");
					}
				}
				else {
					idsToSeqArrays.get(currentId).add(strLine);
				}

			}


			boolean printedGapsWarning = false;

			for(String id : idsToSeqArrays.keySet()) {
				ArrayList<String> seqArray = idsToSeqArrays.get(id);
				Sequence seq = new Sequence();
				seq.setName(id);
				boolean seqHasGaps = seq.addPieces(seqArray, 0);
				if(!printedGapsWarning && seqHasGaps) {
					printedGapsWarning = true;
					System.out.println("#!WARNING: Seq id " + id + " in file " + fileName + " has ~ characters in it; this sequence is being split into multiple sequences without ~ characters.");
				}
				idsToSeqs.put(id, seq);
			}


		}
	catch(Exception e){
		System.err.println("Ooops, problem reading file "+ fileName);
		throw e;
		}
	return idsToSeqs;
	}
}
