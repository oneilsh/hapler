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
 *
 * @author soneil
 */
public class SAMMultipleAlignmentParser extends AbstractMultipleAlignmentParser {

	public SAMMultipleAlignmentParser() {
	}
	

	public ArrayList<MultipleAlignment> openFile(String fileName, String allowGaps) throws Exception {
		ArrayList<MultipleAlignment> multipleAlignmentList = new ArrayList<MultipleAlignment>();
		
		
		//Open the input...
		BufferedReader br = null;
		if(fileName.compareTo("-") != 0) {
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream in = new DataInputStream(fstream);
			br = new BufferedReader(new InputStreamReader(in));
		}
		else {
			br = new BufferedReader(new InputStreamReader(System.in));
		}


		HashMap<String, HashMap<String, HashMap<String, Integer>>> scaffoldsToSeqNamesEtc = new HashMap<String, HashMap<String, HashMap<String, Integer>>>();
		MultipleAlignment alignment = null;

		String strLine;
		while ((strLine = br.readLine()) != null) {
			if(strLine.matches("^@.*")) {
				continue;
			}

			String lineArray[] = strLine.split("\\s+");
			String readName = lineArray[0];
			int flag = Integer.parseInt(lineArray[1]);
			String scaffoldName = lineArray[2];
			int startPosition = Integer.parseInt(lineArray[3])-1;
			String mapQ = lineArray[4];
			String cigar = lineArray[5];
			String mateName = lineArray[6];
			int matePosition = Integer.parseInt(lineArray[7]);
			int insertSize = Integer.parseInt(lineArray[8]);
			String seq = lineArray[9];
			String qualSeq = lineArray[10];

			if(!cigar.equals("*")) {

				String gappedSeq = parseCigar(cigar, seq);
				HashMap<String, Integer> gappedSeqPieces = breakByTildes(gappedSeq, startPosition);

				if(!scaffoldsToSeqNamesEtc.containsKey(scaffoldName)) {
					scaffoldsToSeqNamesEtc.put(scaffoldName, new HashMap<String, HashMap<String, Integer>>());
				}
				if(!scaffoldsToSeqNamesEtc.get(scaffoldName).containsKey(readName)) {
					scaffoldsToSeqNamesEtc.get(scaffoldName).put(readName, new HashMap<String, Integer>());
				}

				for(String piece: gappedSeqPieces.keySet()) {
					int startPos = gappedSeqPieces.get(piece);
					scaffoldsToSeqNamesEtc.get(scaffoldName).get(readName).put(piece, startPos);
					//System.out.println("Adding to scaffold " + scaffoldName + " and read " + readName + " the seq " + piece + " starting at " + startPos);
				}
			}
		}

		boolean printedGapsWarning = false;

		for(String scaffold: scaffoldsToSeqNamesEtc.keySet()) {
			alignment = new MultipleAlignment();
			alignment.setName(scaffold);
			multipleAlignmentList.add(alignment);
			HashMap<String, HashMap<String,Integer>> readNamesToPiecesHash = scaffoldsToSeqNamesEtc.get(scaffold);
			for(String readName: readNamesToPiecesHash.keySet()) {
				HashMap<String, Integer> pieceToStartPosHash = readNamesToPiecesHash.get(readName);

				//If we have gaps but allowGaps is false, choke.
				if(allowGaps.compareTo("false") == 0 && pieceToStartPosHash.size() > 1) {
					throw new Exception("Sorry, the read " + readName + " is in multiple pieces (e.g. is part of a mate-pair). This isn't allowed when --allow-gaps is false. See --help");
				}

				//If we have gaps, but allow them to split, throw a warning and continue
				if(pieceToStartPosHash.size() > 1 && !printedGapsWarning) {
					System.out.println("#!WARNING: alignment " + alignment.getName() + " either has reads containing ~ characters or mate-pair reads. These are being SPLIT into separate reads (e.g., mate pair information is being ignored.)");
					printedGapsWarning = true;
				}

				//Add each piece as a seperate sequence
				int pieceIndex = 0;
				for(String piece: pieceToStartPosHash.keySet()) {
					int startPos = pieceToStartPosHash.get(piece);
					Sequence newSeq = new Sequence();
					if(pieceToStartPosHash.size() > 1) {
						newSeq.setName(readName + "_" + pieceIndex);
					}
					else {
						newSeq.setName(readName);
					}
					newSeq.setAlignment(alignment);

					newSeq.addPiece(piece, startPos);

					pieceIndex = pieceIndex + 1;
					alignment.addSequence(newSeq);
				}

			}
		}


		return multipleAlignmentList;
	}

	/**
	 * Given a sequence and a start position for it, if the sequence contains ~ characters,
	 * split the string by those characters, and return a hashmap of each piece to the start position of
	 * the piece. Eg:
	 * AAAA~~~~TT~GGG~~CC, starting at 5, returns
	 * AAAA => 5, TT => 13, GGG => 16, CC => 21
	 * @param seq
	 * @param startPosition
	 * @return
	 */
	private HashMap<String, Integer>  breakByTildes(String seq, int startPosition) {
		HashMap<String, Integer> piecesToStartPositions = new HashMap<String, Integer>();
		StringBuilder sb = new StringBuilder();

		for(int seqIndex = 0; seqIndex < seq.length(); seqIndex++) {
			if(seq.charAt(seqIndex) == '~') {
				if(sb.length() > 0) { //time to stuff our stringbuilder into the hash
					piecesToStartPositions.put(sb.toString(), seqIndex - sb.length() + startPosition);
					sb = new StringBuilder();
				}
			}
			else {
				sb.append(seq.charAt(seqIndex));
			}
		}

		piecesToStartPositions.put(sb.toString(), seq.length() - sb.length() + startPosition);
		return piecesToStartPositions;
	}


	private String parseCigar(String cigar, String seq) throws Exception {
		StringBuilder sb = new StringBuilder();
		String cigarPairs[] = cigar.split("(?<=[MIDNSHP])");
		int seqIndex = 0;

		for(int i = 0; i < cigarPairs.length; i++) {
			Integer num = Integer.parseInt(cigarPairs[i].substring(0, cigarPairs[i].length()-1));
			String letter = cigarPairs[i].substring(cigarPairs[i].length()-1,cigarPairs[i].length());
			//System.out.println(cigarPairs[i] + " => " + num + " : " + letter);
			if(letter.compareTo("S") == 0) {
				seqIndex = seqIndex + num;
			}
			else if(letter.compareTo("H") == 0) {
				// Do nothing
			}
			else if(letter.compareTo("M") == 0) {
				for(int j = 0; j < num; j++) {
					sb.append(seq.charAt(seqIndex));
					seqIndex = seqIndex + 1;
				}
			}
			else if(letter.compareTo("D") == 0 || letter.compareTo("P") == 0) {
				for(int j = 0; j < num; j++) {
					sb.append("-");
				}
			}
			else if(letter.compareTo("I") == 0) {
				throw new Exception("Oops, the cigar string " + cigar + "contains an I character; this tool only deals with SAM formats that are with respect to a gapped reference.");
			}
			else if(letter.compareTo("N") == 0) {
				for(int j = 0; j < num; j++) {
					sb.append("~");
				}
			}
		}

		//System.out.println(cigar + " : " + seq + " => " + sb.toString());
		return sb.toString();
	}

}
