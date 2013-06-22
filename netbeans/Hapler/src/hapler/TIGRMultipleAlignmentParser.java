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
public class TIGRMultipleAlignmentParser extends AbstractMultipleAlignmentParser {

	public TIGRMultipleAlignmentParser() {
	}

	public ArrayList<MultipleAlignment> openFile(String fileName, String allowGaps) throws Exception {
		ArrayList<MultipleAlignment> multipleAlignmentList = new ArrayList<MultipleAlignment>();
		//Parsing code here...

		//try {
			BufferedReader br = null;
			if(fileName.compareTo("-") != 0) {
				FileInputStream fstream = new FileInputStream(fileName);
				DataInputStream in = new DataInputStream(fstream);
				br = new BufferedReader(new InputStreamReader(in));
			}
			else {
				br = new BufferedReader(new InputStreamReader(System.in));
			}

			// gettingRead will be true if we are reading a read, false if we are
			// reading a consensus
			boolean gettingRead = false;
			MultipleAlignment alignment = null;
			int currentStartPosition = 0;
			String currentName = "";

			// A temporary structure for building up the sequence string
			ArrayList<String> sequenceArray = null;

			String strLine;

			boolean printedGapWarning = false;

			while ((strLine = br.readLine()) != null) {

				// If we are getting a new multiple alignment
				if(strLine.matches("^##.*")) {
					// If this isn't our first contig, we must have a read to stuff into an already
					// existing sequence, which we need to stuff into the old alignment
					if(alignment != null && gettingRead == true) { // our sequence is finished
						boolean hasGaps = this.splitAndAddToAlignment(sequenceArray, currentStartPosition, alignment, allowGaps, currentName);
						if(hasGaps && !printedGapWarning) {
							System.out.println("ALLOWGAPS IS "+ allowGaps+"#!WARNING: alignment " + alignment.getName() + " either has reads containing ~ characters or mate-pair reads. These are being SPLIT into separate reads (e.g., mate pair information is being ignored.)");
							printedGapWarning = true;
						}

					}
					gettingRead = false;
					alignment = createNewMultipleAlignment(strLine);
					multipleAlignmentList.add(alignment);
					//sequence = new Sequence();
					sequenceArray = new ArrayList<String>();
					//System.out.println(strLine);
				}
				// If we are getting a new read
				else if(strLine.matches("^#.*")) {
					// Here we must be done reading some sequence... now, did we just get a read or a consensus
					if(!gettingRead) {
						alignment.setGivenConsensus(join(sequenceArray,""));
					}
					else { // our last sequence is finished
						boolean hasGaps = this.splitAndAddToAlignment(sequenceArray, currentStartPosition, alignment, allowGaps, currentName);
						if(hasGaps && !printedGapWarning) {
							System.out.println("ALLOWGAPS IS "+ allowGaps+"#!WARNING: alignment " + alignment.getName() + " either has reads containing ~ characters or mate-pair reads. These are being SPLIT into separate reads (e.g., mate pair information is being ignored.)");
							printedGapWarning = true;
						}
					}
					gettingRead = true;
					//sequence = createNewSequence(strLine);
					currentName = this.extractNameFromTIGRSeqLine(strLine);
					currentStartPosition = this.extractStartPositionFromTIRGSeqLine(strLine);
					sequenceArray = new ArrayList<String>();
					//System.out.println(strLine + " read def");
				}
				// We are getting a line of either the consensus or the a read
				else {
					sequenceArray.add(strLine);
				}
			}
		//}
		/*catch (Exception e) {
			System.err.println("#############################");
			System.err.println("Error: " + e.getMessage());
			System.err.println("#############################");
			System.err.println();
			e.printStackTrace();
		}*/


		return multipleAlignmentList;
	}


	/**
	 * Given an array of strings which when contatenated represents a sequence, that
	 * may contain ~'s, cats and splits the sequence by ~'s adding them one at a time to the alignment.
	 * @param sequenceArray
	 * @param currentStartPosition
	 * @param alignment
	 * @return
	 */
	private boolean splitAndAddToAlignment(ArrayList<String> sequenceArray, int currentStartPosition, MultipleAlignment alignment, String allowGaps, String currentName) throws Exception {
		boolean hasGaps = false;
		
		HashMap<String, Integer> sequencesToStartPositions = this.splitSequence(sequenceArray, currentStartPosition);

		//If we have gaps but allowGaps is false, choke.
		if(allowGaps.compareTo("false") == 0 && sequencesToStartPositions.size() > 1) {
			throw new Exception("Sorry, the read " + currentName + " is in multiple pieces (e.g. is part of a mate-pair). This isn't allowed when --allow-gaps is false. See --help");
		}

		//If we have gaps, but allow them to split, set the warning and continue
		if(sequencesToStartPositions.size() > 1 && !hasGaps) {
			hasGaps = true;
		}

		Sequence sequence = null;

		// Now we can add each piece as a seperate sequence
		int pieceIndex = 0;
		for(String piece: sequencesToStartPositions.keySet()) {
			int startPos = sequencesToStartPositions.get(piece);
			sequence = new Sequence();
			if(sequencesToStartPositions.size() > 1) {
				sequence.setName(currentName + "_" + pieceIndex);
			}
			else {
				sequence.setName(currentName);
			}

			sequence.addPiece(piece, startPos);
			sequence.setAlignment(alignment);
			pieceIndex = pieceIndex + 1;
			alignment.addSequence(sequence);
		}

		return hasGaps;
	}


	/**
	 * Given a line from a TIGR format assembly file starting with ##
	 * (defining a multiple alignment) create a shell of one that we can
	 * stick reads in
	 * @param inLine
	 * @return
	 */
	private MultipleAlignment createNewMultipleAlignment(String inLine) {
		MultipleAlignment alignment = new MultipleAlignment();
		String[] lineArray = inLine.split("[ \t\n\r]+");
		String name = lineArray[0].replaceFirst("^##", "");
		//int numReads = Integer.parseInt(lineArray[1]);
		//int numBases = Integer.parseInt(lineArray[2]);

		alignment.setName(name);
		return alignment;
	}


	/**
	 * Given a line from a TIGR format assembly file starting with #
	 * (defining a sequence) return the name of the sequence
	 * @param inLine
	 * @return
	 */
	private String extractNameFromTIGRSeqLine(String inLine) {
		//Sequence seq = new Sequence();
		String[] lineArray = inLine.split("[\\(\\)]");
		String name = lineArray[0].replaceFirst("^#","");
		//int startPosition = Integer.parseInt(lineArray[1]);

		//seq.setName(name);
		//seq.setStartPosition(startPosition);
		return name;
	}

	/**
	 * Given a line from a TIGR format assembly file starting with #
	 * (defining a sequence) return the name of the sequence
	 * @param inLine
	 * @return
	 */
	private int extractStartPositionFromTIRGSeqLine(String inLine) {
		//Sequence seq = new Sequence();
		String[] lineArray = inLine.split("[\\(\\)]");
		String name = lineArray[0].replaceFirst("^#","");
		int startPosition = Integer.parseInt(lineArray[1]);

		//seq.setName(name);
		//seq.setStartPosition(startPosition);
		return startPosition;
	}


	/**
	 * Given an arraylist of strings, joins them together to one big string.
	 * @param stringArrayList
	 * @param separator
	 * @return
	 */
	private String join(ArrayList<String> stringArrayList, String separator) {
		StringBuffer sb = new StringBuffer();
		boolean first = true;
		for(String piece : stringArrayList) {
			if(!first) sb.append(separator);
			sb.append(piece);
			first = false;
		}
		return sb.toString();
	}



	/*
	 * Given an arrayList of strings, which when concatted together represents
	 * a single sequence, we need to split it up into pieces by removing segments of ~'s
	 * and add each piece to the seq
	 * returns true if there were any ~ characters that were split on
	 */
	public HashMap<String,Integer> splitSequence(ArrayList<String> seqArray, int startPosition) {
		boolean toRet = false;
		StringBuilder fullSeq = new StringBuilder(join(seqArray,""));
		StringBuilder pieceBuffer = null;
		int pieceIndex = -1;

		HashMap<String,Integer> sequencesToStartPositions = new HashMap<String,Integer>();

		for(int i = 0; i < fullSeq.length(); i++) {
			if(fullSeq.charAt(i) == '~') {
				if(pieceBuffer != null) { // we had a seq that we are ending
				toRet = true;
					// add the piece with the appropriate index and start a new one
					sequencesToStartPositions.put(pieceBuffer.toString(), pieceIndex+startPosition);
					pieceIndex = -1;
					pieceBuffer = null;
				}
				else {
					// do nothing, we're just dropping ~ chars
				}
			}
			else {
				if(pieceBuffer != null) {
					// we're continuing the pieceBuffer
					pieceBuffer.append(fullSeq.charAt(i));
				}
				else {
					// here we're starting a new piece, so we need to record the index
					// and start a new pieceBuffer, not forgetting to add the current char
					pieceBuffer = new StringBuilder();
					pieceIndex = i;
					pieceBuffer.append(fullSeq.charAt(i));
				}
			}
		}
		if(fullSeq.charAt(fullSeq.length()-1) != '~') {
			// if we didnt' end on a ~, we need to add the final crated pieceBuffer
			sequencesToStartPositions.put(pieceBuffer.toString(), pieceIndex+startPosition);
		}
		return sequencesToStartPositions;
	}


}
