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
import java.util.HashSet;

/**
 * A Sequence is a string of characters in the multiple alignment with a start position and length.
 * It _may_ have gap characters '~'.
 *
 * In the case of paired end reads, the number of consecutive gap characters may be very large.
 *
 * For this reason, internally Sequences are represented as a number of "pieces" each with a start position.
 * These pieces should _not_ overlap or abut each other; ~ characters fill the gaps between them.
 * @author soneil
 */
public class Sequence {
	private HashMap<String, Integer> sequencesToStartPositions;
	private String name;
	private Integer startPosition;
	private HashSet<SNP> snpsCovered;
	private boolean masked;
	// masking is a list of sequences which this sequences _covers_ and doesn't conflict with...
	private ArrayList<Sequence> masking;
	private MultipleAlignment alignment;




	public Sequence() {
		snpsCovered = new HashSet<SNP>();
		startPosition = null;
		sequencesToStartPositions = new HashMap<String, Integer>();
		name = "";
		masked = false;
		masking = new ArrayList<Sequence>();
		alignment = null;
		// Empty Constructor
	}


	public void addMaskedSequence(Sequence seq) {
		masking.add(seq);
	}

	public void removeMaskedSequence(Sequence seq) {
		masking.remove(seq);
	}

	/**
	 * Throw a hissy fit if the SNP is outside the range of this sequence
	 * If a snp already exists at this position, do nothing
	 * @param theSNP
	 * @throws Exception
	 */
	public void addCoveredSNP(SNP theSNP) throws Exception {
		int snpPosition = theSNP.getPosition();
		
		if(!this.containsPosition(snpPosition)) {
			throw new Exception("Sorry, you are trying to add a SNP at position "
					  + theSNP.getPosition()
					  + " in read "
					  + this.getName()
					  + " which spans positions "
					  + this.getStartPosition() + " to " + this.endPosition());
		}
		if(!snpsCovered.contains(theSNP)) {
			snpsCovered.add(theSNP);
		}
	}

	/**
	 * Adds multiple SNPs via addCoveredSNP
	 * @param theSNPs
	 * @throws Exception
	 */
	public void addCoveredSNPs(ArrayList<SNP> theSNPs) throws Exception {
		for(int i = 0; i < theSNPs.size(); i++) {
			addCoveredSNP(theSNPs.get(i));
		}
	}

	/**
	 * Adds a piece to this sequence. This method will throw up if
	 * the piece added overlaps or abuts another piece already contained;
	 * it will also throw up if the piece contains a ~ character (and thus
	 * should be added as two seperate pieces.)
	 * @param seq
	 * @param pieceStartPos
	 */
	public void addPiece(String piece, int pieceStartPos) throws Exception {
		//System.out.println("Adding... piece starting at " + pieceStartPos + " to seq " + this.getName() + ", start position is now " + startPosition);
		if(piece.indexOf("~") != -1) {
			throw new Exception("Sorry, you are trying to add a piece to a sequence starting at position "
					  + pieceStartPos
					  + " which has a ~ character in it. You should attempt to add this piece as two seperate pieces.");
		}
		int thisPieceEnd = pieceStartPos + piece.length() - 1;
		for(String otherPiece: sequencesToStartPositions.keySet()) {
			int otherStart = sequencesToStartPositions.get(otherPiece);
			int otherEnd = otherStart + otherPiece.length() - 1;
			if(otherEnd < pieceStartPos-1 || otherStart > thisPieceEnd + 1) { // ok to add according to this piece
			}
			else { // This piece overlaps or abuts some other piece
				/*throw new Exception("Sorry, you are trying to add a piece to a sequence starting at position "
						+ pieceStartPos
						+ " which either overlaps or abuts another piece which starts at position "
						+ otherStart
						+ ", these two should be added as one piece. (If they don't conflict.)");*/
				System.out.println("#!WARNING: sequence " + this.getName() + " has pieces (e.g. mate-pair ends) which overlap or abut each other. This might break some stuff...");
				//alignment.maskSequence(this);
			}
		}
		//System.out.println("Adding piece to seq " + name + " at position " + pieceStartPos);
		if(startPosition == null) {
			startPosition = pieceStartPos;
		}
		else if(pieceStartPos < startPosition) {
			startPosition = pieceStartPos;
			//System.out.println("Startposition for " + name + " is now " + startPosition);
		}
		sequencesToStartPositions.put(piece, pieceStartPos);
	}

	/**
	 * Returns whether the given SNP is contained as a SNP in the snpsCovered list
	 */
	public boolean containsSNP(SNP theSNP) {
		return snpsCovered.contains(theSNP);
	}



	public String sequenceToString(boolean padded) {
		int myTotalLength = this.length();
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < startPosition; i++) {
			if(padded) {
				sb.append("~");
			}
			else {
				sb.append(" ");
			}
		}
		for(int i = startPosition; i < startPosition+myTotalLength; i++) {
			sb.append("~");
		}
		//sb.append(sequence);
		for(String seq: sequencesToStartPositions.keySet()) {
			int startPos = sequencesToStartPositions.get(seq);
			sb.replace(startPos, startPos+seq.length(), seq);
		}
		return sb.toString();
	}



	public String sequenceToStringReverseCase(boolean padded) {
		String seqString = sequenceToString(padded);
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < seqString.length(); i++) {
			char thisChar = seqString.charAt(i);
			if(Character.toUpperCase(thisChar) == thisChar) {
				thisChar = Character.toLowerCase(thisChar);
			}
			else {
				thisChar = Character.toUpperCase(thisChar);
			}
			sb.append(thisChar);
		}
		
		return sb.toString();
	}


	public String toString() {
		//return name;
		StringBuilder sb = new StringBuilder();
		sb.append(name + "\t");
		sb.append(getStartPosition() + "\t");
		sb.append(endPosition() + "\t");
		ArrayList<SNP> snpsCoveredList = new ArrayList<SNP>(snpsCovered);
		Collections.sort(snpsCoveredList, new SNPPosComparator());
		for(int i = 0; i < snpsCoveredList.size(); i++) {
			SNP snpi = snpsCoveredList.get(i);
			if(i != snpsCoveredList.size()-1) {
				sb.append(snpi.getPosition() + ",");
			}
			else {
				sb.append(snpi.getPosition() + "\t");
			}
		}
		sb.append(sequenceToString(false) + "\t");
		return sb.toString();
	}

	/**
	 * Returns the base at position i in the alignment position i
	 * If the position is outside the read, '~' is returned
	 * @param i
	 * @return
	 */
	public char baseAtAlignmentPosition(int i) {
		for(String seq: sequencesToStartPositions.keySet()) {
			int seqStart = sequencesToStartPositions.get(seq);
			if(i >= seqStart && i < seqStart + seq.length()) {
				return seq.charAt(i-seqStart);
			}
		}
		return '~';
	}

	/**
	 * Returns a substring of the read from alignment position startPos to alignment
	 * position endPos. If this would cover before the read or after the read, ~'s are
	 * used in those places. (Well, whatever baseAtAlignmentPosition sends, actually).
	 * @param startPos
	 * @param endPos
	 * @return
	 */
	public String basesFromAlignmentPositionTo(int startPos, int endPos) {
		StringBuffer sb = new StringBuffer();
		for(int i = startPos; i <= endPos; i++) {
			sb.append(baseAtAlignmentPosition(i));
		}
		return sb.toString();
	}

	
	/**
	 * Returns whether this sequence contains position pos as a non-gap character in the multiple alignment
	 * @param pos
	 * @return
	 */
	public boolean containsPosition(int pos) {
		for(String seq: sequencesToStartPositions.keySet()) {
			int seqStart = sequencesToStartPositions.get(seq);
			if(pos >= seqStart && pos < seqStart + seq.length()) {
				return true;
			}
		}
		return false;
	}



	/**
	 * DEPRECATED
	 * Returns true if the other sequence is contained within this one
	 * @param other
	 * @return
	 */
	/*public boolean covers(Sequence other) {
		boolean covers = false;
		if(this.getStartPosition() <= other.getStartPosition() && this.endPosition() >= other.endPosition()) {
			covers = true;
		}
		return covers;
	}*/

	
	/**
	 * Returns true if any character of the sequence is '~'
	 * @return
	 */
	public boolean hasGaps() {
		if(sequencesToStartPositions.size() > 1) {
			return true;
		}
		return false;
	}


	/**
	 * Returns true if this sequence starts before the other sequence in the multiple alignment.
	 * @param other
	 * @return
	 */
	public boolean startsBefore(Sequence other) {
		if(this.getStartPosition() < other.getStartPosition()) {
			return true;
		}
		return false;
	}


	/**
	 * Returns true if this sequence agrees with the other at all the unmasked SNPs that this
	 * sequence covers.
	 * @param other
	 * @return
	 */
	public boolean agreesWithAtCoveredSNPs(Sequence other) {
		for(SNP theSnp : snpsCovered) {
		//for(int i = 0; i < snpsCovered.size(); i++) {
			//SNP theSnp = snpsCovered.get(i);
			//if(!theSnp.isMasked()) {
				int pos = theSnp.getPosition();
				if(this.baseAtAlignmentPosition(pos) != '~' && other.baseAtAlignmentPosition(pos) != '~') {
					if(this.baseAtAlignmentPosition(pos) != other.baseAtAlignmentPosition(pos)) {
						return false;
					}
				}
			//}
		}
		return true;
	}


	/**
	 * Returns true if this sequence overlaps the other at all the unmasked SNPs that this
	 * sequence covers.
	 * @param other
	 * @return
	 */
	public boolean overlapsAtCoveredSNPs(Sequence other) {
		for(SNP theSnp : snpsCovered) {
		//for(int i = 0; i < snpsCovered.size(); i++) {
			//SNP theSnp = snpsCovered.get(i);
			//if(!theSnp.isMasked()) {
				int pos = theSnp.getPosition();
				if(other.baseAtAlignmentPosition(pos) == '~') {
					return false;
				}
			//}
		}
		return true;
	}

	/**
	 * returns true if this sequence overlaps the other, as determined by start and end positions
	 * of both
	 * @param other
	 * @return
	 */
	public boolean overlaps(Sequence other) {
		boolean toRet = true;
		if(this.getStartPosition() > other.endPosition() || other.getStartPosition() > this.endPosition()) {
			toRet = false;
		}
		return toRet;
	}

	/**
	 * Returns the first SNP, in sorted order
	 * If there are no SNPs, returns null
	 * @return
	 */
	public SNP firstSNP() {
		SNP toRet = null;
		int firstSnpPosition = Integer.MAX_VALUE;
		for(SNP snpi : snpsCovered) {
			if(firstSnpPosition < snpi.getPosition()) {
				firstSnpPosition = snpi.getPosition();
				toRet = snpi;
			}
		}
		//Collections.sort(snpsCovered, new SNPPosComparator());
		//if(snpsCovered.size() > 0) toRet = snpsCovered.get(0);
		return toRet;
	}

	/**
	 * DEPRECATED
	 * Returns whether two sequences share a column in the multiple alignment, ie, if there are any positions
	 * for which both are not ~.
	 * @param other
	 * @return
	 */
	/*public boolean overlaps(Sequence other) {
		boolean overlaps = true;
		if(this.endPosition() < other.getStartPosition() || this.getStartPosition() > other.endPosition()) {
			overlaps = false;
		}

		return overlaps;
	}*/

	/**
	 * If the two sequences overlap, and conflict* at some SNP position (of those
	 * defined in this sequence), then return true, else return false.
	 *
	 * * Two sequences conflict at a position i if seq1[i] != seq2[i] AND neither of these are ~
	 * @param other
	 * @return
	 */
	public boolean conflictsWithAtAnySNP(Sequence other) {
		for(SNP snp : snpsCovered) {
		//for(int i = 0; i < snpsCovered.size(); i++) {
			//SNP snp = snpsCovered.get(i);
			int pos = snp.getPosition();
			if(this.baseAtAlignmentPosition(pos) != '~' && other.baseAtAlignmentPosition(pos) != '~') {
				if(this.baseAtAlignmentPosition(pos) != other.baseAtAlignmentPosition(pos)) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Returns true if this sequence conflicts with the other at any snp in the given SNP list.
	 * @param other
	 * @param snpList
	 * @return
	 */
	public boolean conflictsWithAtAnySNP(Sequence other, ArrayList<SNP> snpList) {
		for(int i = 0; i < snpList.size(); i++) {
			SNP snp = snpList.get(i);
			int pos = snp.getPosition();
			if(this.baseAtAlignmentPosition(pos) != '~' && other.baseAtAlignmentPosition(pos) != '~') {
				if(this.baseAtAlignmentPosition(pos) != other.baseAtAlignmentPosition(pos)) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Returns true if an unmasked SNP is contained in both sequences
	 * @param other
	 * @return
	 */
	public boolean sharesSNPWith(Sequence other) {
		for(SNP snp : snpsCovered) {
		//for(int i = 0; i < snpsCovered.size(); i++) {
			//SNP snp = snpsCovered.get(i);
			if(other.containsSNP(snp)) {
				return true;
			}
		}

		return false;
	}

		/**
	 * If the two sequences overlap, and conflict* at some SNP position (of those
	 * defined in this sequence), then return true, else return false.
	 *
	 * * Two sequences conflict at a position i if seq1[i] != seq2[i] AND neither of these are ~
	 * @param other
	 * @return
	 */
	public boolean conflictsWithAtSNP(Sequence other) {
			for(SNP snp : snpsCovered) {
			//for(int i = 0; i < snpsCovered.size(); i++) {
				//SNP snp = snpsCovered.get(i);
				int pos = snp.getPosition();
				if(this.baseAtAlignmentPosition(pos) != '~' && other.baseAtAlignmentPosition(pos) != '~') {
					if(this.baseAtAlignmentPosition(pos) != other.baseAtAlignmentPosition(pos)) {
						//if(!snp.isMasked()) {
							return true;
						//}
					}
				}
			}
		return false;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}


	public HashSet<SNP> getSNPsCovered() {
		return snpsCovered;
	}

	/**
	 * DEPRECATED
	 * @return
	 */
	/*public void setSequence(String sequence) {
		this.sequence = sequence;
	}*/


	public int getStartPosition() {
		return startPosition;
	}

	public int getEndPosition() {
		return startPosition + length() - 1;
	}

	public int length() {
		//return sequence.length();
		Integer firstStart = null;
		Integer lastEnd = null;
		for(String seq: sequencesToStartPositions.keySet()) {
			int seqStart = sequencesToStartPositions.get(seq);
			if(firstStart == null) {
				firstStart = seqStart;
			}
			else if(seqStart < firstStart) {
				firstStart = seqStart;
			}
			int seqEnd = seqStart + seq.length() - 1;
			if(lastEnd == null) {
				lastEnd = seqEnd;
			}
			else if(seqEnd > lastEnd) {
				lastEnd = seqEnd;
			}
		}
		return lastEnd - firstStart + 1;
	}


	/**
	 * Returns the number of bases (non-~ characters) this sequence has
	 * @return
	 */
	public int numBases() {
		int count = 0;
		for(String seq: sequencesToStartPositions.keySet()) {
			count = count + seq.length();
		}
		return count;
	}




	public int endPosition() {
		Integer lastEnd = null;
		//System.out.println("Number of peices in  " + this.name + " : " + sequencesToStartPositions.size());
		for(String seq: sequencesToStartPositions.keySet()) {
			int seqStart = sequencesToStartPositions.get(seq);

			int seqEnd = seqStart + seq.length() - 1;
			if(lastEnd == null) {
				lastEnd = seqEnd;
			}
			else if(seqEnd > lastEnd) {
				lastEnd = seqEnd;
			}
		}
		return lastEnd;
	}
	

	/**
	 * Returns true if this sequence should mask the other: ie, if the other's SNPs covered loci
	 * are a subset of this one, and they agree at all those alleles.
	 * @param other
	 * @return
	 */
	public boolean shouldMask(Sequence other) throws Exception {
		if(other == this) return false; // a sequence shouldn't mask itself...
		boolean toRet = true;
		HashSet<SNP> otherSNPs = other.getSNPsCovered();
		for(SNP snp : otherSNPs) {
		//for(int i = 0; i < otherSNPs.size(); i++){
		//	SNP snp = otherSNPs.get(i);
			int pos = snp.getPosition();
			if(other.baseAtAlignmentPosition(pos) != '~') {
				if(!this.containsPosition(pos)){
					toRet = false;
				}
				else if(this.baseAtAlignmentPosition(pos) != other.baseAtAlignmentPosition(pos)){
					toRet = false;
				}
			}
			else {
				// This shouldn't ever happen... if the sequence covers a SNP then it should not be a ~ at that position
				throw new Exception("Sequence.java: shouldMask: this shouldn't happen... Seq " + other.getName() + " at position " + pos + " is " + other.baseAtAlignmentPosition(pos));
			}
		}

		
		return toRet;
	}





	/*
	 * Given an arrayList of strings, which when concatted together represents
	 * a single sequence, we need to split it up into pieces by removing segments of ~'s
	 * and add each piece to the seq
	 * returns true if there were any ~ characters that were split on
	 */
	public boolean addPieces(ArrayList<String> seqArray, int startPosition) throws Exception {
		boolean toRet = false;
		StringBuilder fullSeq = new StringBuilder(join(seqArray,""));
		StringBuilder pieceBuffer = null;
		int pieceIndex = -1;
		for(int i = 0; i < fullSeq.length(); i++) {
			if(fullSeq.charAt(i) == '~') {
				if(pieceBuffer != null) { // we had a seq that we are ending
				toRet = true;
					// add the piece with the appropriate index and start a new one
					this.addPiece(pieceBuffer.toString(), pieceIndex+startPosition);
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
			this.addPiece(pieceBuffer.toString(), pieceIndex+startPosition);
		}
		return toRet;
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





	
	public boolean isMasked() {
		return masked;
	}

	/**
	 * returns true if the status was changed
	 * @param masked
	 */
	public boolean setMasked(boolean masked) {
		boolean toRet = false;
		if(masked != this.masked) toRet = true;
		this.masked = masked;
		return toRet;
	}

	public ArrayList<Sequence> getMasking() {
		return masking;
	}

	public MultipleAlignment getAlignment() {
		return alignment;
	}

	public void setAlignment(MultipleAlignment alignment) {
		this.alignment = alignment;
	}

}
