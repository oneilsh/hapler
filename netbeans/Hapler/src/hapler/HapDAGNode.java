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
 * This class wraps a (haplotype) sequence and is used in reconstructing contigs
 * from the minimum weight path in the overlap DAG of haplotype sequences.
 * @author soneil
 */
public class HapDAGNode {
	private Haplotype hap;
	private double minWeight;
	private HapDAGNode backPointer;
	private ArrayList<HapDAGNode> connectsTo;
	private double weight;


	public HapDAGNode(Haplotype theHap) {
		hap = theHap;
		connectsTo = new ArrayList<HapDAGNode>();
		weight = 0.0;
		minWeight = Double.POSITIVE_INFINITY;
		backPointer = null;
	}


	public void addConnection(HapDAGNode theNode) {
		connectsTo.add(theNode);
	}

	/**
	 * returns the start position of the sequence wrapped
	 * @return
	 */
	public int startPosition() {
		return hap.startPos();
	}

	/**
	 * returns the end position of the sequence wrapped
	 * @return
	 */
	public int endPosition() {
		return hap.endPos();
	}


	/**
	 * Returns true if the sequence wrapped starts before the other sequence wrapped
	 * @param other
	 * @return
	 */
	public boolean startsBefore(HapDAGNode other) {
		return hap.startsBefore(other.getHap());
	}

	/**
	 * Returns true if the sequence wrapped ends before the other sequence wrapped
	 * @param other
	 * @return
	 */
	public boolean endsBefore(HapDAGNode other) {
		return hap.endsBefore(other.getHap());
	}

	/**
	 * returns true if the sequence wrapped overlaps the other sequence wrapped
	 * @param other
	 * @return
	 */
	public boolean overlaps(HapDAGNode other) {
		return hap.overlaps(other.getHap());
	}


	public double getWeight() {
		return weight;
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}

	public Haplotype getHap() {
		return hap;
	}

	public void setHap(Haplotype seq) {
		this.hap = seq;
	}


	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("%1.4f", weight) + "\t");
		sb.append(startPosition() + "\t");
		sb.append(endPosition() + "\t");
		sb.append("Connections:" + connectsTo.size() + "\t");
		sb.append("\t|\t");
		sb.append(hap.consensusHumanReadable());
		return sb.toString();
	}


	public HapDAGNode getBackPointer() {
		return backPointer;
	}

	public void setBackPointer(HapDAGNode backPointer) {
		this.backPointer = backPointer;
	}

	public double getMinWeight() {
		return minWeight;
	}

	public void setMinWeight(double minWeight) {
		this.minWeight = minWeight;
	}

	public ArrayList<HapDAGNode> getConnectsTo() {
		return connectsTo;
	}

	public ArrayList<HapDAGNode> getMinWeightPathToHere() {
		ArrayList<HapDAGNode> toRet = null;
		if(backPointer != null) {
			toRet = backPointer.getMinWeightPathToHere();
		}
		else {
			toRet = new ArrayList<HapDAGNode>();
		}
		toRet.add(this);
		return toRet;
	}


}
