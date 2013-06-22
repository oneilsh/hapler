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
import java.util.HashSet;

/**
 * A data container for the Reconstructor class to put finished product in.
 * @author soneil
 */
public class ReconstructedConsensus {
	private String consensusAsString;
	private HashSet<Haplotype> uniqueHapsUsed;
	//private HashSet<Haplotype> uniqueLowCoverageHapsUsed;
	private int SNPAlleleCoverage;

	// minSnpCompatibilityList is the array of backpointers that reconstructs the consensus, one per snp
	private ArrayList<SnpDAGNode> minSnpCompatibilityList;
	private ArrayList<SnpDAGNode> crossOverPoints;

	public String getConsensusString() {
		return consensusAsString;
	}
	
	
	
	
	/*
	 * Returns a string of the haplotypes used in the minimum reconstruction,
	 * in order, as well as their average coverages.
	 */
	public String hapsUsedAsString(boolean includeCoverage) throws Exception {
		StringBuilder sb = new StringBuilder();
		if(minSnpCompatibilityList.size() == 0) {
			throw new Exception("I'm attempting to list for you what haplotypes were used to reconstruct a consensus sequence, buuut, there is no haplotype as being associated with any snp. Were any snps called?");
		}
		SnpDAGNode firstNode = minSnpCompatibilityList.get(0);
		sb.append(firstNode.getCompatibleHap().fullName());
		sb.append("(");
		sb.append("0" + ";");
		sb.append( String.format("%.2f", firstNode.getCompatibleHap().averageCoverage()));
		sb.append(")");
		for(int i = 1; i < crossOverPoints.size(); i++) {
			SnpDAGNode nodei = crossOverPoints.get(i);
			sb.append(",");
			sb.append( nodei.getCompatibleHap().fullName());
			sb.append("(");
			if(includeCoverage) {
				sb.append( nodei.getSnp().getPosition() + ";");
				sb.append(String.format("%.2f", nodei.getCompatibleHap().averageCoverage()));
			}
			sb.append(")");
		}

		return sb.toString();
	}
	

	public void setConsensus(String consensus) {
		this.consensusAsString = consensus;
	}

	public ArrayList<SnpDAGNode> getCrossOverPoints() {
		return crossOverPoints;
	}

	public void setCrossOverPoints(ArrayList<SnpDAGNode> crossOverPoints) {
		this.crossOverPoints = crossOverPoints;
	}

	public ArrayList<SnpDAGNode> getMinSnpCompatibilityList() {
		return minSnpCompatibilityList;
	}

	public void setMinSnpCompatibilityList(ArrayList<SnpDAGNode> minSnpCompatibilityList) {
		this.minSnpCompatibilityList = minSnpCompatibilityList;
	}

	public HashSet<Haplotype> getUniqueHapsUsed() {
		return uniqueHapsUsed;
	}

	public void setUniqueHapsUsed(HashSet<Haplotype> uniqueHapsUsed) {
		this.uniqueHapsUsed = uniqueHapsUsed;
	}

	/*public HashSet<Haplotype> getUniqueLowCoverageHapsUsed() {
		return uniqueLowCoverageHapsUsed;
	}

	public void setUniqueLowCoverageHapsUsed(HashSet<Haplotype> uniqueLowCoverageHapsUsed) {
		this.uniqueLowCoverageHapsUsed = uniqueLowCoverageHapsUsed;
	}*/



	public int getSNPAlleleCoverage() {
		return SNPAlleleCoverage;
	}

	public void setSNPAlleleCoverage(int alleleCoverage) {
		this.SNPAlleleCoverage = alleleCoverage;
	}

	public ReconstructedConsensus() {
	}
}
