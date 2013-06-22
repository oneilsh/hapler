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
 * This class is used to associate an allele at a SNP with a list of haplotypes it is
 * compatible with. Useful for dynamic programming min-chimeric-cover of a consensusHumanReadable.
 * @author soneil
 */
public class SnpDAGNode {
	private Haplotype compatibleHap;
	private SNP snp;
	private int minChimerisms; // minimum number of chimerisms to go through (between snp sites) to get here
	private int minLowCount; // minimum number of low-coverage haps to go through to get here (within the constraint of minChimerisms)
	private HashSet<Haplotype> minUsage; // minimum set of haps used to get here (within above two constraints)
	private SnpDAGNode backPointer;
	private int maxCoverage; // maximum haplotype specific coverage used up until this point


	public SnpDAGNode(SNP theSNP, Haplotype theCompatibleHap) {
		snp = theSNP;
		minChimerisms = Integer.MAX_VALUE;
		minLowCount = Integer.MAX_VALUE;
		minUsage = new HashSet<Haplotype>();
		minUsage.add(compatibleHap);
		compatibleHap = theCompatibleHap;
		backPointer = null;
		maxCoverage = theCompatibleHap.coverageAtPosition(theSNP.getPosition());
	}



	public boolean startsBefore(SnpDAGNode other) {
		if(this.getPosition() < other.getPosition()) {
			return true;
		}
		else {
			return false;
		}
	}



	public int getPosition() {
		return snp.getPosition();
	}

	public int getMinChimerisms() {
		return minChimerisms;
	}

	public void setMinChimerisms(int minChim) {
		this.minChimerisms = minChim;
	}



	public Haplotype getCompatibleHap() {
		return compatibleHap;
	}

	public void setCompatibleHap(Haplotype compatibleHap) {
		this.compatibleHap = compatibleHap;
	}

	public SnpDAGNode getBackPointer() {
		return backPointer;
	}

	public void setBackPointer(SnpDAGNode backPointer) {
		this.backPointer = backPointer;
	}

	public SNP getSnp() {
		return snp;
	}

	public void setSnp(SNP snp) {
		this.snp = snp;
	}

	public int getMinLowCount() {
		return minLowCount;
	}

	public void setMinLowCount(int minLowCount) {
		this.minLowCount = minLowCount;
	}

	public HashSet<Haplotype> getMinUsage() {
		return minUsage;
	}

	public void setMinUsage(HashSet<Haplotype> minUsage) {
		this.minUsage = minUsage;
	}

	public ArrayList<SnpDAGNode> backPointersList() {
		if(backPointer == null) {
			ArrayList<SnpDAGNode> toRet = new ArrayList<SnpDAGNode>();
			toRet.add(this);
			return toRet;
		}
		else {
			ArrayList<SnpDAGNode> toRet = backPointer.backPointersList();
			toRet.add(this);
			return toRet;
		}
	}

	public int getMaxCoverage() {
		return maxCoverage;
	}

	public void setMaxCoverage(int maxCoverage) {
		this.maxCoverage = maxCoverage;
	}
}
