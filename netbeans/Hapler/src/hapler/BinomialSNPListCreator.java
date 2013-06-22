/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package hapler;

import java.util.ArrayList;
import JSci.maths.statistics.*;
import java.util.Collections;
import java.util.HashMap;

/**
 *
 * @author soneil
 */
public class BinomialSNPListCreator extends AbstractSNPListCreator {
	HashMap<Integer, BinomialDistribution> binomialDists;
	double errorRate;
	double alpha;


	public BinomialSNPListCreator(double e, double a) {
		// Empty constructor
		binomialDists = new HashMap<Integer, BinomialDistribution>();
		errorRate = e;
		alpha = a;
	}


	public ArrayList<SNP> computeSNPList(MultipleAlignment alignment) {
		ArrayList<SNP> snpList = new ArrayList<SNP>();
		double correctedAlpha = alpha/alignment.getLength();
		//System.err.println("corrected alpha is: " + correctedAlpha);
		for(int pos = 0; pos < alignment.getLength(); pos++) {
			HashMap<Character, Integer> variantsHash = alignment.variantsHash(pos);

			if(variantsHash.size() >= 2) {
				ArrayList<Integer> alleleCounts = new ArrayList<Integer>();
				alleleCounts.addAll(variantsHash.values());
				Collections.sort(alleleCounts);
				Collections.reverse(alleleCounts);
				//System.err.println("AlleleCount(0): " + alleleCounts.get(0) + ", AlleleCount(1): " + alleleCounts.get(1));

				int majorityCount = alleleCounts.get(0);
				int minorityCount = alleleCounts.get(1);
				int depth = 0;
				for(int i = 0; i < alleleCounts.size(); i++) {depth = depth + alleleCounts.get(i);}

				double correctedErrorRate = errorRate/(alignment.getSymbolDictionary().size() - 1.0);

				if(!binomialDists.containsKey(depth)) binomialDists.put(depth, new BinomialDistribution(depth, correctedErrorRate));
				double inverse = binomialDists.get(depth).inverse(1-correctedAlpha);
				//System.err.println(" Inverse: " + inverse);
				if(minorityCount > inverse) {
					//System.err.println(" MinorityCount: " + minorityCount + " > " + inverse + " ... Calling a SNP!");
					snpList.add(new SNP(pos, alignment));
				}
			}


			// compute the coverage of this position, as well as the coverage of the majority allele
			/*int depth = 0;
			int maxAlleleCount = -1;
			for(Character key : variantsHash.keySet()) {
				int count = variantsHash.get(key);
				//if(variantsHash.size() >= 2) System.err.print(" " + pos + "::: " + key + " : " + count);
				depth = depth + count;
				if(count > maxAlleleCount) maxAlleleCount = count;
			}
			//if(variantsHash.size() >= 2) System.err.println();
			if(depth >= 2) {
				int numDiffs = depth - maxAlleleCount;
				if(!binomialDists.containsKey(depth)) binomialDists.put(depth, new BinomialDistribution(depth, errorRate));

				double inverse = binomialDists.get(depth).inverse(1-correctedAlpha);
				if(numDiffs > 0) {
					//System.err.println("Pos is " + pos + ", Depth is " + depth + ", numDiffs is " + numDiffs + ", e is " + errorRate + ", inverse CDF binomial is " + inverse);
					if(numDiffs > inverse) {
						//System.err.println("  Stat sig!");
						snpList.add(new SNP(pos, alignment));
					}
				}
			}*/
		}

		return snpList;
	}
}
