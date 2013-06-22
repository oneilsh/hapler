/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package hapler;

import java.util.Comparator;

/**
 *
 * @author soneil
 */
public class SequenceEndPosLengthComparator implements Comparator<Sequence> {


	/**
	 * Returns which sequence ends last. If they end at the same point, returns the shorter one.
	 * If they end at the same point and are the same length, they are "equal"
	 * @param seqi
	 * @param seqj
	 * @return
	 */
	public int compare(Sequence seqi, Sequence seqj) {
		if(seqi.getEndPosition() > seqj.getEndPosition()) {
			return -1;
		}
		else if(seqi.getEndPosition() < seqj.getEndPosition()) {
			return 1;
		}
		else {
			if(seqi.length() < seqj.length()) {
				return -1;
			}
			else if(seqi.length() > seqj.length()) {
				return 1;
			}
			else {
				return 0;
			}
		}
	}



}
