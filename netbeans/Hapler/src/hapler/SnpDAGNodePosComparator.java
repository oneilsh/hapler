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

import java.util.Comparator;

/**
 *
 * @author soneil
 */
public class SnpDAGNodePosComparator implements Comparator<SnpDAGNode> {
	public int compare(SnpDAGNode nodei, SnpDAGNode nodej) {
		if(nodei.getPosition() < nodej.getPosition()) {
			return -1;
		}
		else if(nodej.getPosition() < nodei.getPosition()) {
			return 1;
		}
		else {
			return 0;
		}
	}
}
