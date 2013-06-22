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

/**
 *
 * @author soneil
 */
public class DisjointSetNode {

	private DisjointSetNode parent;
	private int rank;
	private Object data;


	public DisjointSetNode() {
		rank = 0;
		parent = this;
		data = null;
	}

	public DisjointSetNode(Object theData) {
		rank = 0;
		parent = this;
		data = theData;
	}

	public void union(DisjointSetNode other) {
		DisjointSetNode myRoot = this.find();
		DisjointSetNode otherRoot = other.find();
		if(myRoot.getRank() > otherRoot.getRank()) {
			otherRoot.setParent(myRoot);
		}
		else if(myRoot.getRank() < otherRoot.getRank()) {
			myRoot.setParent(otherRoot);
		}
		else {
			otherRoot.setParent(myRoot);
			myRoot.setRank(myRoot.getRank()+1);
		}
	}


	public DisjointSetNode find() {
		if(parent == this) {
			return this;
		}
		else {
			parent = parent.find();
			return parent;
		}
	}

	public Object getData() {
		return data;
	}

	public void setData(Object data) {
		this.data = data;
	}


	public int getRank() {
		return rank;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}

	public DisjointSetNode getParent() {
		return parent;
	}

	public void setParent(DisjointSetNode parent) {
		this.parent = parent;
	}
}
