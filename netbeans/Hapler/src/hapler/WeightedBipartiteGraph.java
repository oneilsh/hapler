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
import java.util.Set;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.jgrapht.graph.SimpleWeightedGraph;

/**
 * Given a List<Sequence> which appear as a connected component in the conflict graph,
 * build a bipartite representation of the directed, transitive, compatibility graph.
 *
 * After the bipartite representation is constructed, we can run the weighted matching algorithm.
 *
 * Note that we must sort the sequences by start position; then, any sequence X which
 * starts before sequence Y indicates a directed edge from X->Y in the compatibility graph.
 * If the sequences are all gapless, then the compatibility graph will be transitive--allowing for
 * clique cover (solving min color on the conflict compliment) to be solved by path cover, which
 * we do in the bipartite representation using the matchine.
 * @author soneil
 */
public class WeightedBipartiteGraph {

	private SimpleWeightedGraph<GraphNode, DefaultWeightedEdge> graph;
	private ArrayList<GraphNode> leftNodes;
	private ArrayList<GraphNode> rightNodes;

	public SimpleWeightedGraph<GraphNode, DefaultWeightedEdge> getGraph() {
		return graph;
	}

	public void setGraph(SimpleWeightedGraph<GraphNode, DefaultWeightedEdge> graph) {
		this.graph = graph;
	}

	public ArrayList<GraphNode> getLeftNodes() {
		return leftNodes;
	}

	public void setLeftNodes(ArrayList<GraphNode> leftNodes) {
		this.leftNodes = leftNodes;
	}

	public ArrayList<GraphNode> getRightNodes() {
		return rightNodes;
	}

	public void setRightNodes(ArrayList<GraphNode> rightNodes) {
		this.rightNodes = rightNodes;
	}




	/**
	 * A representation of a weighted bipartite graph
	 * @param sequences
	 */

	public WeightedBipartiteGraph(SimpleWeightedGraph<GraphNode, DefaultWeightedEdge> theGraph, ArrayList<GraphNode> theLeftNodes, ArrayList<GraphNode> theRightNodes) {
		graph = theGraph;
		leftNodes = theLeftNodes;
		rightNodes = theRightNodes;
	}

	/**
	 * Creates a weighted bipartite graph which represents the input DAG. Each node in the 
	 * dag is present on both the left and the right, and a connection from A to B in the dag 
	 * is shown as an edge from A left to B right in the bipartite graph.
	 * 
	 * If smallEdgesLeftToRight is used, A left to A right are connected with a small edge (1/(n+1)), for each pair A left and A right
	 * TODO: This is currently O(n^2), could probably be made faster
	 * @param dag
	 * @param smallEdgesLeftToRight
	 */
	public WeightedBipartiteGraph(SimpleDirectedGraph<GraphNode, DefaultEdge> dag, boolean smallEdgesLeftToRight) {
		Set<GraphNode> dagNodes = dag.vertexSet();
		SimpleWeightedGraph<GraphNode, DefaultWeightedEdge> bpGraph = new SimpleWeightedGraph<GraphNode, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		ArrayList<GraphNode> newLeftNodes = new ArrayList<GraphNode>();
		ArrayList<GraphNode> newRightNodes = new ArrayList<GraphNode>();
		HashMap<GraphNode, GraphNode> dagNodesToLeftNodes = new HashMap<GraphNode, GraphNode>();
		HashMap<GraphNode, GraphNode> dagNodesToRightNodes = new HashMap<GraphNode, GraphNode>();

		for(GraphNode dagNode: dagNodes) {
			GraphNode leftNode = new GraphNode(dagNode.getData());
			GraphNode rightNode = new GraphNode(dagNode.getData());
			newLeftNodes.add(leftNode);
			newRightNodes.add(rightNode);
			dagNodesToLeftNodes.put(dagNode, leftNode);
			dagNodesToRightNodes.put(dagNode, rightNode);
			bpGraph.addVertex(leftNode);
			bpGraph.addVertex(rightNode);
		}

		for(GraphNode dagNodej: dagNodes) {
			for(GraphNode dagNodek: dagNodes) {
				if(dagNodej != dagNodek) {
					// TODO: containsEdge is a memory hotspot (creates a lot of iterators...)
					if(dag.containsEdge(dagNodej, dagNodek)) {
						//System.out.println("Dag contains an edge..");
						bpGraph.addEdge(dagNodesToLeftNodes.get(dagNodej), dagNodesToRightNodes.get(dagNodek));
						bpGraph.setEdgeWeight(bpGraph.getEdge(dagNodesToLeftNodes.get(dagNodej), dagNodesToRightNodes.get(dagNodek)), 1.0);
					}
				}
			}
		}

		if(smallEdgesLeftToRight) {
			int numNodes = dagNodes.size();
			double epsilon = 1.0/(numNodes + 1.0);
			for(GraphNode dagNode: dagNodes) {
				bpGraph.addEdge(dagNodesToLeftNodes.get(dagNode), dagNodesToRightNodes.get(dagNode));
				bpGraph.setEdgeWeight(bpGraph.getEdge(dagNodesToLeftNodes.get(dagNode), dagNodesToRightNodes.get(dagNode)), epsilon);
			}
		}

		graph = bpGraph;
		leftNodes = newLeftNodes;
		rightNodes = newRightNodes;
	}


	/**
	 * Shuffles the order of the left and right node sets, such that they end up in the same relative order
	 * to each other. (EG, ABCD left and ABCD right to ACDB left ACDB right)
	 */
	public void shuffleOrder() {
		HashMap<GraphNode, GraphNode> leftToRightNodes = new HashMap<GraphNode, GraphNode>();

		for(int i = 0; i < leftNodes.size(); i++) {
			GraphNode leftNode = leftNodes.get(i);
			GraphNode rightNode = rightNodes.get(i);
			leftToRightNodes.put(leftNode, rightNode);
		}

		Collections.shuffle(leftNodes);
		rightNodes = new ArrayList<GraphNode>();
		for(int i = 0; i < leftNodes.size(); i++) {
			GraphNode leftNode = leftNodes.get(i);
			rightNodes.add(leftToRightNodes.get(leftNode));
		}

	}


		/**
	 * The bipartitegraph should contain as nodes leftNodes and rightNodes,
	 * where |leftNodes| = |rightNodes| = n
	 * O(n^2) memory, O(n^3) time?
	 * @param bipartiteGraph
	 * @param leftNodes
	 * @param rightNodes
	 * @return
	 */
	public WeightedBipartiteGraph maxWeightedBipartiteMatching() {

		// weightMatrix[leftNodeIndex][rightNodeIndex]
		float[][] weightMatrix = buildWeightMatrix();

		// Find the maximum of the weight matrix that we built
		float maxvalue = weightMatrix[0][0];
		for(int i = 0; i < weightMatrix.length; i++) {
			for(int j = 0; j < weightMatrix[0].length; j++) {
				if(weightMatrix[i][j] > maxvalue) {
					maxvalue = weightMatrix[i][j];
				}
			}
		}

		// Adjust it so we are doing max weighted matching
		for(int i = 0; i < weightMatrix.length; i++) {
			for(int j = 0; j < weightMatrix[0].length; j++) {
				weightMatrix[i][j] = maxvalue - weightMatrix[i][j];
			}
		}

		// assignments[leftnodeIndex][rightnodeIndex]
		int[][] assignments = new Hungarian(weightMatrix).execute();


		// Make a new bipartite graph
		ArrayList<GraphNode> leftNodes = this.getLeftNodes();
		ArrayList<GraphNode> rightNodes = this.getRightNodes();
		SimpleWeightedGraph<GraphNode, DefaultWeightedEdge> newGraph = new SimpleWeightedGraph<GraphNode, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		//for(GraphNode lnode: leftNodes) {
		for(int i = 0; i < leftNodes.size(); i++) {
			GraphNode lnode = leftNodes.get(i);
			newGraph.addVertex(lnode);
		}
		//for(GraphNode rnode: rightNodes) {
		for(int i = 0; i < rightNodes.size(); i++) {
			GraphNode rnode = rightNodes.get(i);
			newGraph.addVertex(rnode);
		}

		for(int i = 0; i < assignments.length; i++) {
			int leftNodeIndex = assignments[i][0];
			int rightNodeIndex = assignments[i][1];

			GraphNode leftNode = leftNodes.get(leftNodeIndex);
			GraphNode rightNode = rightNodes.get(rightNodeIndex);
			// we have an assignment from leftNode to rightNode
			// If the original graph actually had that edge, then put it in the new graph
			if(this.getGraph().getEdge(leftNode, rightNode) != null) {
			double oldWeight = this.getGraph().getEdgeWeight(this.getGraph().getEdge(leftNode, rightNode));
				newGraph.addEdge(leftNode, rightNode);
				newGraph.setEdgeWeight(newGraph.getEdge(leftNode, rightNode), oldWeight);
				Sequence leftSeq = (Sequence)leftNode.getData();
				Sequence rightSeq = (Sequence)rightNode.getData();
			}
		}

		return new WeightedBipartiteGraph(newGraph, leftNodes, rightNodes);
	}

	/**
	 * Given a weighted bipartite graph defined by leftNodes and rightNodes, returns
	 * a weight matrix representation of the graph, with leftNodes in the first index
	 * and rightNodes in the second index, eg, if
	 * leftNode 5 connects to rightNode 3 (in their given orders) with weight 0.5,
	 * then weightMatrix[5][3] = 0.5.
	 * @param bipartiteGraph
	 * @param leftNodes
	 * @param rightNodes
	 * @return
	 */
	private float[][] buildWeightMatrix() {
		// initialize weight matrix
		float[][] weightMatrix = new float[leftNodes.size()][rightNodes.size()];
		for(int i = 0; i < leftNodes.size(); i++) {
			for(int j = 0; j < rightNodes.size(); j++) {
				weightMatrix[i][j] = 0.0f;
			}
		}

		// Let's make a hash set of the leftnodes for fast querying against
		HashSet<GraphNode> leftNodesHashSet = new HashSet<GraphNode>();
		for(int i = 0; i < leftNodes.size(); i++) {
			GraphNode nodei = leftNodes.get(i);
			leftNodesHashSet.add(nodei);
		}

		// Build the weight matrix from
		Set<DefaultWeightedEdge> edgeSet = graph.edgeSet();
		for(DefaultWeightedEdge edge: edgeSet) {
			//  I can't gaurantee which "side" of each edge the source is
			GraphNode nodeA = graph.getEdgeSource(edge);
			GraphNode nodeB = graph.getEdgeTarget(edge);
			float thisWeightFloat = new Float(graph.getEdgeWeight(edge));
			//System.out.println(thisWeightFloat);
			if(leftNodesHashSet.contains(nodeA)) {
				// nodeA is a leftNode
				weightMatrix[leftNodes.indexOf(nodeA)][rightNodes.indexOf(nodeB)] = thisWeightFloat;
			}
			else{
				// nodeB is a leftNode
				weightMatrix[leftNodes.indexOf(nodeB)][rightNodes.indexOf(nodeA)] = thisWeightFloat;
			}
		}
		return weightMatrix;
	}


}
