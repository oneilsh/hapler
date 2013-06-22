/*
 *  Copyright 2010 Shawn Thomas O'Neil
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
import java.util.List;
import java.util.Set;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author soneil
 */
public class MultipleAlignmentTest {

	private MultipleAlignment alignment;
	private Sequence seq1;
	private Sequence seq2;
	private Sequence seq3;
	private Sequence seq4;
	private SNP newSNP;


	// 0 1 2 3 4 5 6 7 8 9
	// A T G C
	//     G C - T
	//         A T ~ T G A
	//           T T T
   public MultipleAlignmentTest() {
		seq1 = new Sequence();
		seq2 = new Sequence();
		seq3 = new Sequence();
		seq4 = new Sequence();
		seq1.setSequence("ATGC"); seq1.setName("seq1"); seq1.setStartPosition(0);
		seq2.setSequence("GC-T"); seq2.setName("seq2"); seq2.setStartPosition(2);
		seq3.setSequence("AT~TGA"); seq3.setName("seq3"); seq3.setStartPosition(4);
		seq4.setSequence("TTT"); seq4.setName("seq4"); seq4.setStartPosition(5);
		alignment = new MultipleAlignment();
		alignment.addSequence(seq1);
		alignment.addSequence(seq2);
		alignment.addSequence(seq3);
		alignment.addSequence(seq4);
		newSNP = new SNP(4);
		alignment.addSNP(newSNP);
    }

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }



	/**
	 * Test of unmaskedColumn method, of class MultipleAlignment.
	 */
	@Test
	public void testGetColumn() {
		System.out.println("* TEST: getColumn");
		char[] col0 = alignment.unmaskedColumn(0);
		char[] col3 = alignment.unmaskedColumn(3);
		char[] col6 = alignment.unmaskedColumn(6);
		char[] col9 = alignment.unmaskedColumn(9);
		assertEquals(4, col0.length);
		assertEquals(4, col3.length);
		assertEquals(4, col6.length);
		assertEquals(4, col9.length);
		int numAs = 0;
		int numCs = 0;
		int numTs = 0;
		int numGaps = 0;
		for(char c : col0) {
			if(c == 'A') numAs = numAs + 1;
			else if(c == 'C') numCs = numCs + 1;
			else if(c == '~') numGaps = numGaps + 1;
			else if(c == 'T') numTs = numTs + 1;
		}
		assertEquals(1, numAs);
		assertEquals(3, numGaps);

		numAs = 0;
		numCs = 0;
		numTs = 0;
		numGaps = 0;
		for(char c : col3) {
			if(c == 'A') numAs = numAs + 1;
			else if(c == 'C') numCs = numCs + 1;
			else if(c == '~') numGaps = numGaps + 1;
			else if(c == 'T') numTs = numTs + 1;

		}
		assertEquals(2, numCs);
		assertEquals(2, numGaps);

		numAs = 0;
		numCs = 0;
		numTs = 0;
		numGaps = 0;
		for(char c : col6) {
			if(c == 'A') numAs = numAs + 1;
			else if(c == 'C') numCs = numCs + 1;
			else if(c == '~') numGaps = numGaps + 1;
			else if(c == 'T') numTs = numTs + 1;

		}
		assertEquals(1, numTs);
		assertEquals(3, numGaps);

		numAs = 0;
		numCs = 0;
		numTs = 0;
		numGaps = 0;
		for(char c : col9) {
			if(c == 'A') numAs = numAs + 1;
			else if(c == 'C') numCs = numCs + 1;
			else if(c == '~') numGaps = numGaps + 1;
			else if(c == 'T') numTs = numTs + 1;

		}
		assertEquals(1, numAs);
		assertEquals(3, numGaps);
	}

	@Test
	/**
	 * Test of the maskSeqeuence and unmaskedColumns methods, of class MultipleAlignment.
	 *                     1 1
	 * 0 1 2 3 4 5 6 7 8 9 0 1
	 * A T G C
	 *     G C - T
	 *         A T ~ T G A
	 *           T T T
	 * Removed 3 and 4
	 */
	public void testRemoveSequence() {
		System.out.println("* TEST: removeSequence");
		alignment.maskSequence(seq3);
		alignment.maskSequence(seq4);


		char[] col4 = alignment.unmaskedColumn(0);
		assertEquals(2, col4.length);
		int numAs = 0;
		int numCs = 0;
		int numTs = 0;
		int numGaps = 0;
		for(char c : col4) {
			if(c == 'A') numAs = numAs + 1;
			else if(c == 'C') numCs = numCs + 1;
			else if(c == '~') numGaps = numGaps + 1;
			else if(c == 'T') numTs = numTs + 1;
		}
		//assertEquals(0, numAs);
		//assertEquals(1, numGaps);

		assertEquals(10, alignment.getLength());
		assertEquals(2, alignment.getNumUnmaskedSequences());
	}

	@Test
	/**
	 * Test of the removeSequence method, of class MultipleAlignment.
	 *                     1 1
	 * 0 1 2 3 4 5 6 7 8 9 0 1
	 * A T G C
	 *     G C - T
	 *         A T ~ T G A
	 *           T T T
	 * Removed 3 and 4
	 */
	public void testFilterRedundantSeqs() {
			alignment.maskUnmaskedRedundantSeqsByUnmasked();
			assertEquals(3, alignment.getNumUnmaskedSequences());
		}


	//@Test
	/**
	 * Test of the conflictGraph method, of class MultipleAlignment.
	 *                     1 1
	 * 0 1 2 3 4 5 6 7 8 9 0 1
	 * A T G C
	 *     G C - T
	 *         A T ~ T G A
	 *           T T T
	 */
	//public void testConflictGraph() {
	//		SimpleGraph<Sequence, DefaultEdge> conflictGraph = alignment.conflictGraph();
	//		assertEquals(true, conflictGraph.containsEdge(seq2, seq3));
	//		assertEquals(false, conflictGraph.containsEdge(seq3, seq4));
	//}

	@Test
	/**
	 * Test of the addSNP method, of the class MultipleAlignemtn
	 *                     1 1
	 * 0 1 2 3 4 5 6 7 8 9 0 1
	 * A T G C
	 *     G C - T
	 *         A T ~ T G A
	 *           T T T
	 * Added: a SNP at position 4
	 */
	public void testAddSNP() {
		System.out.println("* TEST: addSNP");
		assertEquals(true, seq2.containsSNP(newSNP));
		assertEquals(false, seq1.containsSNP(newSNP));
		assertEquals(true, newSNP.coversSequence(seq3));
		assertEquals(false, newSNP.coversSequence(seq1));
	}




	@Test
	/**
	 * Test of the getWeightedBipartiteGraphFromSequences method, of the class MultipleAlignment
	 *                     1 1
	 * 0 1 2 3 4 5 6 7 8 9 0 1
	 * A T G C
	 *     G C - T
	 *         A T ~ T G A
	 *           T T T
	 * Masked: Sequence 4,
	 * Already Added: SNP at position 4
	 */
	public void testGetWeightedBipartiteGraphFromSequences() {
		System.out.println("* TEST: getWeightedBipartiteGraphFromSequences");
		alignment.maskSequence(seq4);
		ArrayList<Sequence> unMasked = new ArrayList<Sequence>();
		unMasked.add(seq1); unMasked.add(seq2); unMasked.add(seq3);

		WeightedBipartiteGraph bpgraph = alignment.getWeightedBipartiteGraphFromSequences(unMasked, true);


		SimpleWeightedGraph wgraph = bpgraph.getGraph();
		Set<GraphNode> vertexSet = wgraph.vertexSet();
		assertEquals(6, vertexSet.size());
		ArrayList<GraphNode> leftNodes = bpgraph.getLeftNodes();
		ArrayList<GraphNode> rightNodes = bpgraph.getRightNodes();
		for(GraphNode lnodei: leftNodes) {
			for(GraphNode lnodej: leftNodes) {
				assertEquals(null, wgraph.getEdge(lnodei, lnodej));
			}
		}
		for(GraphNode rnodei: rightNodes) {
			for(GraphNode rnodej: rightNodes) {
				assertEquals(null, wgraph.getEdge(rnodei, rnodej));
			}
		}

		// NOTE: getEdgeWeight(null) returns 1.0! Durrr
		for(GraphNode lnode: leftNodes) {
			for(GraphNode rnode: rightNodes) {
				Sequence lnodeSeq = (Sequence)lnode.getData();
				Sequence rnodeSeq = (Sequence)rnode.getData();
				if(lnodeSeq == rnodeSeq) {
					assertTrue(wgraph.getEdgeWeight(wgraph.getEdge(lnode,rnode))< 0.5);
				}
				if(lnodeSeq.getName() == "seq1" && rnodeSeq.getName() == "seq2") {
					assertTrue(wgraph.getEdgeWeight(wgraph.getEdge(lnode,rnode)) == 1.0 && wgraph.getEdge(lnode,rnode) != null);
				}
				if(lnodeSeq.getName() == "seq2" && rnodeSeq.getName() == "seq3") {
					assertTrue(null == wgraph.getEdge(lnode, rnode));
				}
			}
		}

	}




	@Test
	/**
	 * Test of the getWeightedBipartiteGraphFromSequences method, of the class MultipleAlignment
	 *                     1 1
	 * 0 1 2 3 4 5 6 7 8 9 0 1
	 * A T G C
	 *     G C - T
	 *         A T ~ T G A
	 *           T T T
	 * Masked: Sequence 4,
	 * Already Added: SNP at position 4
	 */
	public void testMaxWeightedBipartiteMatching() {
		System.out.println("* TEST: maxWeightedBipartiteMatching");
		alignment.maskSequence(seq4);
		ArrayList<Sequence> unMasked = new ArrayList<Sequence>();
		unMasked.add(seq1); unMasked.add(seq2); unMasked.add(seq3);

		WeightedBipartiteGraph bpgraph = alignment.getWeightedBipartiteGraphFromSequences(unMasked, true);

		Hungarian matcher = new Hungarian();
		WeightedBipartiteGraph matchedBipartiteGraph = matcher.maxWeightedBipartiteMatching(bpgraph);

		SimpleWeightedGraph matchedGraph = matchedBipartiteGraph.getGraph();
		ArrayList<GraphNode> leftNodes = matchedBipartiteGraph.getLeftNodes();
		ArrayList<GraphNode> rightNodes = matchedBipartiteGraph.getRightNodes();

		List<Set<GraphNode>> ccs = new ConnectivityInspector(matchedGraph).connectedSets();

		// Number of connected components: 4
		assertEquals(4, ccs.size());

		// No edge should be between nodes on the left or right
		for(GraphNode inode: leftNodes) {
			for(GraphNode jnode: leftNodes) {
				assertEquals(null, matchedGraph.getEdge(inode, jnode));
			}
		}
		for(GraphNode inode: rightNodes) {
			for(GraphNode jnode: rightNodes) {
				assertEquals(null, matchedGraph.getEdge(inode, jnode));
			}
		}

		// two edges
		assertEquals(2, matchedGraph.edgeSet().size());

		// There are two valid maximum matchings for this example,
		// in one seq1(index0) maches to seq2 (index1), or seq1 (index0) matches to seq3 (index2)
		assertTrue(matchedGraph.getEdge(leftNodes.get(0), rightNodes.get(1)) != null || matchedGraph.getEdge(leftNodes.get(0), rightNodes.get(2)) != null);

		if(matchedGraph.getEdge(leftNodes.get(0), rightNodes.get(1)) != null) {
			// If the edge from seq1 to seq2 was selected
			assertTrue(matchedGraph.getEdge(leftNodes.get(0), rightNodes.get(2)) == null);
			assertTrue(matchedGraph.getEdge(leftNodes.get(2), rightNodes.get(2)) != null);
		}
		else if(matchedGraph.getEdge(leftNodes.get(0), rightNodes.get(2)) != null) {
			// If the edge from seq1 to seq3 was selected
			assertTrue(matchedGraph.getEdge(leftNodes.get(0), rightNodes.get(1)) == null);
			assertTrue(matchedGraph.getEdge(leftNodes.get(1), rightNodes.get(1)) != null);
		}
		else
		{
			// one of the above should be true
			assertTrue(false);
		}
	}





	@Test
	/**
	 * Test of the addSequence method, of class MultipleAlignment.
	 *                     1 1
	 * 0 1 2 3 4 5 6 7 8 9 0 1
	 * A T G C
	 *     G C - T
	 *         A T ~ T G A
	 *           T T T
	 * Added:
	 *             T C G A T G
	 *               C G
	 */
	public void testAddSequence() {
		System.out.println("* TEST: addSequence");
		Sequence seq5 = new Sequence();
		seq5.setSequence("TCGATG"); seq5.setName("seq5"); seq5.setStartPosition(6);
		Sequence seq6 = new Sequence();
		seq6.setSequence("CG"); seq6.setName("seq6"); seq6.setStartPosition(7);

		alignment.addSequence(seq5);
		alignment.addSequence(seq6);

		char[] col9 = alignment.unmaskedColumn(9);

		int numAs = 0;
		int numCs = 0;
		int numTs = 0;
		int numGaps = 0;
		for(char c : col9) {
			if(c == 'A') numAs = numAs + 1;
			else if(c == 'C') numCs = numCs + 1;
			else if(c == '~') numGaps = numGaps + 1;
			else if(c == 'T') numTs = numTs + 1;

		}
		assertEquals(2, numAs);
		assertEquals(4, numGaps);

		assertEquals(12, alignment.getLength());
		assertEquals(6, alignment.getNumUnmaskedSequences());
	}


	/**
	 * Test of length method, of class MultipleAlignment.
	 */
	@Test
	public void testGetLength() {
		System.out.println("* TEST: getLength");
		assertEquals(10, alignment.getLength());
	}



}