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
public class HaplotypeTest {
	private MultipleAlignment alignment;
	private Sequence seq1;
	private Sequence seq2;
	private Sequence seq3;
	private Sequence seq4;
	private Sequence seq5;
	private Sequence seq6;
	private Sequence seq7;
	private SNP snp1;
	private SNP snp2;
	private SNP snp3;
	private Haplotype hap1;
	private Haplotype hap2;

	/**
	 * 0 1 2 3 4 5 6 7 8 9 0
	 *     * .   *
	 * ###################   HAP1, total support: 4, pieces: 1, Consensus: ~ATGCTGAT~
	 *   A T G C T G         Supports 1
	 *     T - C T           Supports 1
	 *       G C             Supports 1,2: is masked.
	 *         C T G A T     Supports 1
	 * ###################   HAP2, total support: 3, pieces: 2, Consensus: TACGCAGATA
	 *     C G C A           Supports 2
	 *         C A G A T A   Supports 2
	 *               A T A C Supports 1,2
	 *     * *   *
	 */

    public HaplotypeTest()  {
		seq1 = new Sequence();
		seq2 = new Sequence();
		seq3 = new Sequence();
		seq4 = new Sequence();
		seq5 = new Sequence();
		seq6 = new Sequence();
		seq7 = new Sequence();
		seq1.setSequence("ATGCTG"); seq1.setName("seq1"); seq1.setStartPosition(1);
		seq2.setSequence("T-CT"); seq2.setName("seq2"); seq2.setStartPosition(2);
		seq3.setSequence("GC"); seq3.setName("seq3"); seq3.setStartPosition(3);
		seq4.setSequence("CTGAT"); seq4.setName("seq4"); seq4.setStartPosition(4);
		seq5.setSequence("CGCA"); seq5.setName("seq5"); seq5.setStartPosition(2);
		seq6.setSequence("CAGATA"); seq6.setName("seq6"); seq6.setStartPosition(4);
		seq7.setSequence("ATAC"); seq7.setName("seq7"); seq7.setStartPosition(7);
		alignment = new MultipleAlignment();
		alignment.addSequence(seq1);
		alignment.addSequence(seq2);
		alignment.addSequence(seq3);
		alignment.addSequence(seq4);
		alignment.addSequence(seq5);
		alignment.addSequence(seq6);
		alignment.addSequence(seq7);
		snp1 = new SNP(2); alignment.addSNP(snp1);
		snp2 = new SNP(3); alignment.addSNP(snp2);
		snp3 = new SNP(5); alignment.addSNP(snp3);

		alignment.maskUnmaskedRedundantSeqsByUnmasked();

		hap1 = new Haplotype(alignment);
		hap1.addSequence(seq1);
		hap1.addSequence(seq2);
		//hap1.addSequence(seq3);
		hap1.addSequence(seq4);

	 	hap2 = new Haplotype(alignment);
		hap2.addSequence(seq5);
		hap2.addSequence(seq6);
		hap2.addSequence(seq7);

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
	 * Test of numPieces method, of class Haplotype.
	 */
	@Test
	public void testGetNumPieces() {
		System.out.println("* TEST: getNumPieces");
		assertEquals(1, hap1.numPieces());
		assertEquals(2, hap2.numPieces());
	}

	/**
	 * Test of consensusHumanReadable method, of class Haplotype.
	 */
	@Test
	public void testGetConsensus() {
		System.out.println("* TEST: getConsensus");
		assertTrue("~ATGCTGAT~~".compareTo(hap1.consensusHumanReadable()) == 0);
		assertTrue("~~CGCAGATAC".compareTo(hap2.consensusHumanReadable()) == 0);
		//System.out.println(hap2.consensusHumanReadable());
	}

	/**
	 * Test of numInconsistentSNPs method, of class Haplotype.
	 */
	@Test
	public void testGetNumInconsistentUnmaskedSnps() {
		System.out.println("* TEST: getNumInconsistenUnmaskedSnps");
		assertEquals(1, hap1.numInconsistentSNPs());
		assertEquals(0, hap2.numInconsistentSNPs());
	}

	/**
	 * Test of consistentSNPsCovered method, of class Haplotype.
	 */
	@Test
	public void testGetUnmaskedConsistentSnpsCovered() {
		System.out.println("* TEST: getUnmaskedConsistentSnpsCovered");
		ArrayList<SNP> hap1ConsSnps = hap1.consistentSNPsCovered();
		ArrayList<SNP> hap2ConsSnps = hap2.consistentSNPsCovered();


		assertEquals(hap1ConsSnps.size(), 2);
		assertEquals(hap2ConsSnps.size(), 3);
		assertTrue(hap1ConsSnps.contains(snp1));
		assertTrue(hap1ConsSnps.contains(snp3));
		assertTrue(hap2ConsSnps.contains(snp1));
		assertTrue(hap2ConsSnps.contains(snp2));
		assertTrue(hap2ConsSnps.contains(snp3));
	}

	/**
	 * Test of sequenceConsistentAtSNPs method, of class Haplotype.
	 */
	@Test
	public void testSequenceConsistentAtSnps() {
		System.out.println("* TEST: sequenceConsistentAtSnps");
		ArrayList<SNP> snp13 = new ArrayList<SNP>(); snp13.add(snp1); snp13.add(snp3);
		ArrayList<SNP> snp123 = new ArrayList<SNP>(); snp123.add(snp1); snp123.add(snp2); snp123.add(snp3);

		assertTrue(hap1.sequenceConsistentAtSNPs(seq1, snp13));
		assertFalse(hap2.sequenceConsistentAtSNPs(seq1, snp123));
		assertTrue(hap1.sequenceConsistentAtSNPs(seq2, snp13));
		assertFalse(hap2.sequenceConsistentAtSNPs(seq2, snp123));
		assertTrue(hap1.sequenceConsistentAtSNPs(seq3, snp13));
		assertTrue(hap2.sequenceConsistentAtSNPs(seq3, snp123));
		assertTrue(hap1.sequenceConsistentAtSNPs(seq4, snp13));
		assertFalse(hap2.sequenceConsistentAtSNPs(seq4, snp123));
		assertFalse(hap1.sequenceConsistentAtSNPs(seq5, snp13));
		assertTrue(hap2.sequenceConsistentAtSNPs(seq5, snp123));
		assertFalse(hap1.sequenceConsistentAtSNPs(seq6, snp13));
		assertTrue(hap2.sequenceConsistentAtSNPs(seq6, snp123));
		assertTrue(hap1.sequenceConsistentAtSNPs(seq7, snp13));
		assertTrue(hap2.sequenceConsistentAtSNPs(seq7, snp123));
	}

	/**
	 * Test of numNonRedundantSequences method, of class Haplotype.
	 */
	@Test
	public void testGetNumNonRedundantSequences() {
		System.out.println("* TEST: getNumNonRedundantSequences");
		assertEquals(3, hap1.numNonRedundantSequences());
		assertEquals(3, hap2.numNonRedundantSequences());
	}

	/**
	 * Test of numRedundantSequences method, of class Haplotype.
	 */
	@Test
	public void testGetNumRedundantSequences() {
		System.out.println("* TEST: getNumRedundantSequences");
		assertTrue(seq3.isMasked());
		assertEquals(1, seq1.getMasking().size());
		assertEquals(1, hap1.numRedundantSequences());
		assertEquals(0, hap2.numRedundantSequences());
	}

	/**
	 * Test of numSNPsCovered method, of class Haplotype.
	 */
	@Test
	public void testGetNumUnmaskedSNPsCovered() {
		System.out.println("* TEST: getNumUnmaskedSNPsCovered");
		assertEquals(3, hap1.numSNPsCovered());
		assertEquals(3, hap2.numSNPsCovered());
	}

	/**
	 * Test of startPos method, of class Haplotype.
	 */
	@Test
	public void testGetStartPos() {
		System.out.println("* TEST: getStartPos");
		assertEquals(1, hap1.startPos());
		assertEquals(2, hap2.startPos());
	}

	/**
	 * Test of endPos method, of class Haplotype.
	 */
	@Test
	public void testGetEndPos() {
		System.out.println("* TEST: getEndPos");
		assertEquals(8, hap1.endPos());
		assertEquals(10, hap2.endPos());
	}

	/**
	 * Test of allSequences method, of class Haplotype.
	 */
	@Test
	public void testGetAllSequences() {
		System.out.println("* TEST: getAllSequences");
		assertEquals(4, hap1.allSequences().size());
		assertEquals(3, hap2.allSequences().size());
	}

	/**
	 * Test of length method, of class Haplotype.
	 */
	@Test
	public void testGetLength() {
		System.out.println("* TEST: getAllSequences");
		assertEquals(8, hap1.length());
		assertEquals(9, hap2.length());
	}

}