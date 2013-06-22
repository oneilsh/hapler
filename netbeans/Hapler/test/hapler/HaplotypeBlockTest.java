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
import java.util.HashMap;
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
public class HaplotypeBlockTest {
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
	private HaplotypeBlock block1;

	/**
	 * 0 1 2 3 4 5 6 7 8 9 0
	 *     * .   *
	 * ###################   HAP1, total support: 4, pieces: 1,
	 *   A T G C T G         Supports 1
	 *     T - C T           Supports 1
	 *       G C             Supports 1,2: is masked.
	 *         C T G A T     Supports 1
	 * ###################   HAP2, total support: 3, pieces: 1,
	 *     C G C A           Supports 2
	 *         C A G A T A   Supports 2
	 *               A T A C Supports 1,2
	 *     * *   *
	 */


    public HaplotypeBlockTest() {
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

		block1 = new HaplotypeBlock(alignment);
		block1.addHaplotype(hap1);
		block1.addHaplotype(hap2);
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
	 * Test of numPieces method, of class HaplotypeBlock.
	 */
	@Test
	public void testGetNumPieces() {
		System.out.println("* TEST: getNumPieces");
		assertEquals(3, block1.numPieces());
	}

	/**
	 * Test of numInconsistentSNPs method, of class HaplotypeBlock.
	 */
	@Test
	public void testGetNumInconsistentUnmaskedSNPs() {
		System.out.println("* TEST: getNumInconsistentUnmaskedSNPs");
		assertEquals(1, block1.numInconsistentSNPs());
	}

	/**
	 * Test of lexicographicSize method, of class HaplotypeBlock.
	 */
	@Test
	public void testGetLexicographicSize() {
		System.out.println("* TEST: getLexicographicSize");
		ArrayList<Double> lexSize = block1.lexicographicSize();
		assertEquals(2, lexSize.size());
		assertEquals(3.0, lexSize.get(0), 0.0001);
		assertEquals(4.0, lexSize.get(1), 0.0001);
	}

	/**
	 * Test of haplotypeSupports method, of class HaplotypeBlock.
	 */
	@Test
	public void testGetHaplotypeSupports() {
		System.out.println("* TEST: getHaplotypeSupports");
		HashMap<Haplotype, Double> supports = block1.haplotypeSupports();
		assertEquals(4.0, supports.get(hap1), 0.0001);
		assertEquals(3.0, supports.get(hap2), 0.0001);
	}

	/**
	 * Test of allSequences method, of class HaplotypeBlock.
	 */
	@Test
	public void testGetAllSequences() {
		System.out.println("* TEST: getAllSequences");
		ArrayList<Sequence> seqs = block1.allSequences();
		assertEquals(7, seqs.size());
	}

}