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
public class SequenceTest {
	public Sequence seq1;
	public Sequence seq2;
	public Sequence seq3;
	public Sequence seq4;

	// 0 1 2 3 4 5 6 7 8 9
	// A T G C
	//     G C - T
	//         A T ~ T G A
	//           T T T
    public SequenceTest() {
		 seq1 = new Sequence();
		 seq2 = new Sequence();
		 seq3 = new Sequence();
		 seq4 = new Sequence();
		 seq1.setSequence("ATGC"); seq1.setName("seq1"); seq1.setStartPosition(0);
		 seq2.setSequence("GC-T"); seq1.setName("seq2"); seq2.setStartPosition(2);
		 seq3.setSequence("AT~TGA"); seq1.setName("seq3"); seq3.setStartPosition(4);
		 seq4.setSequence("TTT"); seq1.setName("seq4"); seq4.setStartPosition(5);
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
	 * Test of baseAtAlignmentPosition method, of class Sequence.
	 */
	@Test
	public void testGetBaseAtAlignmentPosition() {
		System.out.println("* TEST: getBaseAtAlignmentPosition");
		assertEquals('A', seq1.baseAtAlignmentPosition(0));
		assertEquals('-', seq2.baseAtAlignmentPosition(4));
		assertEquals('~', seq3.baseAtAlignmentPosition(6));
		assertEquals('~', seq1.baseAtAlignmentPosition(4));
	}

	/**
	 * Test of basesFromAlignmentPositionTo method, of class Sequence.
	 */
	 @Test
	 public void testGetBasesFromAlignmentPositionTo() {
		 System.out.println("* TEST: getBasesFromAlignmentPositionTo");
		 assertEquals("ATGC", seq1.basesFromAlignmentPositionTo(0,3));
		 assertEquals("ATGC~~", seq1.basesFromAlignmentPositionTo(0,5));
		 assertEquals("TGC", seq1.basesFromAlignmentPositionTo(1,3));
		 assertEquals("C-T", seq2.basesFromAlignmentPositionTo(3,5));
	 }

	/**
	 * Test of containsPosition method, of class Sequence.
	 */
	@Test
	public void testContainsPosition() {
		System.out.println("* TEST: containsPosition");
		assertEquals(true, seq1.containsPosition(0));
		assertEquals(true, seq1.containsPosition(3));
		assertEquals(false, seq1.containsPosition(4));
		assertEquals(true, seq3.containsPosition(4));
		assertEquals(true, seq3.containsPosition(6));
		assertEquals(false, seq3.containsPosition(3));
	}

	/**
	 * Test of covers method, of class Sequence.
	 */
	@Test
	public void testCovers() {
		System.out.println("* TEST: covers");
		assertEquals(true, seq3.covers(seq4));
		assertEquals(false, seq4.covers(seq3));
		assertEquals(false, seq1.covers(seq2));
	}

	/**
	 * Test of hasGaps method, of class Sequence.
	 */
	@Test
	public void testHasGaps() {
		System.out.println("* TEST: hasGaps");
		assertEquals(false, seq1.hasGaps());
		assertEquals(false, seq2.hasGaps());
		assertEquals(true, seq3.hasGaps());
	}

	/**
	 * Test of startsBefore method, of class Sequence.
	 */
	@Test
	public void testStartsBefore() {
		System.out.println("* TEST: startsBefore");
		assertEquals(true, seq1.startsBefore(seq2));
		assertEquals(true, seq2.startsBefore(seq3));
		assertEquals(false, seq2.startsBefore(seq1));
		assertEquals(true, seq3.startsBefore(seq4));
		assertEquals(false, seq4.startsBefore(seq3));
	}

	/**
	 * Test of overlaps method, of class Sequence.
	 */
	@Test
	public void testOverlaps() {
		System.out.println("* TEST: overlaps");
		assertEquals(true, seq1.overlaps(seq2));
		assertEquals(true, seq2.overlaps(seq3));
		assertEquals(false, seq1.overlaps(seq3));
		assertEquals(true, seq3.overlaps(seq4));
		assertEquals(true, seq2.overlaps(seq4));
	}


	/**
	 * Test of endPosition method, of class Sequence.
	 */
	@Test
	public void testGetEndPosition() {
		System.out.println("* TEST: getEndPosition");
		assertEquals(3, seq1.endPosition());
		assertEquals(5, seq2.endPosition());
		assertEquals(9, seq3.endPosition());
		assertEquals(7, seq4.endPosition());
	}



}