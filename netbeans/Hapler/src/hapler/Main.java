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
public class Main {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		//System.out.println("Netbeans, baby.");
		//String fileName = args[0];

		IOHandler handler = new IOHandler();

		int parseResult = -1;

		if(args.length == 0) {
			// Use Default Args as specified in IOHandler, if none given
			parseResult = handler.parseOptions(new String[0]);

			//For profiling
			/*String[] testArgs = new String[10];
			testArgs[0] = "--input";
			//testArgs[1] = "/Users/soneil/Documents/research/bioInformatics/kHaplotyping/software/testData/agambiae.1925488427.sam";
			testArgs[1] = "/Users/soneil/Documents/research/bioInformatics/kHaplotyping/software/testData/bigTest.contig";
			testArgs[2] = "--random-repetitions";
			testArgs[3] = "1";
			testArgs[4] = "--allow-gaps";
			testArgs[5] = "split";
			testArgs[6] = "--snp-caller";
			testArgs[7] = "454";
			testArgs[8] = "--alignment-type";
			//testArgs[9] = "sam";
			testArgs[9] = "tigr";
			parseResult = handler.parseOptions(testArgs);*/
		}
		else {
			parseResult = handler.parseOptions(args);
		}

		if(parseResult != 0) {
			return;
		}
		try{
			handler.execAndPrintSummaryStats();
		}
		catch(Exception e) {
			System.err.println("#############################");
			System.err.println("Error: " + e.getMessage());
			System.err.println("#############################");
			System.err.println();
			e.printStackTrace();
		}


		}


}
