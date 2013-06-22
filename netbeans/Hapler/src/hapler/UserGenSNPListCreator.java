/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package hapler;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author soneil
 */
public class UserGenSNPListCreator extends AbstractSNPListCreator {
	private String filename;
	private HashMap<String,HashSet<Integer>> idsToSNPPosLists;


	public UserGenSNPListCreator(String theFilename) throws Exception {
		filename = theFilename;
		idsToSNPPosLists = new HashMap<String,HashSet<Integer>>();
		
		FileInputStream fstream = new FileInputStream(filename);
		DataInputStream in = new DataInputStream(fstream);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));

		HashMap<String,ArrayList<String>> idsToSeqArrays = new HashMap<String,ArrayList<String>>();

		String currentId = null;

		String strLine = null;
		while ((strLine = br.readLine()) != null) {
			if(!strLine.equals("")) {
				String lineArray[] = strLine.split("\\s+");
				if(lineArray.length < 2) throw new Exception("Error while parsing SNP list file: the line '" + strLine + "' should have at least two whitespace separated fields.");
				String id = lineArray[0];
				Integer position = Integer.parseInt(lineArray[1]);
				
				if(!idsToSNPPosLists.containsKey(id)) {
					idsToSNPPosLists.put(id, new HashSet<Integer>());
				}
				HashSet<Integer> currentSNPPositions = idsToSNPPosLists.get(id);
				if(currentSNPPositions.contains(position)) {
					System.err.println("#!WARNING: while parsing SNP list file: you are trying to add a SNP to alignment " + id + " at position " + position + ", but there is already one defined for that position at that alignment.");
				}
				else {
					currentSNPPositions.add(position);
				}
			}
		}
	
	}


	public ArrayList<SNP> computeSNPList(MultipleAlignment alignment) {
		ArrayList<SNP> snpList = new ArrayList<SNP>();

		String id = alignment.getName();
		if(!idsToSNPPosLists.containsKey(id)) {
			System.out.println("#!WARNING: There are apparently no SNP positions listed for the alignment " + id + " in the file " + filename + ", defining no SNPs for this alignment.");
		}
		else {
			HashSet<Integer> positions = idsToSNPPosLists.get(id);
			for(Integer pos : positions) {
				if(pos >= alignment.getLength()) {
					System.out.println("#!WARNING: while parsing the SNP list, we are trying to add a SNP at position " + pos + " to alignment " + alignment.getName() + ", but this alignment is only " + alignment.getLength() + " bases long. NOT adding this SNP.");
				}
				else if(!alignment.containsVariationAtPosition(pos)) {
					System.out.println("#!WARNING: while parsing the SNP list, we are trying to add a SNP at position " + pos + " to alignment " + alignment.getName() + ", but the alignment is non-variant at this position. NOT adding this SNP.");
				}
				else {
					snpList.add(new SNP(pos, alignment));
				}
			}
		}


		return snpList;
	}

}
