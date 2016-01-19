/**************************************************************************
**  OMBlast
**  Software aligning optical maps
**  
**  Version 1.0 -- September 1, 2015
**  
**  Copyright (C) 2015 by Alden Leung, All rights reserved.
**  Contact:  aldenleung@link.cuhk.edu.hk
**  Organization:  Hong Kong Bioinformatics Centre, School of Life Sciences, The
**                 Chinese University of Hong Kong, Shatin, NT,
**                 Hong Kong SAR
**  
**  This file is part of OMBlast.
**  
**  OMBlast is free software; you can redistribute it and/or 
**  modify it under the terms of the GNU General Public License 
**  as published by the Free Software Foundation; either version 
**  3 of the License, or (at your option) any later version.
**  
**  OMBlast is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**  
**  You should have received a copy of the GNU General Public 
**  License along with OMBlast; if not, see 
**  <http://www.gnu.org/licenses/>.
**************************************************************************/


package aldenjava.opticalmapping;

import java.io.IOException;

import org.apache.commons.lang.StringUtils;

/**
 * A class storing Cigar string, which is composed of 'M', 'I' and 'D' for match, insertion and deletion of signal (not segment). A pre-cigar (e.g. 'MMMDMMIM') is the stored rather than the Cigar string ('3M1D2M1I1M'), for the sake of convenience in modification. This cigar accepts characters other than 'M', 'I' and 'D' for future development.
 * 
 * @author Alden
 *
 */
public class Cigar {

	private int indexM = 0;
	private int indexD = 0;
	private int indexI = 0;
	private StringBuilder precigar;

	/**
	 * Creates a blank <code>Cigar</code> without any information.
	 */
	public Cigar() {
		this.precigar = new StringBuilder();
		resetMDPIndex();
	}

	/**
	 * Creates a <code>Cigar</code> according to the specified <code>precigar</code>.
	 * 
	 * @param precigar
	 *            a string containing the precigar
	 */
	public Cigar(String precigar) {
		this.precigar = new StringBuilder(precigar);
		reCalcMDPIndex();
	}

	/**
	 * Constructs a new <code>Cigar</code>, initialized to match the values of the specified <code>Cigar</code>.
	 * 
	 * @param cigar
	 *            the <code>Cigar</code> from which to copy initial values to a newly constructed <code>Cigar</code>
	 */
	public Cigar(Cigar cigar) {
		this.precigar = new StringBuilder(cigar.precigar);
		this.indexM = cigar.indexM;
		this.indexI = cigar.indexI;
		this.indexD = cigar.indexD;
	}

	/**
	 * Import information into <code>Cigar</code> according to the specified <code>cigar</code>.
	 * 
	 * @param cigar
	 *            a string containing the cigar
	 * @see #importPrecigar(String)
	 */
	public void importCigar(String cigar) {
		StringBuilder precigar = new StringBuilder();
		StringBuilder recent = new StringBuilder();
		for (int i = 0; i < cigar.length(); i++) {
			char c = cigar.charAt(i);
			if (c >= '0' && c <= '9')
				recent.append(c);
			else {
				if (recent.length() == 0) {
					this.precigar = new StringBuilder();
					return;
				} else {
					int number = Integer.parseInt(recent.toString());
					for (int j = 0; j < number; j++)
						precigar.append(c);
					recent = new StringBuilder();
				}
			}
		}
		this.precigar = precigar;
		reCalcMDPIndex();
	}
	/**
	 * Imports information into <code>Cigar</code> according to the specified <code>precigar</code>.
	 * 
	 * @param precigar
	 *            a string containing the precigar
	 * @see #importCigar(String)
	 */
	public void importPrecigar(String precigar) {
		this.precigar = new StringBuilder(precigar);
		reCalcMDPIndex();
	}

	/**
	 * Appends any match, insertion or deletion to this <code>Cigar</code>
	 * @param c	Any character 'M', 'I' or 'D'
	 */
	public void append(char c) {
		precigar.append(c);
		addMDPIndex(c);
	}
	/**
	 * Appends a <code>Cigar</Cigar> to this <code>Cigar</code>
	 * @param cigar	another <code>Cigar</Cigar>
	 */
	public void append(Cigar cigar) {
		precigar.append(cigar.precigar);
		indexM += cigar.indexM;
		indexI += cigar.indexI;
		indexD += cigar.indexD;
	}

	private void addMDPIndex(char c) {
		switch (c) {
			case 'M':
				indexM++;
				break;
			case 'I':
				indexI++;
				break;
			case 'D':
				indexD++;
				break;
			default:
				;
		}
	}

	private void addMDPIndex(String precigar) {
		for (char c : precigar.toCharArray())
			addMDPIndex(c);
	}
	
	private void resetMDPIndex() {
		indexM = 0;
		indexI = 0;
		indexD = 0;
	}

	public void reCalcMDPIndex() {
		resetMDPIndex();
		addMDPIndex(this.precigar.toString());
	}

	public String getCigar() {
		char prevc = '~';
		int accumulate = 0;
		StringBuilder cigar = new StringBuilder();
		for (int i = 0; i < precigar.length(); i++) {
			char c = precigar.charAt(i);
			if (prevc != c) {
				if (accumulate > 0) {
					cigar.append(Integer.toString(accumulate));
					cigar.append(prevc);
				}
				prevc = c;
				accumulate = 1;
			} else
				accumulate++;
		}
		if (accumulate > 0) {
			cigar.append(Integer.toString(accumulate));
			cigar.append(prevc);
		}
		return cigar.toString();
	}

	public String getPrecigar() {
		return precigar.toString();
	}

	public Cigar getReverseCigar() {
		Cigar newcigar = new Cigar(this);
		newcigar.reverse();
		return newcigar;
	}

	public int getMatch() {
		return indexM;
	}

	public int getFP() {
		return indexI;
	}

	public int getFN() {
		return indexD;
	}

	public int getNumber(char matchedc) {
		int total = 0;
		for (char c : precigar.toString().toCharArray())
			if (c == matchedc)
				total++;
		return total;
	}

	public boolean equals(Cigar cigar) {
		return this.precigar.toString().equalsIgnoreCase(cigar.precigar.toString());
	}

	public double calcScore(double match, double fpp, double fnp) {
		return getMatch() * match - getFP() * fpp - getFN() * fnp;
	}

	public void removeFlankingUnmatch() {
		String tmp = precigar.toString();
		if (tmp.indexOf('M') == -1)
			this.precigar = new StringBuilder();
		else
			this.precigar = new StringBuilder(tmp.substring(tmp.indexOf('M'), tmp.lastIndexOf('M') + 1));
	}

	public void reverse() {
		precigar.reverse();
	}

	public Cigar reverseRefAndFrag() {
		StringBuilder newprecigar = new StringBuilder();
		for (char c : precigar.toString().toCharArray())
			newprecigar.append(c == 'M' ? 'M' : c == 'I' ? 'D' : 'I');
		Cigar cigar = new Cigar();
		cigar.importPrecigar(precigar.toString());
		return cigar;
	}

	@Override
	public String toString() {
		return this.getCigar();
	}

	public static Cigar newUnmapCigar(int insertion, int deletion) {
		if (insertion < 0)
			throw new IllegalArgumentException("Insertion < 0");
		if (deletion < 0)
			throw new IllegalArgumentException("Deletion < 0");
		Cigar cigar = new Cigar();
		cigar.importCigar(insertion + "I" + deletion + "D");
		return cigar;
	}

	public static String convertCIGAR(String precigarString) {
		char prevc = '~';
		int accumulate = 0;
		StringBuilder cigar = new StringBuilder();
		for (int i = 0; i < precigarString.length(); i++) {
			char c = precigarString.charAt(i);
			if (prevc != c) {
				if (accumulate > 0) {
					cigar.append(Integer.toString(accumulate));
					cigar.append(prevc);
				}
				prevc = c;
				accumulate = 1;
			} else
				accumulate++;
		}
		if (accumulate > 0) {
			cigar.append(Integer.toString(accumulate));
			cigar.append(prevc);
		}
		return cigar.toString();
	}

	public static String convertpreCIGAR(String cigar) {
		StringBuilder precigar = new StringBuilder();
		StringBuilder recent = new StringBuilder();
		for (int i = 0; i < cigar.length(); i++) {
			char c = cigar.charAt(i);
			if (c >= '0' && c <= '9')
				recent.append(c);
			else {
				if (recent.length() == 0)
					return "";
				else {
					int number = Integer.parseInt(recent.toString());
					for (int j = 0; j < number; j++)
						precigar.append(c);
					recent = new StringBuilder();
				}
			}
		}
		return precigar.toString();
	}

	public static int getCertainNumber(String cigar, char matchedc) {
		return getCertainNumberFromPrecigar(convertpreCIGAR(cigar), matchedc);
	}

	public static int getCertainNumberFromPrecigar(String precigar, char matchedc) {
		int total = 0;
		for (char c : precigar.toCharArray())
			if (c == matchedc)
				total++;
		return total;
	}

	public static String reverseRefAndFrag(String cigar) {
		String precigar = Cigar.convertpreCIGAR(cigar);
		StringBuilder newprecigar = new StringBuilder();
		for (char c : precigar.toCharArray())
			newprecigar.append(c == 'M' ? 'M' : c == 'I' ? 'D' : 'I');
		return Cigar.convertCIGAR(newprecigar.toString());
	}

	public static String reverse(String cigar) {
		String precigar = Cigar.convertpreCIGAR(cigar);
		precigar = StringUtils.reverse(precigar);
		return Cigar.convertCIGAR(precigar);
	}

	public static void main(String[] args) throws IOException {

		System.out.println("Please enter the cigar string to convert:");
		StringBuilder cigar = new StringBuilder();
		int c = System.in.read();
		while (c != -1) {
			cigar.append((char) c);
			if (System.in.available() == 0)
				break;
			c = System.in.read();
		}

		String precigar = Cigar.convertpreCIGAR(cigar.toString().trim());
		System.out.println("The precigar is:");
		System.out.println(precigar);
	}
}
