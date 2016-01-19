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


package aldenjava.opticalmapping.mapper.seeding;

import java.util.ArrayList;
import java.util.List;

import aldenjava.opticalmapping.Cigar;

/**
 * A class storing the matched reference and query kmers
 * @author Alden
 *
 */
public class Seed extends Kmer {
	public Kmer kmerpointer; // * Please don't use kmer as direct kmerpointer, create another new object for kmerpointer
	public double rangeUBound = -1;
	public double rangeLBound = -1;

	public Seed(String source, int pos, List<Long> sizelist, Kmer kmerpointer) {
		super(source, pos, sizelist);
		this.kmerpointer = kmerpointer;
	}

	public Seed(Kmer kmer, Kmer kmerpointer) {
		super(kmer);
		this.kmerpointer = kmerpointer;

	}

	public Cigar getCigar(boolean flankingMatch) {
		if (kmerpointer.getErrorNo() > 0 && super.getErrorNo() > 0) {
			System.err.println("Do not support errors in both ref-kmer and frag-kmer.");
			return null;
		}
		List<Integer> errorposlist;
		char c;
		if (kmerpointer.getErrorNo() > 0) {
			errorposlist = kmerpointer.errorposlist;
			c = 'I';
		} else if (super.getErrorNo() > 0) {
			errorposlist = super.errorposlist;
			c = 'D';
		} else {
			errorposlist = new ArrayList<Integer>();
			c = 'M';
		}

		StringBuilder precigar = new StringBuilder();
		int pointer = 0;
		for (int i = 0; i < realsizelist.size(); i++) {
			precigar.append('M');
			while (pointer < errorposlist.size() && errorposlist.get(pointer) == i) {
				precigar.append(c);
				pointer++;
			}
		}
		precigar.append('M');
		if (flankingMatch)
			return new Cigar(precigar.toString());
		else
			return new Cigar(precigar.substring(1, precigar.length() - 1));

	}

	public boolean limitRange(int measure, double ear) {
		double ubound = 1 + ear;
		double lbound = 1 - ear;
		for (int pos = 0; pos < k(); pos++) {
			double newubound = (kmerpointer.get(pos) + measure) / (double) get(pos);
			double newlbound = (kmerpointer.get(pos) - measure) / (double) get(pos);
			if (newubound < ubound)
				ubound = newubound;
			if (newlbound > lbound)
				lbound = newlbound;
		}
		rangeUBound = ubound;
		rangeLBound = lbound;
		return (ubound >= lbound);
	}
}
