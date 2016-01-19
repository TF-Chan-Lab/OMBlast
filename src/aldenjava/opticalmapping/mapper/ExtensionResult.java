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


package aldenjava.opticalmapping.mapper;

import java.util.Comparator;
import java.util.LinkedHashMap;

import aldenjava.opticalmapping.Cigar;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;

/**
 * The <code>ExtensionResult</code> specifies the result of extension of the data from a specific location. This includes location of extension on both reference and data, alignment details (score and CIGAR string). However, strand information of the data is not included.
 * 
 * @author Alden
 *
 */
public class ExtensionResult implements Comparable<ExtensionResult> {
	public final String refName;
	public final int startfinalrefpos;
	public final int startfinalfragmentpos;
	public final int stopfinalrefpos;
	public final int stopfinalfragmentpos;
	public final String precigar;
	public final double score;
	/**
	 * This is the scale used in extension. However, it is not used in the conversion to alignment result.    
	 */
	public final double scale;

	public ExtensionResult(String refName, int startfinalrefpos, int startfinalfragmentpos, int stopfinalrefpos, int stopfinalfragmentpos, String precigar, double score, double scale) {
		this.refName = refName;
		this.startfinalrefpos = startfinalrefpos;
		this.startfinalfragmentpos = startfinalfragmentpos;
		this.stopfinalrefpos = stopfinalrefpos;
		this.stopfinalfragmentpos = stopfinalfragmentpos;
		this.precigar = precigar;
		this.score = score;
		this.scale = scale;
	}

	/**
	 * Converts the extension result to alignment result
	 * 
	 * @param data
	 *            data used in the extension
	 * @param optrefmap
	 *            the reference
	 * @param mappedstrand
	 *            alignment orientation
	 * @return the converted alignment result
	 */
	public OptMapResultNode toAlignment(DataNode data, LinkedHashMap<String, DataNode> optrefmap, int mappedstrand) {
		DataNode ref = optrefmap.get(this.refName);
		Cigar cigar = new Cigar(precigar);
		if (mappedstrand == 1)
			return new OptMapResultNode(data, ref.getGenomicPos(startfinalrefpos, stopfinalrefpos, true), mappedstrand, startfinalrefpos, stopfinalrefpos, startfinalfragmentpos, stopfinalfragmentpos,
					cigar, score, -1);
		else if (mappedstrand == -1)
			return new OptMapResultNode(data, ref.getGenomicPos(startfinalrefpos, stopfinalrefpos, true), mappedstrand, startfinalrefpos, stopfinalrefpos, data.getTotalSegment()
					- startfinalfragmentpos - 1, data.getTotalSegment() - stopfinalfragmentpos - 1, cigar, score, -1);
		else
			return null;

	}

	@Override
	public int compareTo(ExtensionResult result) {
		if (result == null)
			return 1;
		else
			return (Double.compare(this.score, result.score));

	}

	public static Comparator<ExtensionResult> comparator() {
		return new Comparator<ExtensionResult>() {
			@Override
			public int compare(ExtensionResult e1, ExtensionResult e2) {
				int now = e1.refName.compareTo(e2.refName);
				if (now != 0)
					return now;
				else {
					now = Integer.valueOf(e1.startfinalrefpos).compareTo(Integer.valueOf(e2.startfinalrefpos));
					if (now != 0)
						return now;
					else
						return Integer.valueOf(e1.stopfinalrefpos).compareTo(Integer.valueOf(e2.stopfinalrefpos));
				}
			}
		};
	}

}
