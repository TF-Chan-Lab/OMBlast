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


package aldenjava.opticalmapping.mapper.clustermodule;

import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;

/**
 * A filter class to check if the partial alignment can be joined based on the relationship (not based on overlapping)
 * Indel: sameStrand == true; oppositeStrand == false
 * Inversion: sameStrand == false; oppositeStrand == true
 * Indel+Inversion OR Translocation: sameStrand == false; oppositeStrand == false
 * @author Alden
 *
 */
public class PathBuilderFilter {
	public boolean sameStrand;
	public boolean oppositeStrand;
	public boolean correctOrder;
	public long closeReference;
	public long closeFragment;
	
	
	public PathBuilderFilter(boolean sameStrand, boolean oppositeStrand, boolean correctOrder, long closeReference, long closeFragment) {
		if (sameStrand && oppositeStrand)
			throw new IllegalArgumentException("Cannot require both same strand and opposite strand!");
		this.sameStrand = sameStrand;
		this.oppositeStrand = oppositeStrand;
		this.correctOrder = correctOrder;
		this.closeReference = closeReference;
		this.closeFragment = closeFragment;
	}

	public boolean checkPass(OptMapResultNode map1, OptMapResultNode map2) {
		if ((!map1.isUsed()) || !map2.isUsed())
			return false;
		// same strand
		if (sameStrand)
			if (map1.mappedstrand != map2.mappedstrand)
				return false;
		if (oppositeStrand)
			if (map1.mappedstrand == map2.mappedstrand)
				return false;
		// ordering;		
		if (correctOrder) {
			if (map1.mappedstrand == 1)
				if (map1.mappedRegion.start > map2.mappedRegion.start)
					return false;
			if (map1.mappedstrand == -1)
				if (map2.mappedRegion.start > map1.mappedRegion.start)
					return false;
		}
		// close Reference Distance
		if (closeReference >= 0)
			if (!map1.isRefClose(map2, closeReference))
				return false;
		// close fragment distance
		if (closeFragment >= 0)
			if (!map1.isFragClose(map2, closeFragment))
				return false;
		//Criteria passed.
		return true;
	}
}