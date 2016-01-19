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

import java.util.List;

import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;

/**
 * The joined/clustered alignment results
 * @author Alden
 *
 */
public class ClusteredResult implements Comparable<ClusteredResult> {
	public List<OptMapResultNode> updatedResult = null;
	public Double score = null;

	public ClusteredResult() {
	}

	public void process(VirtualMapProcessor vmProcessor) {
		score = vmProcessor.calcScore(updatedResult);
	}

	public void importUpdatedResult(List<OptMapResultNode> updatedResult) {
		this.updatedResult = updatedResult;
	}

	@Override
	public int compareTo(ClusteredResult cr) {
		return Double.compare(this.score, cr.score);
	}

	public double getFragRatio() {
		double total = 0;
		for (OptMapResultNode result : updatedResult)
			total += result.getSubFragRatio();
		return total;
	}

	public double getMapSigRatio() {
		double total = 0;
		for (OptMapResultNode result : updatedResult)
			total += result.getMapSigRatio();
		return total;
	}

}