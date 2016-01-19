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
 * 
 * This class contains a pair of joined parital alignments. This pair forms a node in the graph.  
 * @author Alden
 *
 */
public class ClusterPathNode {
	// Trimmed
	public OptMapResultNode firstMap;
	public OptMapResultNode secondMap;
	private int firstTrimmed;
	private int secondTrimmed;
	private OptMapResultNode firstTrimmedMap;
	private OptMapResultNode secondTrimmedMap;
	private double bestPrevScore = Double.NEGATIVE_INFINITY;
	private double currentScore = Double.NEGATIVE_INFINITY;
	public ClusterPathNode previousPath = null;
	
	public ClusterPathNode(OptMapResultNode firstMap, OptMapResultNode secondMap,
			int firstTrimmed, int secondTrimmed, OptMapResultNode firstTrimmedMap, OptMapResultNode secondTrimmedMap) {
		this.firstMap = firstMap;
		this.secondMap = secondMap;
		this.firstTrimmed = firstTrimmed;
		this.secondTrimmed = secondTrimmed;
		this.firstTrimmedMap = firstTrimmedMap;
		this.secondTrimmedMap = secondTrimmedMap;
	}
	/**
	 * Checks if the current instance can link to the target <code>ClusterPathNode</code>. Four criteria are listed. 
	 * <li> The current instance cannot be linked to itself </li>
	 * <li> Second map of the current instance points at the first map of the target </li>
	 * <li> Minimum match of a trimmed map must be larger or equal to <code>minMatch</code> </li>
	 * <li> The score of the link between this instance and target <code>ClusterPathNode</code> is greater than the score of existing <code>previousPath</code> </li> 
	 * @param cp The target <code>ClusterPathNode</code>
	 * @param minMatch Minimum match of a trimmed map
	 * @return <code>true</code> if the two <code>ClusterPathNode</code> can be linked 
	 */
	public boolean canLink(ClusterPathNode cp, int minMatch) {
		if (this == cp) // cannot link to itself
			return false;
		if (!(secondMap == cp.firstMap))
			return false;
		if (secondMap.cigar.getMatch() - cp.firstTrimmed - this.secondTrimmed < minMatch)
			return false;
		if (currentScore <= cp.bestPrevScore)
			return false;
		return true;
	}
	
	/**
	 * Assigns <code>previousPath</code> according to the target <code>ClusterPathNode</code>. Once should only do this if <code>cp.canLink(this)</code> returns <code>true</code>
	 * 
	 * @param cp the target <code>ClusterPathNode</code>
	 */
	public void assignPreviousPath(ClusterPathNode cp) {
		this.previousPath = cp;
		this.bestPrevScore = cp.currentScore;
	}
	
	/**
	 * Update the current score according to the scoring scheme in <code>vmProcessor</code>. 
	 * The score calculation requires either the <code>firstMap</code> is <code>null</code>, i.e. this is the start of the graph,
	 * or the previous path is not <code>null</code>. Score is updated to negative infinity if the above criteria is not fulfilled.   
	 * @param vmProcessor an instance of <code>VirtualMapProcessor</code> which provides the scoring scheme
	 */
	public void updateScore(VirtualMapProcessor vmProcessor) {
		// Trim1       ====================
		// Trim2    ====================
		// Trim1+2     =================
		// NoTrim   =======================
		// Score(Trim1+2) ~= Trim1 + Trim2 - NoTrim	(Scaling error after trimming is neglected)	
		if (firstMap == null) {
			// Previous path must be null
			currentScore = secondTrimmedMap.mappedscore;
			currentScore -= vmProcessor.calcLocalPenaltyFront(secondTrimmedMap);
		}
		else
			if (previousPath != null)
				if (secondMap == null) {
					currentScore = bestPrevScore + firstTrimmedMap.mappedscore - firstMap.mappedscore;
					currentScore -= vmProcessor.calcLocalPenaltyBack(firstTrimmedMap);
				}
				else {
					currentScore = bestPrevScore + firstTrimmedMap.mappedscore + secondTrimmedMap.mappedscore - firstMap.mappedscore;
					currentScore -= vmProcessor.calcPenalty(firstTrimmedMap, secondTrimmedMap);
				}
			else
				currentScore = Double.NEGATIVE_INFINITY; 
	}
	public double getScore() {
		return currentScore;
	}
	public OptMapResultNode getFullTrimmedFirstMap(VirtualMapProcessor vmProcessor) {
		if (firstMap == null)
			return null;
		OptMapResultNode result = new OptMapResultNode(firstMap);
		result.trimResult(this.firstTrimmed * -1 * this.firstMap.mappedstrand, vmProcessor.optrefmap);
		result.trimResult(previousPath.secondTrimmed * this.firstMap.mappedstrand, vmProcessor.optrefmap);		
		vmProcessor.calcScore(result);
		return result;
	}
	public boolean isStart() {
		return firstMap == null;
	}
	public boolean isEnd() {
		return secondMap == null;
	}
}