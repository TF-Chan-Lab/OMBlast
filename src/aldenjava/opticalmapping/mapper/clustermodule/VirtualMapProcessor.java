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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import aldenjava.common.SimpleLocation;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;

public class VirtualMapProcessor {

	public final LinkedHashMap<String, DataNode> optrefmap;
	public final int indelPenalty;
	public final int inversionPenalty;
	public final int translocationPenalty;
	public final PathBuilderFilter indelFilter;
	public final PathBuilderFilter inversionFilter;	
	private final int match;
	private final int fpp;
	private final int fnp;
	private final boolean localPenalty; // Penalty for local alignment
	public VirtualMapProcessor(LinkedHashMap<String, DataNode> optrefmap, int indelPenalty, int inversionPenalty, int translocationPenalty, PathBuilderFilter indelFilter, PathBuilderFilter inversionFilter, int match, int fpp, int fnp, boolean localPenalty)
	{
		super();
		this.optrefmap = optrefmap;
		this.indelPenalty = indelPenalty;
		this.inversionPenalty = inversionPenalty;
		this.translocationPenalty = translocationPenalty;
		this.indelFilter = indelFilter;
		this.inversionFilter = inversionFilter;
		this.match = match;
		this.fpp = fpp;
		this.fnp = fnp;
		this.localPenalty = localPenalty;
	}
	
	public int calcBasicPenalty(OptMapResultNode map1, OptMapResultNode map2)
	{
		int penalty;
		if (indelFilter.checkPass(map1, map2))
			penalty = indelPenalty;
		else
			if (inversionFilter.checkPass(map1, map2))
				penalty = inversionPenalty;
			else
				penalty = translocationPenalty;
		return penalty;
	}
	public int calcPenalty(OptMapResultNode map1, OptMapResultNode map2) {
		// map1, map2 should not overlap each other
		// map1 has a smaller subfragstart and subfragstop
		int penalty;
		if (indelFilter.checkPass(map1, map2))
			penalty = indelPenalty;
		else
			if (inversionFilter.checkPass(map1, map2))
				penalty = inversionPenalty;
			else
				penalty = translocationPenalty;

		int insertion = 0;
		int deletion = 0; 
		insertion = Math.min(map2.subfragstart, map2.subfragstop) - Math.max(map1.subfragstart, map1.subfragstop) - 2;
		if (indelFilter.checkPass(map1, map2) || (inversionFilter.checkPass(map1, map2))) // Cannot account for deletion in translocation
			deletion = Math.max(map1.subrefstart, map2.subrefstart) - Math.min(map1.subrefstop, map2.subrefstop) - 2;

		penalty += fpp * insertion + fnp * deletion;
		return penalty;
				

	}

	public double calcScore(OptMapResultNode vmm)
	{
//		if (vmm instanceof VirtualResultNode)
//			return calcVirtualScore((VirtualResultNode) vmm);
		vmm.updateScore(optrefmap, match, fpp, fnp);
		return vmm.mappedscore;
	}
//	public double calcScore(OptMapResultNode... groupedResult)
//	{
//		return calcScore(Arrays.asList(groupedResult));
//	}
	public double calcScore(List<OptMapResultNode> groupedResult)
	{
		double score = 0;
		OptMapResultNode prevResult = null;
		for (OptMapResultNode result : groupedResult)
		{
			score += calcScore(result);
			if (prevResult != null)
				score -= calcPenalty(prevResult, result);
			prevResult = result;
		}
		if (localPenalty)
		{
			score -= this.calcLocalPenaltyFront(groupedResult.get(0));
			score -= this.calcLocalPenaltyBack(groupedResult.get(groupedResult.size() - 1));
		}
		return score;
	}
//	public double calcVirtualScore(VirtualResultNode vmm)
//	{
//		double totalScore = 0;
//		for (OptMapResultNode map : vmm.mapList)
//		{
//			if (map.cigar == null)
//				totalScore += map.mappedscore;
//			else
//			{
//				map.updateScore(optrefmap, match, fpp, fnp);
//				totalScore += map.mappedscore;
//			}
//		}
//		OptMapResultNode prevMap = null;
//		double penalty = 0;
//		for (OptMapResultNode map : vmm.mapList)
//		{
//			// Give penalty to every SVs
//			if (prevMap != null)
//				penalty += this.calcPenalty(prevMap, map);
//			prevMap = map;
//		}
//		totalScore -= penalty;
//		
//		if (localPenalty)
//		{
//			totalScore -= this.calcLocalPenaltyFront(vmm.mapList.get(0));
//			totalScore -= this.calcLocalPenaltyBack(vmm.mapList.get(vmm.mapList.size() - 1));
//		}
//		return totalScore;
//	}

	
	public int calcLocalPenalty(OptMapResultNode... maparray) {
		return calcLocalPenalty(Arrays.asList(maparray));
	}
	public int calcLocalPenalty(List<OptMapResultNode> mapList)
	{
		if (localPenalty)
			if (!mapList.isEmpty())
			{
				List<SimpleLocation> slist = new ArrayList<SimpleLocation>();
				for (OptMapResultNode map : mapList)
					slist.add(new SimpleLocation(map.subfragstart, map.subfragstop));
				slist = SimpleLocation.mergeLocationList(slist);
				int extra = 0;
				int start = 1;
				for (SimpleLocation loc : slist)
				{
					extra += loc.min - start;
					start = loc.max + 1;
				}
				
				extra += mapList.get(0).parentFrag.getTotalSegment() - 1 - start;
				//temp
	//			extra = 0;
	//			extra += slist.get(0).min - 1;
	//			extra += vmm.mapList.get(0).parentFrag.getTotalSegment() - 1 - slist.get(slist.size() - 1).max - 1;
				return extra * fpp / 2;
			}
			else
			{
				System.err.println("mapList is empty when calculating virtual penalty");
				return -1; 
			}
		else
			return 0;

	}
	/*
	public int calcLocalPenalty(VirtualResultNode vmm)
	{
		if (localPenalty)
			if (!vmm.mapList.isEmpty())
			{
				return calcLocalPenalty(vmm.mapList);
			}
			else
			{
				int extra = vmm.parentFrag.getTotalSegment() - 1 - 1;
				return extra * fpp / 2; 
			}
		else
			return 0;
	}
	 */
	
	public int calcLocalPenaltyFront(OptMapResultNode map) 
	{
		if (!localPenalty)
			return 0;
		int fragpos;
		if (map.mappedstrand == 1)
			fragpos = map.subfragstart - 1;
		else
			fragpos = map.subfragstop - 1;
		DataNode ref = optrefmap.get(map.mappedRegion.ref);
		int insertion = fragpos;
		long displacement = map.length(1, fragpos);
		int index = ref.findRefpIndex((long) (map.mappedRegion.start - displacement));
		assert index >= 0 && index <= ref.refp.length;
		int deletion = map.subrefstart - index;
		// Partial maps tend to go to the start of reference to eliminate the deletions, we have to add extra penalty
		if (index == 0)
		{
			long extralen = 0;
			if (index + 1 <= map.subrefstart - 1)
				extralen = displacement - ref.length(index + 1, map.subrefstart - 1);
			int xpos = map.parentFrag.findRefpIndex(extralen);
			xpos--;
			if (xpos >= 0) // Virtual deletion, equals to extra insertions
				deletion += (xpos + 1) * 2;
		}
			
		return (insertion * fpp + deletion * fnp) / 2;
	}
	public int calcLocalPenaltyBack(OptMapResultNode map) 
	{		
		if (!localPenalty)
			return 0;
		int fragpos;
		if (map.mappedstrand == 1)
			fragpos = map.subfragstop + 1;
		else
			fragpos = map.subfragstart + 1;
		DataNode ref = optrefmap.get(map.mappedRegion.ref);
		int insertion = map.getTotalSegment() - 1 - 1 - fragpos + 1;
		long displacement = map.length(fragpos, map.getTotalSegment() - 1 - 1);
		long refcor;
		if (fragpos >= map.getTotalSegment() - 1)
			refcor = map.mappedRegion.stop; 
		else
			refcor = (long) (map.mappedRegion.stop + displacement);
		int index = ref.findRefpIndex(refcor);
		assert index >= 0 && index <= ref.refp.length;
		int deletion;
		if (index == ref.refp.length)
		{
			deletion = (index - 1) - (map.subrefstop + 1) + 1;
			long extralen = 0;
			if (index + 1 <= map.subrefstart - 1)
				if (map.mappedstrand == 1)
					extralen = displacement - ref.length(map.subrefstop + 1, index - 1);
			long startCor = map.parentFrag.size - extralen + 1;
			int xpos = map.parentFrag.findRefpIndex(startCor);
			deletion += (map.parentFrag.getTotalSegment() - 1 - xpos + 1) * 2;
			
			// Unknown if bug-free
		}
		else
		{
			deletion = (index - (ref.refp[index]==refcor?0:1)) - (map.subrefstop + 1) + 1;
		}
		return (insertion * fpp + deletion * fnp) / 2;
	}
	
}