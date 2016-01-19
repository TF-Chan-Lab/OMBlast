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
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import joptsimple.OptionSet;
import aldenjava.common.SimpleLongLocation;
import aldenjava.opticalmapping.GenomicPosNode;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;
import aldenjava.opticalmapping.miscellaneous.SelectableMode;


public class ResultClusterModule implements SelectableMode{
 // replacing the original ResultCluster class
	
	private LinkedHashMap<String, DataNode> optrefmap;
	
	private int clustermode = 0; 
	
	private long closeReference = 250000;
	private long closeFragment = 250000;
	private int minMatch;
	private int maxTrim;

	private double trimear;
	
	private int match;
	private int fpp;
	private int fnp;
	
	private int indelPenalty;
	private int inversionPenalty;
	private int translocationPenalty;
	private boolean localPenalty = false;
	private PathBuilderFilter indelFilter;
	private PathBuilderFilter inversionFilter;
	
	
	private int minClusterScore = 0;
	private double minconf = 0;
	private double minClusterFragRatio = 0;
	private double minClusterMapSigRatio = 0;
	private int confItemCount = 10;
	private boolean overlapCluster = true;
	private int maxClusterItem = -1;
		
	
	public ResultClusterModule(LinkedHashMap<String, DataNode> optrefmap)
	{
		this.optrefmap = optrefmap;
//		anchorAligner = new AnchorAligner(); 
	}
	@Override
	public void setMode(OptionSet options)
	{
		this.setMode((int) options.valueOf("clustermode"));
	}
	@Override
	public void setMode(int mode) {
		this.clustermode = mode;
	}
	@Override
	public int getMode() {
		return clustermode;
	}


	public void setParameters(OptionSet options)
	{
		this.setParameters((long) options.valueOf("closeref"),
				(long) options.valueOf("closefrag"),
				(int) options.valueOf("minmatch"),
				(int) options.valueOf("maxtrim"),

				(double) options.valueOf("trimear"),
				
				(int) options.valueOf("match"), 
				(int) options.valueOf("fpp"), 				
				(int) options.valueOf("fnp"),
				
				(int) options.valueOf("indelp"),
				(int) options.valueOf("invp"),
				(int) options.valueOf("transp"),
				(boolean) options.valueOf("clusterlocalpenalty"),

				
				(int) options.valueOf("minclusterscore"),
				(double) options.valueOf("minconf"),
				(double) options.valueOf("minclusterfragratio"),
				(double) options.valueOf("minclustersigratio"),
				(boolean) options.valueOf("overlapcluster"),
				(int) options.valueOf("maxclusteritem")
				);
		
				

	}
	public void setParameters(long closeReference, long closeFragment, int minMatch, int maxTrim, 
			double trimear, 
			int match, int fpp, int fnp, // For score updating
			int indelPenalty, int inversionPenalty, int translocationPenalty, boolean localPenalty,
			int minClusterScore, double minconf, double minClusterFragRatio, double minClusterMapSigRatio, boolean overlapCluster, int maxClusterItem) // cluster filtering
	{
		this.closeReference = closeReference;
		this.closeFragment = closeFragment;
		this.indelFilter = new PathBuilderFilter(true, false, true, closeReference, closeFragment);
		this.inversionFilter = new PathBuilderFilter(false, true, false, closeReference, closeFragment); // Note that this inversion filter is different from pbFilter

		this.minMatch = minMatch;
		this.maxTrim = maxTrim;
		
		this.trimear = trimear;
//		this.anchorAligner.setParameters(ear, meas, true);
		
		this.match = match;
		this.fpp = fpp;
		this.fnp = fnp;
		
		this.indelPenalty = indelPenalty;
		this.inversionPenalty = inversionPenalty;
		this.translocationPenalty = translocationPenalty;
		this.localPenalty = localPenalty;
		
		this.minClusterScore = minClusterScore;
		this.minconf = minconf;
		this.minClusterFragRatio = minClusterFragRatio;
		this.minClusterMapSigRatio = minClusterMapSigRatio;
		this.overlapCluster = overlapCluster;
		this.maxClusterItem = maxClusterItem;
	}
	
	/**
	 * Group the partial maps
	 * 
	 * 
	 * @param mapList
	 * @param strandSpecific
	 * @param closeReference
	 * @param closeFragment
	 * @return
	 */
	public List<List<OptMapResultNode>> group(List<OptMapResultNode> mapList, boolean strandSpecific, long closeReference, long closeFragment)
	{
		if (strandSpecific)	{
			List<OptMapResultNode> forwardList = new ArrayList<OptMapResultNode>();
			List<OptMapResultNode> reverseList = new ArrayList<OptMapResultNode>();
			List<List<OptMapResultNode>> unknownListList = new ArrayList<List<OptMapResultNode>>(); // for virtual fragment
			for (OptMapResultNode map : mapList) {
				if (map.mappedstrand == 1)
					forwardList.add(map);
				else 
					if (map.mappedstrand == -1)
						reverseList.add(map);
					else {
						List<OptMapResultNode> unknownList = new ArrayList<OptMapResultNode>();
						unknownList.add(map);
						unknownListList.add(unknownList);
					}
			}
			List<List<OptMapResultNode>> clusteredList = new ArrayList<List<OptMapResultNode>>();
			clusteredList.addAll(group(forwardList, false, closeReference, closeFragment));
			clusteredList.addAll(group(reverseList, false, closeReference, closeFragment));
			clusteredList.addAll(unknownListList);
			return clusteredList;
		}
		else {
			// first group according to reference
			List<List<OptMapResultNode>> clusteredList = new ArrayList<List<OptMapResultNode>>(); 
			Collections.sort(mapList, OptMapResultNode.mappedstartcomparator);
			List<List<OptMapResultNode>> groupRefList = new ArrayList<List<OptMapResultNode>>();
			if (closeReference < 0)
				clusteredList.add(mapList);
			else {
				GenomicPosNode region = null;
				List<OptMapResultNode> recentList = null;
				for (OptMapResultNode map : mapList) {
					if (region == null || map.mappedRegion.start == -1 || !map.isClose(region, closeReference)) {
						recentList = new ArrayList<OptMapResultNode>();
						groupRefList.add(recentList);
						if (map.mappedRegion.start == -1)
							region = null;
						else
							region = new GenomicPosNode(map.mappedRegion.ref, map.mappedRegion.start, map.mappedRegion.stop);
					}
					recentList.add(map);
					if (map.mappedRegion.stop > region.stop)
						region = new GenomicPosNode(region.ref, region.start, map.mappedRegion.stop);
				}
				clusteredList = groupRefList;
			}
			
			if (closeFragment < 0)
				return clusteredList;
			else {
				List<List<OptMapResultNode>> groupFragList = new ArrayList<List<OptMapResultNode>>();
				for (List<OptMapResultNode> tmpMapList : clusteredList)	{
					Collections.sort(tmpMapList, OptMapResultNode.subfragstartstopcomparator);
					SimpleLongLocation loc = null;
					List<OptMapResultNode> recentList = null;
					for (OptMapResultNode map : tmpMapList)	{
						SimpleLongLocation maploc = new SimpleLongLocation(map.length(0, map.subfragstart - 1), map.length(0, map.subfragstop + 1));
						if (loc == null || (!loc.overlap(maploc, closeFragment)))
						{
							recentList = new ArrayList<OptMapResultNode>();
							groupFragList.add(recentList);
							loc = maploc;
						}
						recentList.add(map);
						if (maploc.max > loc.max)
							loc.max = maploc.max;
						if (maploc.min < loc.min)
							loc.min = maploc.min;
					}
				}
				clusteredList = groupFragList;

				return clusteredList;
			}
			
		}
	}
	
	
	/**
	 * Check if two partial maps can be connected.
	 * Trimming is done on overlapping partial maps
	 * 
	 * It is essential to use trimming.
	 * Yet, there are two drawbacks using trimming
	 * 1. Trimming is expensive.
	 * 2. Trimming may lead to inappropriate results (This is done by validation using trimear) 
	 * 
	 * 
	 * @param optrefmap
	 * @param maxTrim
	 * @param result1
	 * @param result2
	 * @param vmProcessor
	 * @param ear
	 * @return
	 */
	private TrimResult trimOverlap(LinkedHashMap<String, DataNode> optrefmap, OptMapResultNode[] trim1Result, OptMapResultNode[] trim2Result, PathBuilderFilter pbFilter, VirtualMapProcessor vmProcessor) {
		int penalty = vmProcessor.calcBasicPenalty(trim1Result[0], trim2Result[0]);
		double bestscore = Double.NEGATIVE_INFINITY;
		int bestTrim1 = -1;
		int bestTrim2 = -1;
		boolean[][] trimPair = new boolean[maxTrim + 1][maxTrim + 1];
		// temp increasing maxTrim greatly:
		// 1. if the two score addition is already lower than any one of the score, discard
		// 2. maxTrim
		for (int trimmed = 0; trimmed <= maxTrim; trimmed++)
		{
			for (int trim1 = 0; trim1 <= trimmed; trim1++) {
				int trim2 = trimmed - trim1;
				if (trimPair[trim1][trim2]) // this trim pair should no longer be done
					continue;				
				if ((trim1Result[0].cigar.getMatch() - trim1 >= minMatch) && (trim2Result[0].cigar.getMatch() - trim2 >= minMatch)) {
					OptMapResultNode tresult1 = trim1Result[trim1];
					OptMapResultNode tresult2 = trim2Result[trim2];
					// trim1Result[trim1 - 1] must exist
					if (tresult1 == null) {

						tresult1 = new OptMapResultNode(trim1Result[trim1 - 1]);
						tresult1.trimResult(-1 * trim1Result[trim1 - 1].mappedstrand, optrefmap);
						tresult1.updateScore(optrefmap, match, fpp, fnp);
						trim1Result[trim1] = tresult1;
					}
					if (tresult2 == null) {
						tresult2 = new OptMapResultNode(trim2Result[trim2 - 1]);
						tresult2.trimResult(1 * trim2Result[trim2 - 1].mappedstrand, optrefmap);
						tresult2.updateScore(optrefmap, match, fpp, fnp);
						trim2Result[trim2] = tresult2;
					}
					
					boolean withinTrimear = true;
					withinTrimear = Math.abs(tresult1.getMapScale() - 1) < trimear && Math.abs(tresult2.getMapScale() - 1) < trimear;
					boolean removeLaterTrim = false;
				
					if (pbFilter.checkPass(tresult1, tresult2)) {
						if (withinTrimear && !tresult1.overlap(tresult2)) {
							if (tresult1.mappedscore + tresult2.mappedscore - vmProcessor.calcPenalty(tresult1, tresult2) > bestscore) {
								bestscore = tresult1.mappedscore + tresult2.mappedscore - vmProcessor.calcPenalty(tresult1, tresult2);
								bestTrim1 = trim1;
								bestTrim2 = trim2;
							}
							removeLaterTrim = true;
						}
						else
							if (withinTrimear && tresult1.mappedscore + tresult2.mappedscore - penalty < bestscore) // We could not calcPenalty directly when tresult1.overlap(tresult2)
								removeLaterTrim = true;
					}
					else
						removeLaterTrim = true;
					
					if (removeLaterTrim)
						for (int i = trim1; i <= maxTrim; i++)
							for (int j = trim2; j <= maxTrim; j++)
								trimPair[i][j] = true;
				}
			}
		}

		if (bestTrim1 == -1)
			return TrimResult.FailedTrimResult();
		else
			return TrimResult.SuccessfulTrimResult(trim1Result[bestTrim1], trim2Result[bestTrim2]);
	}
	
	/**
	 * Build paths for two partial maps if they have the following properties
	 * 1. No overlap (passing trimOverlap)
	 * 2. passing PathBuilderFilter
	 * 
	 * @param groupedMap
	 * @param pbFilter
	 * @param vmProcessor
	 * @param maxTrim
	 * @param ear
	 * @return
	 */
	private List<ClusterPathNode> buildPath(List<OptMapResultNode> groupedMap, PathBuilderFilter pbFilter, VirtualMapProcessor vmProcessor) {
		// resulting path will not form loop (i.e. is acyclic), no A->B B->C C->A occurs at any time. But should handle null start and null stop properly
		
		Collections.sort(groupedMap, OptMapResultNode.subfragstartstopcomparator);
		List<ClusterPathNode> clusterPathList = new ArrayList<ClusterPathNode>();
		OptMapResultNode[][] trim1Results = new OptMapResultNode[groupedMap.size()][maxTrim + 1];
		OptMapResultNode[][] trim2Results = new OptMapResultNode[groupedMap.size()][maxTrim + 1];
		for (int i = 0; i < groupedMap.size(); i++) {
			trim1Results[i][0] = groupedMap.get(i);
			trim2Results[i][0] = groupedMap.get(i);
		}
		
		// head of the path
		for (OptMapResultNode map : groupedMap)
			clusterPathList.add(new ClusterPathNode(null, map, 0, 0, null, map));
			
		// start building paths for map to map
		for (int i = 0; i < groupedMap.size(); i++)	{
			int j = i - 1;
			while (j >= 0) {
				if (pbFilter.checkPass(groupedMap.get(j), groupedMap.get(i)))	{
					TrimResult trimResult = trimOverlap(optrefmap, trim1Results[j], trim2Results[i], pbFilter, vmProcessor);
					if  (trimResult.successful) {
						OptMapResultNode map1 = trimResult.result1;
						OptMapResultNode map2 = trimResult.result2;
						int trim1 = groupedMap.get(j).cigar.getMatch() - map1.cigar.getMatch();
						int trim2 = groupedMap.get(i).cigar.getMatch() - map2.cigar.getMatch();
						ClusterPathNode cp = new ClusterPathNode(groupedMap.get(j), groupedMap.get(i), trim1, trim2, map1, map2);
						clusterPathList.add(cp);
					}
				}
				j--;
			}
		}
		// tail of the path
		for (OptMapResultNode map : groupedMap)			
			clusterPathList.add(new ClusterPathNode(map, null, 0, 0, map, null));

		for (ClusterPathNode path : clusterPathList)
			path.updateScore(vmProcessor);
		return clusterPathList;
	}
	
	/**
	 * Links the PathNode and constructs edges in <code>ClusterPathNode - previousPath</code>
	 * @param clusterPathList
	 * @param vmProcessor
	 */
	private void linkPath(List<ClusterPathNode> clusterPathList, VirtualMapProcessor vmProcessor) {
		// important property: all cluster must be sorted by second map.
		
		// Build a map for fast access to the respective ClusterPathNode according to the first map
		Map<OptMapResultNode, List<ClusterPathNode>> mapForFirst = new LinkedHashMap<OptMapResultNode, List<ClusterPathNode>>();
		for (int i = 0; i < clusterPathList.size(); i++) {
			if (!mapForFirst.containsKey(clusterPathList.get(i).firstMap))
				mapForFirst.put(clusterPathList.get(i).firstMap, new ArrayList<ClusterPathNode>());
			mapForFirst.get(clusterPathList.get(i).firstMap).add(clusterPathList.get(i));
		}
			
		for (int i = 0; i < clusterPathList.size(); i++) {
			ClusterPathNode cp1 = clusterPathList.get(i);
			mapForFirst.get(cp1.firstMap).remove(cp1);
			if ((cp1.previousPath == null && cp1.firstMap == null) || (cp1.previousPath != null)) {
				cp1.updateScore(vmProcessor);
				if (mapForFirst.containsKey(cp1.secondMap))
					for (ClusterPathNode cp2 : mapForFirst.get(cp1.secondMap))
						if (cp1.canLink(cp2, minMatch))
							cp2.assignPreviousPath(cp1);
			}
		}
	}
	
	private List<ClusterPathNode> getBestPath(List<ClusterPathNode> clusterPathList) {
		ClusterPathNode bestPathNode = null;
		double bestScore = Double.NEGATIVE_INFINITY;
		
		// Find out the best last ClusterPathNode
		for (ClusterPathNode cp : clusterPathList) {
			if (cp.getScore() >= bestScore && cp.isEnd()) {
				bestPathNode = cp;
				bestScore = cp.getScore();
			}
		}
		
		// Trace back the path
		List<ClusterPathNode> bestPath = new ArrayList<ClusterPathNode>();
		ClusterPathNode recentPath = bestPathNode;
		while (recentPath != null) {
			bestPath.add(recentPath);
			recentPath = recentPath.previousPath;
		}
		Collections.reverse(bestPath);
		
		return bestPath;
	}
	
	private List<OptMapResultNode> convertPathToMapList(List<ClusterPathNode> clusterPathList, VirtualMapProcessor vmProcessor) {
		List<OptMapResultNode> mapList = new ArrayList<OptMapResultNode>();
		for (ClusterPathNode cp : clusterPathList)
			if (!cp.isStart())
				mapList.add(cp.getFullTrimmedFirstMap(vmProcessor));
		return mapList;
	}
	
	/**
	 * A module to check the connection state of two partial maps 
	 * 1. size difference (indel)
	 * 2. too many extra / missing signals
	 * 
	 * If two partial maps connection belongs to case 2, attempt to join the partial maps to single partial map 
	 *  
	 * Module is temporarily unused.  
	 * @param resultList
	 * @param vmProcessor
	 */
	/*
	private void checkDirectLink(List<OptMapResultNode> resultList, VirtualMapProcessor vmProcessor)
	{
		if (true)
			return;
//		int threshold = 1000;
		if (resultList != null)
			for (int i = resultList.size() - 1; i > 0; i--)
			{
				OptMapResultNode result1 = resultList.get(i - 1);
				OptMapResultNode result2 = resultList.get(i);
				if (indelFilter.checkPass(result1, result2))
				{
					
					int strand = result1.mappedstrand;
//					System.out.println("======");
					if (strand == -1)
					{
						result1.reverse();
						result2.reverse();
						OptMapResultNode r = result1;
						result1 = result2;
						result2 = r;
					}
//					long refDis = (result2.mappedRegion.start - result1.mappedRegion.stop + 1);
//					long moleDis = mole.refp[result2.getSubMoleSigStart()] - mole.refp[result1.getSubMoleSigStop()] + 1;
					DataNode ref = optrefmap.get(result1.mappedRegion.ref);
					DataNode mole = result1.parentFrag;
					
					Cigar cigar = anchorAligner.align(ref, mole, result1.getSubRefSigStart(), result2.getSubRefSigStop(), result1.getSubMoleSigStart(), result2.getSubMoleSigStop());
					if (cigar != null)
					{
						OptMapResultNode newResult = new OptMapResultNode(mole, null, result1.mappedstrand, result1.subrefstart, result2.subrefstop, result1.subfragstart, result2.subfragstop, cigar, -1, -1);
						newResult.updateMappedRegion(optrefmap.get(result1.mappedRegion.ref));
						newResult.updateScore(optrefmap, match, fpp, fnp);
						if (newResult.getMatch() >= result1.getMatch() + result2.getMatch())
							if (result1.mappedscore + result2.mappedscore - vmProcessor.calcPenalty(result1, result2) < newResult.mappedscore) 
							{
								resultList.remove(i);
								resultList.remove(i - 1);
								if (strand == -1)
									newResult.reverse();
								resultList.add(i - 1, newResult);
							}
					}
//					moleDis *= (2 / (result1.getMapScale() + result2.getMapScale()));
//					if (Math.abs(refDis - moleDis) <= threshold)
//					{
//						Cigar midc = AnchorMapperT.align(optrefmap.get(result1.mappedRegion.ref), mole, result1.getSubRefSigStop(), result2.getSubRefSigStart(), result1.getSubMoleSigStop(), result2.getSubMoleSigStart(), threshold, false);
//						Cigar c = new Cigar();
//						c.append(result1.cigar);
//						c.append(midc);
//						c.append(result2.cigar);
//						OptMapResultNode newResult = new OptMapResultNode(mole, null, result1.mappedstrand, result1.subrefstart, result2.subrefstop, result1.subfragstart, result2.subfragstop, c, -1, -1);
//						newResult.updateMappedRegion(optrefmap.get(result1.mappedRegion.ref));
//						newResult.updateScore(optrefmap, match, fpp, fnp);
//						resultList.remove(i);
//						resultList.remove(i - 1);
//						if (strand == -1)
//							newResult.reverse();
//						resultList.add(i - 1, newResult);
//					}
					if (strand == -1)
					{
						result1.reverse();
						result2.reverse();
					}
					
				}
			}
//		if (x > 0)
//			System.exit(0);

	}
	*/
	/**
	 * 
	 * Cluster the results by using 
	 * 1. Path builder filter which limits indel, inversion and translocation 
	 * 2. vmProcessor to calculate the score
	 * 
	 * The logic includes two steps
	 * 1. Grouping
	 * 2. Build paths between partial maps (which may undergo trimming)
	 * 3. Link paths to paths (not partial maps to partial maps)
	 * 4. Get the path with best score and convert it back to partial maps 
	 *  
	 * @param mapList
	 * @param pbFilter
	 * @param vmProcessor
	 * @return
	 */
	private List<ClusteredResult> cluster(List<OptMapResultNode> mapList, PathBuilderFilter pbFilter, VirtualMapProcessor vmProcessor)
	{
		// exit if no result is input
		if (mapList == null || mapList.size() == 0)
			return new ArrayList<ClusteredResult>(); 
		
		
		// Pre-update the results
		for (OptMapResultNode result : mapList)
			vmProcessor.calcScore(result); 
		
		
		
		List<ClusteredResult> clusteredResultList = new ArrayList<ClusteredResult>();
		
		// Grouping close results to be 1. clustered effectively 2. output one clustered result per group
		List<List<OptMapResultNode>> groupedMapList = group(mapList, pbFilter.sameStrand, pbFilter.closeReference, pbFilter.closeFragment);
		
		// In each result group
		for (List<OptMapResultNode> groupedMap : groupedMapList) {
			// Build a path whenever the results can be joined
			// Please see ClusterPathNode for details
			List<ClusterPathNode> clusterPathList = buildPath(groupedMap, pbFilter, vmProcessor);
			// Link the paths
			linkPath(clusterPathList, vmProcessor);
			// Get the best path and extract the final result maps
			List<OptMapResultNode> clusteredGroupMap = convertPathToMapList(getBestPath(clusterPathList), vmProcessor);
//			if (false)			
//			{
//				for (ClusterPathNode cp : clusterPathList) {
//					if (cp.currentScore < 0)
//						continue;
//					if (cp.secondMap != null)
//						continue;
//					List<ClusterPathNode> bestPath = new ArrayList<ClusterPathNode>();
//					ClusterPathNode recentPath = cp;
//					while (recentPath != null)
//					{
//						bestPath.add(recentPath);
//						recentPath = recentPath.previousPath;
//					}
//					Collections.reverse(bestPath);
//					if (convertPathToMapList(bestPath).isEmpty()) 
//						continue;
//					System.out.println("PATH - " + cp.currentScore);
//					for (OptMapResultNode result : convertPathToMapList(bestPath))
//						System.out.println(result);
//				}
//			}
			ClusteredResult cr = new ClusteredResult();
			// We want to retain the original alignments
			List<OptMapResultNode> expandedGroupedMap = new ArrayList<OptMapResultNode>(); // Saving the original results
			for (OptMapResultNode map : groupedMap)
				expandedGroupedMap.addAll(map.getRealMap());

			// We want to import the modified alignments
			List<OptMapResultNode> expandedClusterGroupMap = new ArrayList<OptMapResultNode>(); // Saving the updated results
			for (OptMapResultNode map : clusteredGroupMap)
				expandedClusterGroupMap.addAll(map.getRealMap());
			
//			checkDirectLink(expandedClusterGroupMap, vmProcessor);
			cr.importUpdatedResult(expandedClusterGroupMap);
			if (cr.updatedResult.size() > 0) {
				
				cr.process(vmProcessor);
				clusteredResultList.add(cr);
			}
			
		}
		return clusteredResultList;
	}
	
	private void removeOverlap(List<ClusteredResult> clusteredResultList)
	{
		Collections.sort(clusteredResultList);
		Collections.reverse(clusteredResultList);
		for (int i = clusteredResultList.size() - 1; i >= 1; i--) {
			for (int j = i - 1; j >= 0; j--)
				for (OptMapResultNode result1 : clusteredResultList.get(i).updatedResult)
					for (OptMapResultNode result2 : clusteredResultList.get(j).updatedResult)
						if (result1.overlap(result2)) {
							clusteredResultList.remove(i);
							break;
						}
		}
					
	}
	/*
	@Deprecated
	private List<ClusteredResult> multiPhaseCluster(List<OptMapResultNode> mapList, List<PathBuilderFilter> filterList, List<VirtualMapProcessor> vmProcessorList)
	{
		PathBuilderFilter lastFilter = filterList.get(filterList.size() - 1);
		List<OptMapResultNode> currentList = mapList;
		for (int i = 0; i < filterList.size(); i++)
		{
//			System.out.println("Phase " + Integer.toString(i));
			PathBuilderFilter pbFilter = filterList.get(i);
			VirtualMapProcessor vmProcessor = vmProcessorList.get(i);
			List<ClusteredResult> clusteredResult = cluster(currentList, pbFilter, vmProcessor);
			if (pbFilter != lastFilter)
			{
				currentList = new ArrayList<OptMapResultNode>();
				for (ClusteredResult result : clusteredResult)
				{
					VirtualResultNode m = new VirtualResultNode(result.updatedResult, vmProcessor);
					currentList.add(m);
				}
				
			}
			else
				return clusteredResult;
		}
		return null; // no cluster filter applied
	}
	*/
	
	/**
	 * Another confidence scheme
	 */
	/*
	private void tempConfidence(List<ClusteredResult> crList)
	{
		for (ClusteredResult cr : crList)
		{
			
			for (OptMapResultNode result : cr.updatedResult)
			{
				double chisum = 0;
				int direction = result.mappedstrand;
				int refpos = result.subrefstart - 1;
				int fragpos = result.subfragstart - direction;
				int lastfragpos = -1;
				int lastrefpos = -1;
				int matchedRegion = 0;
				for (char c : result.cigar.getPrecigar().toCharArray())
				{
					switch (c)
					{
						 case 'M':						 
							 
							 if (lastrefpos != -1)
							 {
								 
								 long fraglen;
								 if (direction == 1) 
									 fraglen = result.parentFrag.length(lastfragpos, fragpos);
								 else
									 fraglen = result.parentFrag.length(fragpos, lastfragpos);
								 long reflen = optrefmap.get(result.mappedRegion.ref).length(lastrefpos, refpos);
								 
								 double sumlen = (fraglen + reflen) * (fraglen + reflen);
								 double difflen = (fraglen - reflen) * (fraglen - reflen);
								 double prodlen = (fraglen * reflen);
								 chisum += 0.2 * difflen / sumlen 
										 + 0.04 * difflen / prodlen 
										 + 0.002 * difflen / reflen; 
								 matchedRegion++;
							 }
							 lastfragpos = fragpos;
							 lastrefpos = refpos;
							 fragpos += direction;
							 refpos++;
							 break;						 
						 case 'I':
							 fragpos += direction;
							 break;
						 case 'D':
							 refpos++;
							 break;
					}
				}
				ChiSquaredDistribution d = new ChiSquaredDistribution(Math.max(matchedRegion, 1));
				System.out.println(chisum);
						
				result.confidence = d.cumulativeProbability(chisum);
//				result.confidence = Math.max(result.confidence, 1 - result.confidence);
			}
		}
	}
	*/
	/**
	 * Process the confidence, a measure of uniqueness of score
	 * @param crList
	 */
	private void processConfidence(List<ClusteredResult> crList) {
		int itemCount = confItemCount;
		if (itemCount == -1 || itemCount > crList.size())
			itemCount = crList.size();
			
		double minscore = Double.MAX_VALUE;
		double bestscore = Double.MIN_VALUE;
		for (int i = 0; i < itemCount; i++)
		{
			ClusteredResult cr = crList.get(i);
			if (cr.score < minscore)
				minscore = cr.score;
			if (cr.score > bestscore)
				bestscore = cr.score;
		}
		int base = 1;
		double total = 0;
		if (minscore == bestscore)
			// evenly distributed, rare case
			for (int i = 0; i < itemCount; i++)
			{
				ClusteredResult cr = crList.get(i);
				for (OptMapResultNode result : cr.updatedResult)
					result.confidence = (1 / (double) crList.size());
			}
		else
		{
			for (int i = 0; i < itemCount; i++)
			{
				ClusteredResult cr = crList.get(i);
				total += (cr.score - minscore) + base;
			}
			
			for (int i = 0; i < itemCount; i++)
			{
				ClusteredResult cr = crList.get(i);
				for (OptMapResultNode result : cr.updatedResult)
					result.confidence = ((cr.score - minscore + base) / (double) total);
			}
		}

		// Confidence not calculated within itemCount range is set as 0
		for (int i = itemCount; i < crList.size(); i++)
		{
			ClusteredResult cr = crList.get(i);
			for (OptMapResultNode result : cr.updatedResult)
				result.confidence = 0;
		}
		

	}

	private List<ClusteredResult> filterClusteredResult(List<ClusteredResult> crList) {
		List<ClusteredResult> filteredList = new ArrayList<ClusteredResult>();
		for (ClusteredResult cr : crList)
			if (!cr.updatedResult.isEmpty())
				if (cr.updatedResult.get(0).confidence >= minconf
						&& cr.score >= minClusterScore 
						&& (minClusterFragRatio == -1 || cr.getFragRatio() >= minClusterFragRatio)
						&& (minClusterMapSigRatio == -1 || cr.getMapSigRatio() >= minClusterMapSigRatio))
					filteredList.add(cr);
		return filteredList;
	}
	
	private List<ClusteredResult> standardindelcluster(List<OptMapResultNode> mapList) {
		PathBuilderFilter pbFilter1 = new PathBuilderFilter(true, false, true, closeReference, closeFragment);
		VirtualMapProcessor vmProcessor1 = new VirtualMapProcessor(optrefmap, indelPenalty, inversionPenalty, translocationPenalty, indelFilter, inversionFilter, match, fpp, fnp, localPenalty);
		
		return cluster(mapList, pbFilter1, vmProcessor1);
	}
	private List<ClusteredResult> standardinversioncluster(List<OptMapResultNode> mapList) {
		PathBuilderFilter pbFilter1 = new PathBuilderFilter(false, false, false, closeReference, closeFragment); // Both can be used
		VirtualMapProcessor vmProcessor1 = new VirtualMapProcessor(optrefmap, indelPenalty, inversionPenalty, translocationPenalty, indelFilter, inversionFilter, match, fpp, fnp, localPenalty);
		
		return cluster(mapList, pbFilter1, vmProcessor1);
	}
	private List<ClusteredResult> standardtransloccluster(List<OptMapResultNode> mapList) {		
		PathBuilderFilter pbFilter1 = new PathBuilderFilter(false, false, false, -1, closeFragment);
		VirtualMapProcessor vmProcessor1 = new VirtualMapProcessor(optrefmap, indelPenalty, inversionPenalty, translocationPenalty, indelFilter, inversionFilter, match, fpp, fnp, localPenalty);
		
		return cluster(mapList, pbFilter1, vmProcessor1);
	}
	private List<ClusteredResult> standardnocluster(List<OptMapResultNode> mapList)	{
		List<ClusteredResult> clusteredResultList = new ArrayList<ClusteredResult>();
//		VirtualMapProcessor vmProcessor1 = new VirtualMapProcessor(optrefmap, indelPenalty, inversionPenalty, translocationPenalty, indelFilter, inversionFilter, match, fpp, fnp, localPenalty);
		
		for (OptMapResultNode result : mapList) {
			List<OptMapResultNode> rmapList = new ArrayList<OptMapResultNode>();
			rmapList.add(result);
			ClusteredResult cr = new ClusteredResult();
//			cr.importOriginalResult(rmapList);
			cr.importUpdatedResult(rmapList);
			cr.score = result.mappedscore;
//			cr.process();
//			if (cr.score < mins)
//				mins = cr.score;

			clusteredResultList.add(cr);
		}
		return clusteredResultList;
		
	}
	
	public List<ClusteredResult> standardcluster(List<OptMapResultNode> mapList) {
		return standardcluster(mapList, false);
	}
	public List<ClusteredResult> standardcluster(List<OptMapResultNode> mapList, boolean processConfidence)	{
		if (mapList.isEmpty())
			return new ArrayList<ClusteredResult>();
		List<ClusteredResult> clusteredResultList;
		switch (clustermode) {
			case 0: clusteredResultList = standardnocluster(mapList); break;
			case 1: clusteredResultList = standardindelcluster(mapList); break;
			case 2: clusteredResultList = standardinversioncluster(mapList); break;
			case 3: clusteredResultList = standardtransloccluster(mapList); break;
			default: clusteredResultList = standardnocluster(mapList);
		}
		Collections.sort(clusteredResultList);
		Collections.reverse(clusteredResultList);
		if (processConfidence)
			processConfidence(clusteredResultList);
		clusteredResultList = filterClusteredResult(clusteredResultList);
		if (!overlapCluster)
			this.removeOverlap(clusteredResultList);
		if (maxClusterItem != -1)
			while (clusteredResultList.size() > maxClusterItem)
				clusteredResultList.remove(clusteredResultList.size() - 1);
		return clusteredResultList;
	}

	public static void assignOptions(ExtendOptionParser parser, int level)
	{
		parser.addHeader("Result Clustering Options", level);
		parser.accepts("clustermode", "The clustering mode. 0: No clustering. 1: Standard indel cluster. 2: Standard indel-inv cluster. 3. Standard cluster transloc").withRequiredArg().ofType(Integer.TYPE).defaultsTo(0);
		parser.accepts("closeref", "The max distance (reference) between two results to be considered at same cluster.").withOptionalArg().ofType(Long.class).defaultsTo((long) 250000);
		parser.accepts("closefrag", "The max distance (fragment) between two results to be considered at same cluster.").withOptionalArg().ofType(Long.class).defaultsTo((long) 250000);
		parser.accepts("minmatch", "Min match to be consider as a valid partial map.").withOptionalArg().ofType(Integer.class).defaultsTo(3);
		parser.accepts("maxtrim", "Maximum trim of a partial map. ").withOptionalArg().ofType(Integer.class).defaultsTo(5);

		parser.accepts("trimear", "Error acceptable range.").withOptionalArg().ofType(Double.class).defaultsTo(0.1);
		
		parser.accepts("match", "Score for one label match").withOptionalArg().ofType(Integer.class).defaultsTo(5);
		parser.accepts("fpp", "False Positive Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
		parser.accepts("fnp", "False Negative Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
		
		parser.accepts("indelp", "Indel Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(10);
		parser.accepts("invp", "Inversion Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(30);
		parser.accepts("transp", "Translocation Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(50);
		parser.accepts("clusterlocalpenalty", "Local Penalty after clustering").withOptionalArg().ofType(Boolean.class).defaultsTo(false);
		
		parser.accepts("minclusterscore", "Min cluster score.").withOptionalArg().ofType(Integer.class).defaultsTo(30);
		parser.accepts("minconf", "Minimum confidence").withOptionalArg().ofType(Double.class).defaultsTo(0.4);
		parser.accepts("minclusterfragratio", "Min clustered fragment ratio").withOptionalArg().ofType(Double.class).defaultsTo(-1.0);
		parser.accepts("minclustersigratio", "Min clustered mapped signal ratio").withOptionalArg().ofType(Double.class).defaultsTo(-1.0);
		parser.accepts("overlapcluster", "Allow overlapping clusters at results.").withOptionalArg().ofType(Boolean.class).defaultsTo(true);
		parser.accepts("maxclusteritem", "Maximum items output. -1: no limit").withOptionalArg().ofType(Integer.class).defaultsTo(1);
	}

	public ResultClusterModule copy() {
		ResultClusterModule rcm = new ResultClusterModule(optrefmap);
		rcm.setMode(clustermode);
		rcm.setParameters(closeReference, closeFragment, minMatch, maxTrim, trimear, match, fpp, fnp, indelPenalty, inversionPenalty, translocationPenalty, localPenalty, minClusterScore, minconf, minClusterFragRatio, minClusterMapSigRatio, overlapCluster, maxClusterItem);
		return rcm;
	}

	
}

class TrimResult {
	public final boolean successful;
	public final OptMapResultNode result1;
	public final OptMapResultNode result2;
	private TrimResult(boolean successful, OptMapResultNode result1, OptMapResultNode result2) {
		super();
		this.successful = successful;
		this.result1 = result1;
		this.result2 = result2;
	}
	public static TrimResult SuccessfulTrimResult(OptMapResultNode result1, OptMapResultNode result2) {
		return new TrimResult(true, result1, result2);
	}
	public static TrimResult FailedTrimResult() {
		return new TrimResult(false, null, null);
	}
	
}