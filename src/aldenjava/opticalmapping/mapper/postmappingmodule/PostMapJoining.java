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


package aldenjava.opticalmapping.mapper.postmappingmodule;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

import joptsimple.OptionSet;
import aldenjava.common.SimpleLocation;
import aldenjava.opticalmapping.Cigar;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;
import aldenjava.opticalmapping.mapper.ExtensionResult;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;
import aldenjava.opticalmapping.miscellaneous.SelectableMode;

/**
 * This includes several post-mapping functions 1. score filtering 2. Overlapped regions merging 3. Close regions joining (Redundant Alignment Removal)
 * 
 * @author Alden
 */
public class PostMapJoining implements SelectableMode {
	private int postjoinmode = 0;
	private int matchscore = 5;
	private int falseppenalty = 2;
	private int falsenpenalty = 2;
	private boolean allowLocalAlignment = false;
	private LinkedHashMap<String, DataNode> optrefmap;

	public PostMapJoining(LinkedHashMap<String, DataNode> optrefmap) {
		this.optrefmap = optrefmap;
	}

	@Override
	public void setMode(OptionSet options) {
		setMode((int) options.valueOf("postjoinmode"));
	}

	@Override
	public void setMode(int mode) {
		this.postjoinmode = mode;
	}

	@Override
	public int getMode() {
		return this.postjoinmode;
	}

	public void setParameters(OptionSet options) {
		this.setParameters((int) options.valueOf("match"), (int) options.valueOf("fpp"), (int) options.valueOf("fnp"));
	}

	public void setParameters(int matchscore, int falseppenalty, int falsenpenalty) {
		this.matchscore = matchscore;
		this.falseppenalty = falseppenalty;
		this.falsenpenalty = falsenpenalty;
		PostJoinPathNode.matchscore = matchscore;
		PostJoinPathEdge.falseppenalty = falseppenalty;
		PostJoinPathEdge.falsenpenalty = falsenpenalty;
//		this.reflimit = reflimit;
//		this.fraglimit = fraglimit;
	}

	public List<OptMapResultNode> join(List<OptMapResultNode> resultList) {
		switch (postjoinmode) {
			case 0:
				return resultList;
			case 1:
				return removeRepeat(resultList);
			case 2:
				return joinFragmentMap(resultList);
			default:
				return resultList;
		}

	}

	public List<OptMapResultNode> removeRepeat(List<OptMapResultNode> resultList) {
		// remove the repeat result from the resultList
		if (resultList == null)
			return null;
		List<OptMapResultNode> newResultList = new ArrayList<OptMapResultNode>();
		Collections.sort(resultList, OptMapResultNode.mappedstartcomparator);

		List<OptMapResultNode> removeList = new ArrayList<OptMapResultNode>();
		for (int i = 0; i < resultList.size(); i++) {
			OptMapResultNode result = resultList.get(i);
			boolean removed = false;
			for (int j = i + 1; j < resultList.size(); j++) {
				OptMapResultNode target = resultList.get(j);
				// Escape the loop as will never find any equal result
				if (target.subrefstart > result.subrefstart)
					break;
				if (target.parentFrag.name.equals(result.parentFrag.name) && target.mappedRegion.ref.equals(result.mappedRegion.ref) && target.subfragstart == result.subfragstart
						&& target.subfragstop == result.subfragstop && target.subrefstart == result.subrefstart && target.subrefstop == result.subrefstop && target.cigar.equals(result.cigar)) {
					removeList.add(result);
					removed = true;
					break;
				}
			}
			if (!removed)
				newResultList.add(result);
		}
		return newResultList;
	}

	public List<OptMapResultNode> joinFragmentMap(List<OptMapResultNode> fragmentmaplist) {
		if (fragmentmaplist.size() <= 1)
			return new ArrayList<OptMapResultNode>(fragmentmaplist);

		List<OptMapResultNode> forwardfragmentmaplist = new ArrayList<OptMapResultNode>();
		List<OptMapResultNode> reversefragmentmaplist = new ArrayList<OptMapResultNode>();
		for (OptMapResultNode fragmentmap : fragmentmaplist)
			if (fragmentmap.mappedstrand == 1)
				forwardfragmentmaplist.add(fragmentmap);
			else if (fragmentmap.mappedstrand == -1)
				reversefragmentmaplist.add(fragmentmap);
		List<OptMapResultNode> finallist = new ArrayList<OptMapResultNode>();
		finallist.addAll(joinFragmentMap(forwardfragmentmaplist, 1));
		finallist.addAll(joinFragmentMap(reversefragmentmaplist, -1));
		for (OptMapResultNode result : finallist)
			result.updateScore(optrefmap, matchscore, falseppenalty, falsenpenalty);
		return finallist;
	}

	private List<OptMapResultNode> joinFragmentMap(List<OptMapResultNode> fragmentmaplist, int direction) {
		if (fragmentmaplist.size() <= 1)
			return new ArrayList<OptMapResultNode>(fragmentmaplist);
		Collections.sort(fragmentmaplist, OptMapResultNode.mappedstartcomparator);
		List<List<OptMapResultNode>> overlaplistlist = extractOverlapList(fragmentmaplist);
		List<OptMapResultNode> finallist = new ArrayList<OptMapResultNode>();
		for (List<OptMapResultNode> overlaplist : overlaplistlist) {
			if (overlaplist.size() > 1) {

				List<List<PostJoinPathNode>> directedGraph = getGraph(overlaplist, direction);
				List<PostJoinPathNode> startNodeList = getStartNode(directedGraph);
				List<PostJoinPathNode> bestStartNodeList = extractBestStartNode(directedGraph, startNodeList);
				List<OptMapResultNode> extractedfragmentmap = getFragmentMapFromGraph(overlaplist, bestStartNodeList, direction);
				finallist.addAll(extractedfragmentmap);
			} else
				finallist.addAll(overlaplist);
		}
		return finallist;
	}

	private List<List<OptMapResultNode>> extractOverlapList(List<OptMapResultNode> fragmentmaplist) {
		List<List<OptMapResultNode>> overlaplistlist = new ArrayList<List<OptMapResultNode>>();
		List<OptMapResultNode> tmpoverlaplist = new ArrayList<OptMapResultNode>();
		SimpleLocation loc = null;
		String ref = null;
		for (int i = 0; i < fragmentmaplist.size(); i++) {
			OptMapResultNode recentfragmentmap = fragmentmaplist.get(i);

			SimpleLocation recentloc = new SimpleLocation(recentfragmentmap.subrefstart, recentfragmentmap.subrefstop);
			if (loc == null) {
				ref = recentfragmentmap.mappedRegion.ref;
				loc = recentloc;
				tmpoverlaplist = new ArrayList<OptMapResultNode>();
				tmpoverlaplist.add(recentfragmentmap);
			} else {
				if (recentfragmentmap.mappedRegion.ref.equalsIgnoreCase(ref) && recentloc.overlap(loc)) {
					if (recentloc.max > loc.max)
						loc.max = recentloc.max;
					tmpoverlaplist.add(recentfragmentmap);
				} else {
					overlaplistlist.add(tmpoverlaplist);
					tmpoverlaplist = new ArrayList<OptMapResultNode>();
					tmpoverlaplist.add(recentfragmentmap);
					ref = recentfragmentmap.mappedRegion.ref;
					loc = recentloc;
				}
			}
		}
		if (!tmpoverlaplist.isEmpty())
			overlaplistlist.add(tmpoverlaplist);
		return overlaplistlist;

	}

	private List<List<PostJoinPathNode>> getGraph(List<OptMapResultNode> overlaplist, int direction) {
		int start = overlaplist.get(0).subrefstart - 1;
		int stop = 0;
		for (OptMapResultNode fragmentmap : overlaplist)
			if (fragmentmap.subrefstop > stop)
				stop = fragmentmap.subrefstop;
		List<List<PostJoinPathNode>> directedGraph = new ArrayList<List<PostJoinPathNode>>();
		for (int refpos = start; refpos <= stop; refpos++)
			directedGraph.add(new ArrayList<PostJoinPathNode>());
		for (OptMapResultNode fragmentmap : overlaplist) {
			int refpos = fragmentmap.subrefstart - 1;
			int fragpos = fragmentmap.subfragstart - direction;
			String precigar = fragmentmap.cigar.getPrecigar();
			if (!allowLocalAlignment)
				precigar = precigar.substring(precigar.indexOf('M'), precigar.lastIndexOf('M') + 1); // since global alignment can start
			char[] cigararray = precigar.toCharArray();
			PostJoinPathNode savedpjpn = null;
			StringBuilder localprecigar = new StringBuilder();

			for (int j = 0; j < cigararray.length; j++) {
				char c = cigararray[j];
				if (c == 'M') {
					PostJoinPathNode recentpjpn = null;
					for (PostJoinPathNode pjpn : directedGraph.get(refpos - start))
						if (pjpn.fragpos == fragpos) {
							recentpjpn = pjpn;
							break;
						}
					if (recentpjpn == null) {
						recentpjpn = new PostJoinPathNode(refpos, fragpos);
						directedGraph.get(refpos - start).add(recentpjpn);
					}

					recentpjpn.support++;
					if (j != 0)
						recentpjpn.isStart = false;
					if (savedpjpn != null) {
						PostJoinPathEdge edge = new PostJoinPathEdge(localprecigar.toString(), recentpjpn);
						savedpjpn.addEdge(edge); // repeat check done in the module
					}
					savedpjpn = recentpjpn;
					localprecigar = new StringBuilder();

					fragpos += direction;
					refpos++;
				} else {
					if (c == 'I') {
						localprecigar.append('I');
						fragpos += direction;
					} else if (c == 'D') {
						localprecigar.append('D');
						refpos++;
					}
				}

			}
		}
		return directedGraph;
	}

	private List<PostJoinPathNode> getStartNode(List<List<PostJoinPathNode>> directedGraph) {
		List<PostJoinPathNode> startNodeList = new ArrayList<PostJoinPathNode>();
		for (List<PostJoinPathNode> pjpnlist : directedGraph)
			for (PostJoinPathNode pjpn : pjpnlist)
				if (pjpn.isStart)
					startNodeList.add(pjpn);

		return startNodeList;
	}

	private List<PostJoinPathNode> extractBestStartNode(List<List<PostJoinPathNode>> directedGraph, List<PostJoinPathNode> startNodeList) {
		List<Double> startNodeBestScore = new ArrayList<Double>();
		LinkedHashMap<Integer, List<Integer>> passmap = new LinkedHashMap<Integer, List<Integer>>(); // contains the startNode being passed
		for (int j = 0; j < startNodeList.size(); j++) {
			PostJoinPathNode startNode = startNodeList.get(j);
			List<Integer> passlist = startNode.process(j);
			passmap.put(j, new ArrayList<Integer>());
			for (int pass : passlist)
				if (pass != j) {

					passmap.get(j).add(pass);
					passmap.get(pass).add(j);
				}
			startNodeBestScore.add(startNode.bestscore);
		}

		List<Boolean> visitlist = new ArrayList<Boolean>();
		for (int j = 0; j < startNodeList.size(); j++)
			visitlist.add(false);

		List<Integer> bestclist = new ArrayList<Integer>();
		List<Double> bestscorelist = new ArrayList<Double>();

		for (int j = 0; j < startNodeList.size(); j++)
			if (!visitlist.get(j)) {
				List<Integer> samecluster = new ArrayList<Integer>();
				samecluster.add(j);
				samecluster.addAll(passmap.get(j));
				double bestscore = Double.NEGATIVE_INFINITY;
				int bestscorec = -1;
				for (int c : samecluster) {
					visitlist.set(c, true);
					if (startNodeBestScore.get(c) > bestscore) {
						bestscorec = c;
						bestscore = startNodeBestScore.get(c);
					}
				}
				bestclist.add(bestscorec);
				bestscorelist.add(bestscore);
			}
		List<PostJoinPathNode> finalStartNodeList = new ArrayList<PostJoinPathNode>();
		for (int j : bestclist)
			finalStartNodeList.add(startNodeList.get(j));
		return finalStartNodeList;
	}

	private List<OptMapResultNode> getFragmentMapFromGraph(List<OptMapResultNode> overlaplist, List<PostJoinPathNode> startNodeList, int direction) {
		List<OptMapResultNode> fragmentmaplist = new ArrayList<OptMapResultNode>();
		for (PostJoinPathNode pjpn : startNodeList)
			fragmentmaplist.add(getFragmentMapFromGraph(overlaplist, pjpn, direction));
		return fragmentmaplist;
	}

	private OptMapResultNode getFragmentMapFromGraph(List<OptMapResultNode> overlaplist, PostJoinPathNode startNode, int direction) {
		OptMapResultNode fragmentmapsource = overlaplist.get(0);
		StringBuilder newprecigar = new StringBuilder();
		int subfragstart = startNode.fragpos + direction;
		int subrefstart = startNode.refpos + 1;
		PostJoinPathNode pjpn = startNode;
		PostJoinPathEdge edge;
		newprecigar.append('M');
		while ((edge = pjpn.bestedge) != PostJoinPathEdge.endEdge) {
			newprecigar.append(edge.localprecigar);
			newprecigar.append('M');
			pjpn = pjpn.bestedge.nextNode;
		}
		int subfragstop = pjpn.fragpos;
		int subrefstop = pjpn.refpos;
		if (!allowLocalAlignment) // regenerate left / right indel after joining at global mapping
		{
			Cigar leftcigar = null;
			Cigar rightcigar = null;
			for (OptMapResultNode fragmentmap : overlaplist) {
				if (fragmentmap.subrefstart == subrefstart && fragmentmap.subfragstart == subfragstart) {
					String oldprecigar = fragmentmap.cigar.getPrecigar();
					Cigar recentCigar = new Cigar(oldprecigar.substring(0, oldprecigar.indexOf('M')));
					if (leftcigar == null || recentCigar.calcScore(matchscore, falseppenalty, falsenpenalty) >= leftcigar.calcScore(matchscore, falseppenalty, falsenpenalty)) {
						leftcigar = recentCigar;
					}
				}
				if (fragmentmap.subrefstop == subrefstop && fragmentmap.subfragstop == subfragstop)
				{
					String oldprecigar = fragmentmap.cigar.getPrecigar();
					Cigar recentCigar = new Cigar(oldprecigar.substring(oldprecigar.lastIndexOf('M') + 1, oldprecigar.length()));
					if (rightcigar == null || recentCigar.calcScore(matchscore, falseppenalty, falsenpenalty) >= rightcigar.calcScore(matchscore, falseppenalty, falsenpenalty))
						rightcigar = recentCigar;
				}
			}
			newprecigar = new StringBuilder((leftcigar == null ? "" : leftcigar.getPrecigar()) + newprecigar.toString() + (rightcigar == null ? "" : rightcigar.getPrecigar()));
		}
		DataNode fragment = fragmentmapsource.parentFrag;
		double refinedratio = Math.abs(fragment.length(subfragstart, subfragstop)) / (double) optrefmap.get(fragmentmapsource.mappedRegion.ref).length(subrefstart, subrefstop);
		double score = (new Cigar(newprecigar.toString())).calcScore(matchscore, falseppenalty, falsenpenalty);

		ExtensionResult ext;
		if (fragmentmapsource.mappedstrand == 1)
			ext = new ExtensionResult(fragmentmapsource.mappedRegion.ref, subrefstart, subfragstart, subrefstop, subfragstop, newprecigar.toString(), score, refinedratio);
		else
			ext = new ExtensionResult(fragmentmapsource.mappedRegion.ref, subrefstart, fragment.getTotalSegment() - subfragstart - 1, subrefstop, fragment.getTotalSegment() - subfragstop - 1,
					newprecigar.toString(), score, refinedratio);
		return ext.toAlignment(fragment, optrefmap, fragmentmapsource.mappedstrand);
	}

	public PostMapJoining copy() {
		PostMapJoining pmj = new PostMapJoining(optrefmap);
		pmj.setMode(postjoinmode);
		pmj.setParameters(matchscore, falseppenalty, falsenpenalty);
		return pmj;
	}

	public static void assignOptions(ExtendOptionParser parser, int level) {
		parser.addHeader("Post-Join Module Options", level);
		parser.accepts("postjoinmode", "Post-joing mode: 0: no postjoin; 1: removing repeat; 2: full post-joining").withRequiredArg().ofType(Integer.class).defaultsTo(2);
		parser.accepts("match", "Score for one label match").withOptionalArg().ofType(Integer.class).defaultsTo(5);
		parser.accepts("fpp", "False Positive Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
		parser.accepts("fnp", "False Negative Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
	}
}