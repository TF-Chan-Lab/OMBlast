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
import java.util.List;

public class PostJoinPathNode
{
	// the following should be unused, now I and D are represented in edge
	// refpos = -1 fragpos != -1: I
	// fragpos = -1 refpos != -1: D
	// refpos != -1 and fragpos != -1 M
	// refpos = -1 fragpos = -1 N
	
	public static int matchscore = 5;
	
	
	
	public final int refpos;
	public final int fragpos;
	public int support;
	public List<PostJoinPathEdge> nextEdge;
	public boolean isStart; // not a child of any node
	
	// for use in graph
//	public int pass; // being passthrough, useful in analyzing the graph
	public List<Integer> passlist;
	public double bestscore;
	public PostJoinPathEdge bestedge;
	public PostJoinPathNode(int refpos, int fragpos)
	{
		this.refpos = refpos;
		this.fragpos = fragpos;
		this.support = 0;
		this.nextEdge = new ArrayList<PostJoinPathEdge>();
		
		this.isStart = true;
		this.passlist = new ArrayList<Integer>();
		this.bestscore = matchscore;
		this.bestedge = null;
	}
	public void addEdge(PostJoinPathEdge edge)
	{
		boolean added = false;
		for (PostJoinPathEdge pastedge : nextEdge)
			if (pastedge.equals(edge))
				added = true;
		if (!added)
			nextEdge.add(edge);
	}
	public List<Integer> process(int pass)
	{
		if (!passlist.contains(pass))
			passlist.add(pass);
		if (bestedge != null)
			return passlist; // don't reprocess
		// find best edge
	
//		System.out.println("PathNode process...");
//		System.out.print(refpos);
//		System.out.print(" ");
//		System.out.println(fragpos);
		for (PostJoinPathEdge edge : nextEdge)
		{
			List<Integer> edgepasslist = edge.process(pass);
			for (int edgepass : edgepasslist)
				if (!passlist.contains(edgepass))
					passlist.add(edgepass);
			
			if (edge.getScore() + matchscore > bestscore)
			{
				this.bestscore = edge.getScore() + matchscore;
				this.bestedge = edge;
			}
		}
//		System.out.println("PathNode process ends.");
//		System.out.print(refpos);
//		System.out.print(" ");
//		System.out.println(fragpos);
//		System.out.println();
//		
		
		if (nextEdge.isEmpty() || bestedge == null) // In case no further edge can be extended, or the future edge points to negative score
		{
			this.bestedge = PostJoinPathEdge.endEdge;
			this.bestscore = matchscore;
		}

		return passlist;
	}
}
