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
import java.util.LinkedHashMap;
import java.util.List;

import joptsimple.OptionSet;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;
import aldenjava.opticalmapping.miscellaneous.Copyable;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;
import aldenjava.opticalmapping.miscellaneous.SelectableMode;

/**
 * Filters alignment results
 * 
 * @author Alden
 *
 */
public class Filter implements Copyable<Filter>, SelectableMode {

	// Filter mode
	private int filtermode;

	// Normal filter
	private int minmatch;
	private int maxfp;
	private int maxfn;
	private double maxfpr;
	private double maxfnr;
	private double minscore;
	private double minsubfragratio;
	private double minsigratio;

	// filter with trimming
	private LinkedHashMap<String, DataNode> optrefmap;
	private int trimmode = 0;
	private int maxTrim = 0;
	private int match = 5;
	private int fpp = 2;
	private int fnp = 2;

	/**
	 * Constructs a new <code>Filter</code> with reference information provided 
	 * @param optrefmap the reference information
	 */
	public Filter(LinkedHashMap<String, DataNode> optrefmap) {
	}

	/**
	 * Constructs a new <code>Filter</code> with no reference information 
	 */
	public Filter() {
		this(null);
	}

	@Override
	public void setMode(OptionSet options) {
		this.setMode((int) options.valueOf("filtermode"));
	}

	@Override
	public void setMode(int mode) {
		this.filtermode = mode;
	}

	@Override
	public int getMode() {
		return filtermode;
	}

	public void setParameters(OptionSet options) {
		this.setParameters((int) options.valueOf("minmatch"), (int) options.valueOf("maxfp"), (int) options.valueOf("maxfn"), (double) options.valueOf("maxfpr"), (double) options.valueOf("maxfnr"),
				(double) options.valueOf("minscore"), (double) options.valueOf("minsubfragratio"), (double) options.valueOf("minsigratio"), (int) options.valueOf("trimmode"),
				(int) options.valueOf("maxtrim"), (int) options.valueOf("match"), (int) options.valueOf("fpp"), (int) options.valueOf("fnp"));

	}

	public void setParameters(int minmatch, int maxfp, int maxfn, double maxfpr, double maxfnr, double minscore, double minsubfragratio, double minsigratio, int trimmode, int maxTrim, int match,
			int fpp, int fnp) {
		this.minmatch = minmatch;
		this.maxfp = maxfp;
		this.maxfn = maxfn;
		this.maxfpr = maxfpr;
		this.maxfnr = maxfnr;
		this.minscore = minscore;
		this.minsubfragratio = minsubfragratio;
		this.minsigratio = minsigratio;
		this.trimmode = trimmode;
		this.maxTrim = maxTrim;
		this.match = match;
		this.fpp = fpp;
		this.fnp = fnp;
	}

	/**
	 * Checks if a result passes the filtering criteria
	 * 
	 * @param result
	 *            the result to be checked
	 * @return <code>true</code> if the result passes the filtering criteria
	 */
	private boolean checkPass(OptMapResultNode result) {
		return (result.mappedscore >= minscore 
		&& result.getSubFragRatio() >= minsubfragratio 
		&& result.getMatch() >= minmatch // this match is based on total match of signal. Not used for calculating match score, which equals to (matched signal - 1) * matchscore
		&& result.getFP() <= maxfp 
		&& result.getFN() <= maxfn 
		&& result.getFPRate() <= maxfpr 
		&& result.getFNRate() <= maxfnr 
		&& result.getMapSigRatio() >= minsigratio);
	}

	/**
	 * Checks if the result passes some selected filtering criteria which cannot be rescued even by trimming.
	 * 
	 * @param result 
	 * 				the result to be checked
	 * @return <code>true</code> if the result passes the selected filtering criteria
	 */
	private boolean checkPassChanceWithTrim(OptMapResultNode result) {
		return (result.mappedscore >= minscore 
				&& result.getSubFragRatio() >= minsubfragratio 
				&& result.getMatch() >= minmatch 
				&& result.getMapSigRatio() >= minsigratio);
	}

	/**
	 * Trims the result so that it can pass the filtering criteria. 
	 * 
	 * @param result
	 *              the result to be trimmed and checked
	 * @return A new trimmed result or <code>null</code> if the result cannot pass the filtering criteria after trimming
	 */
	private OptMapResultNode trimToPassFilter(OptMapResultNode result) {
		if (optrefmap == null)
			; // Should do something on it like throwing an exception
		
		int maxTrim = this.maxTrim;
		if (maxTrim > result.getMatch() - minmatch)
			maxTrim = result.getMatch() - minmatch;
		if (maxTrim < 0)
			return null;
		boolean[][] trimPair = new boolean[maxTrim + 1][maxTrim + 1];
		OptMapResultNode finalMap = null;
		for (int trim = 0; trim <= maxTrim; trim++) {
			for (int leftTrim = 0; leftTrim <= trim; leftTrim++) {
				int rightTrim = trim - leftTrim;
				if (trimPair[leftTrim][rightTrim])
					continue;
				OptMapResultNode map = new OptMapResultNode(result);
				map.trimResult(leftTrim, optrefmap);
				map.trimResult(rightTrim * -1, optrefmap);
				map.updateMappedRegion(optrefmap.get(map.mappedRegion.ref));
				map.updateScore(optrefmap, match, fpp, fnp);
				if (checkPass(map)) {
					if (finalMap == null || finalMap.mappedscore < map.mappedscore)
						finalMap = map;
					// No need to trim further
					for (int i = leftTrim + 1; i <= maxTrim; i++)
						for (int j = rightTrim + 1; j <= maxTrim; j++)
							trimPair[i][j] = true;
				}
				if (!checkPassChanceWithTrim(map)) {
					// No more chance
					for (int i = leftTrim + 1; i <= maxTrim; i++)
						for (int j = rightTrim + 1; j <= maxTrim; j++)
							trimPair[i][j] = true;
				}
			}
		}
		return finalMap;
	}
	/**
	 * Checks if a result passes the score filtering criterion. Other criteria are neglected. 
	 * 
	 * @param result
	 * @return <code>true</code> if the result passes the score filtering criterion
	 */
	private boolean checkPassScore(OptMapResultNode result) {
		return result.mappedscore >= minscore;
	}
	
	/**
	 * Filters the results. If trimming is used, the result instance in the retained list is not the original result instance
	 * 
	 * @param resultlist
	 *              a list of result to be filtered
	 * @return A retained list of results or an empty list if all results are filtered out
	 */
	public List<OptMapResultNode> filter(List<OptMapResultNode> resultlist) {
		List<OptMapResultNode> newlist = new ArrayList<OptMapResultNode>();
		for (OptMapResultNode result : resultlist) {
			if (result.isUsed())
				switch (filtermode) {
					case 0:
						newlist.add(result);
						break;
					case 1:
						switch (trimmode) {
							case 0:
								if (checkPass(result))
									newlist.add(result);
								break;
							case 1:
								OptMapResultNode newresult = trimToPassFilter(result);
								if (newresult != null)
									newlist.add(newresult);
								break;
							default:
								throw new IllegalStateException("Selected trim mode does not exist");
						}
						break;
					case 2:
						if (checkPassScore(result))
							newlist.add(result);
						break;
					default:
						throw new IllegalArgumentException("Selected filter mode does not exist");
				}
		}
		return newlist;
	}

	/**
	 * Restricts the filtering to score filtering only.
	 * 
	 * @return
	 */
	public boolean restrictToScoreOnly() {
		if (filtermode == 1) {
			filtermode = 2;
			trimmode = 0;
			return true;
		}
		return false;
	}

	@Override
	public Filter copy() {
		Filter filter = new Filter(optrefmap);
		filter.setMode(filtermode);
		filter.setParameters(minmatch, maxfp, maxfn, maxfpr, maxfnr, minscore, minsubfragratio, minsigratio, trimmode, maxTrim, match, fpp, fnp);
		return filter;
	}

	public static void assignOptions(ExtendOptionParser parser, int level) {
		parser.addHeader("Result Filter Options", level);
		parser.accepts("filtermode", "Filter Mode. 0: No filter; 1: Normal Filter; 2: Filter Score Only").withOptionalArg().ofType(Integer.class).defaultsTo(0);
		parser.accepts("minmatch", "Minumum matches of a map").withOptionalArg().ofType(Integer.class).defaultsTo(3);
		parser.accepts("maxfp", "Maximum false positives of a map").withOptionalArg().ofType(Integer.class).defaultsTo(10000);
		parser.accepts("maxfn", "Maximum false negatives of a map").withOptionalArg().ofType(Integer.class).defaultsTo(10000);
		parser.accepts("maxfpr", "Maximum false positive rate of a map").withOptionalArg().ofType(Double.class).defaultsTo(1e-4);
		parser.accepts("maxfnr", "Maximum false negative rate of a map").withOptionalArg().ofType(Double.class).defaultsTo(0.5);
		parser.accepts("minscore", "Minimum score of a map").withOptionalArg().ofType(Double.class).defaultsTo(0.0);
		parser.accepts("minsubfragratio", "Minimum subfragment ratio of a map").withOptionalArg().ofType(Double.class).defaultsTo(0.0);
		parser.accepts("minsigratio", "Minimum mapped signal ratio of a map").withOptionalArg().ofType(Double.class).defaultsTo(0.0);

		parser.accepts("trimmode", "Trim Mode. 0: Trim mode disabled; 1: Trim mode enabled").withOptionalArg().ofType(Integer.class).defaultsTo(0);
		parser.accepts("maxtrim", "Maximum trim of a partial map. ").withOptionalArg().ofType(Integer.class).defaultsTo(5);
		parser.accepts("match", "Score for one label match").withOptionalArg().ofType(Integer.class).defaultsTo(5);
		parser.accepts("fpp", "False Positive Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
		parser.accepts("fnp", "False Negative Penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
	}

}
