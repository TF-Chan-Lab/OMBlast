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


package aldenjava.opticalmapping.mapper.omblastmapper;

import java.util.LinkedHashMap;

import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.mapper.ExtensionResult;
import aldenjava.opticalmapping.mapper.seeding.Kmer;
import aldenjava.opticalmapping.mapper.seeding.Seed;
import aldenjava.opticalmapping.miscellaneous.Copyable;

/**
 * A class for extending the data from a matched seed
 * 
 * @author Alden
 * 
 */
public class SeedExtension implements Copyable<SeedExtension> {
	private final LinkedHashMap<String, DataNode> optrefmap;
	private int measure = 500;
	private int matchscore = 5;
	private int falseppenalty = 2;
	private int falsenpenalty = 2;
	private int falselimit = 5;
	private double ear = 0.05;
	private boolean allowLocalAlignment;

	/**
	 * Constructs a <code>SeedExtension</code> module with the reference information
	 * 
	 * @param optrefmap
	 *            the reference information
	 */
	public SeedExtension(LinkedHashMap<String, DataNode> optrefmap) {
		this.optrefmap = optrefmap;
	}

	public void setParameters(int measure, double ear, int matchscore, int falseppenalty, int falsenpenalty, int falselimit, boolean allowLocalAlignment) {
		this.measure = measure;
		this.matchscore = matchscore;
		this.falseppenalty = falseppenalty;
		this.falsenpenalty = falsenpenalty;
		this.falselimit = falselimit;
		this.ear = ear;
		this.allowLocalAlignment = allowLocalAlignment;
	}

	/**
	 * Extends the data on reference using a scale <code>ratio</code> from the initial position
	 * 
	 * @param ref
	 *            the reference information
	 * @param data
	 *            alignment of <code>data</code> to be extended according to <code>seed</code>
	 * @param initialrefpos
	 *            the initial reference position
	 * @param initialdatapos
	 *            the initial data position
	 * @param direction
	 *            progressing steps
	 * @param scale
	 *            the scale for extension
	 * @return the result of extension
	 */
	public ExtensionResult extendCore(DataNode ref, DataNode data, int initialrefpos, int initialdatapos, int direction, double scale) {
		int refstartpos = -1;
		int fragmentstartpos = -1;
		int score = matchscore; // Now any matching signals count, not any
								// matching subfragment counts
		int highestscore = 0;
		int refpos = initialrefpos;
		int datapos = initialdatapos;
		int err = 0;
		double cumulatefragmentlen = 0;
		double cumulatereflen = 0;

		StringBuilder precigar = new StringBuilder();
		StringBuilder tmpcigar = new StringBuilder();
		do {
			if (Math.abs(cumulatereflen - cumulatefragmentlen) <= measure) {
				score += matchscore;
				tmpcigar.append("M");
				// double tshift = cumulatefragmentlen - cumulatereflen;
				// longer is allowed!
				if (score > highestscore || !allowLocalAlignment) // go ahead!! Onlyhigher score and longer is allowed!
				{
					refstartpos = refpos;
					fragmentstartpos = datapos;
					highestscore = score;
					precigar.append(tmpcigar);
					tmpcigar = new StringBuilder();
				}
				datapos += direction;
				refpos += direction;
				cumulatefragmentlen += data.getRefl(datapos) + 1;
				cumulatereflen += (ref.getRefl(refpos) + 1) * scale;
				err = 0;
			} else {
				if (cumulatereflen > cumulatefragmentlen) // falsep
				{
					err++;
					datapos += direction;
					score -= falseppenalty;
					cumulatefragmentlen += data.getRefl(datapos) + 1;
					tmpcigar.append("I");
				} else // falsen
				{
					err++;
					refpos += direction;
					score -= falsenpenalty;
					cumulatereflen += (ref.getRefl(refpos) + 1) * scale;
					tmpcigar.append("D");
				}
				if (!allowLocalAlignment) {
					highestscore = score;
					precigar.append(tmpcigar);
					tmpcigar = new StringBuilder();
				}
			}

		} while ((datapos >= 1 && datapos < data.getTotalSegment() - 1) && (refpos >= 1 && refpos < ref.refp.length) && (err <= falselimit || !allowLocalAlignment));

		// Need to resolve wrong cigar problem:
		// Direction -->
		// --R----------R----R
		// --S-------S--S----S
		// W: M-------M--I----M
		// C: M-------I--M----M
		// to minimizing measuring error,
		// DP:
		// cumulatefragmentlen = 0;
		//
		//
		// int refstart = initialrefpos;
		// int refstop = refstartpos;
		// int fragstart = initialfragmentpos;
		// int fragstop = fragmentstartpos;
		//
		// List<List<Integer>> fragToRefList = new ArrayList<List<Integer>>();
		//
		// for (int f = fragstart; f <= fragstop; f++)
		// {
		// cumulatefragmentlen += fragment.fragment[f] + 1;
		// cumulatereflen = 0;
		// fragToRefList.add(new ArrayList<Integer>());
		// for (int r = refstart; r <= refstop; r++)
		// {
		// cumulatereflen += (ref.refl[refpos] + 1) * ratio;
		// if (Math.abs(cumulatereflen - cumulatefragmentlen) <= measure)
		// fragToRefList.get(f - fragstart).add(r);
		// }
		// }
		//
		// class MatchPair implements Comparable<MatchPair>
		// {
		// int lastref;
		// int lastfrag;
		// MatchPair prev;
		// double meas;
		// int fn;
		// int fp;
		// public MatchPair(int lastref, int lastfrag, MatchPair prev, double
		// meas, int fn, int fp) {
		// this.lastref = lastref;
		// this.lastfrag = lastfrag;
		// this.prev = prev;
		// this.meas = meas;
		// }
		// public MatchPair(MatchPair m)
		// {
		// this.lastref = m.lastref;
		// this.lastfrag = m.lastfrag;
		// this.prev = m.prev;
		// this.meas = m.meas;
		// }
		// public int compareTo(MatchPair m) {
		// int result = 0;
		// result = Integer.compare(this.lastfrag, m.lastfrag);
		// return Integer.compare(this.lastref, m.lastref);
		//
		// }
		//
		// }
		// MatchPair initial = new MatchPair(-1, -1, null, 0, 0, 0);
		// List<MatchPair> matchPairList = new ArrayList<MatchPair>();
		// for (int f = fragstart; f <= fragstop; f++)
		// {
		// List<MatchPair> newMatchPairList = new ArrayList<MatchPair>();
		// for (MatchPair m : matchPairList)
		// {
		// MatchPair fpmatch = new MatchPair(m);
		// fpmatch.fp++;
		// newMatchPairList.add(fpmatch);
		//
		//
		// }
		//
		// }
		//

		// if (!allowLocalAlignment)
		// {
		// if (direction == -1)
		// while (fragpos >= 1)
		// {
		// tmpcigar.append("I");
		// score -= falseppenalty;
		// highestscore = score;
		// fragmentstartpos = fragpos;
		// fragpos += direction;
		//
		// }
		// if (direction == 1)
		// while (fragpos < fragment.getTotalSegment() - 1)
		// {
		// tmpcigar.append("I");
		// score -= falseppenalty;
		// highestscore = score;
		// fragmentstartpos = fragpos;
		// fragpos += direction;
		// }
		// precigar.append(tmpcigar);
		// }
		if (direction == -1)
			precigar.reverse();
		return new ExtensionResult(ref.name, initialrefpos, initialdatapos, refstartpos, fragmentstartpos, precigar.toString(), highestscore, scale);

	}

	/**
	 * Attempts to find a better scaling factor by extending with 3 scaling factors at <code>startscale - ear</code>, <code>startscale</code> and <code>startscale + ear</code>. The scaling factor which yields the best extension results is used for next round of <code>extensionLoop</code> (recursion) with smaller <code>ear</code>
	 * 
	 * @param data
	 *            alignment of <code>data</code> to be extended according to <code>seed</code>
	 * @param seed
	 *            as the start point of extension
	 * @param startscale
	 *            the central point of the scale
	 * @param ear
	 *            the flanking range of the scale
	 * @param times
	 *            Number of rounds remained for recursion to get a better scaling factor
	 * @return
	 */
	public ExtensionResult extensionLoop(DataNode data, Seed seed, double startscale, double ear, int times) {
		Kmer refKmer = seed;
		Kmer dataKmer = seed.kmerpointer;
		DataNode ref = optrefmap.get(refKmer.source);
		double highestscore = 0;
		ExtensionResult combinedExtension = null;
		double highestratio = -1;
		for (double ratio = startscale - ear; ratio <= startscale + ear; ratio += ear / 2) {
			ExtensionResult leftExtension = this.extendCore(ref, data, refKmer.pos, dataKmer.pos, -1, ratio);
			ExtensionResult rightExtension = this.extendCore(ref, data, refKmer.pos + refKmer.k() + refKmer.getErrorNo() - 1, dataKmer.pos + dataKmer.k() + dataKmer.getErrorNo() - 1, 1, ratio);
			// module for refining scale at global alignment; no longer used as we now define the starting position to be the first matched signal
			double refinedratio = data.length(leftExtension.stopfinalfragmentpos, rightExtension.stopfinalfragmentpos)
					/ (double) ref.length(leftExtension.stopfinalrefpos, rightExtension.stopfinalrefpos);

			double finalscore = (leftExtension.score + rightExtension.score + seed.getCigar(false).calcScore(matchscore, falseppenalty, falsenpenalty)) * (1 - Math.abs((1 - refinedratio)));
			if (combinedExtension == null
					|| finalscore > highestscore
					|| (finalscore == highestscore && (rightExtension.stopfinalfragmentpos - leftExtension.stopfinalfragmentpos + 1) > (combinedExtension.stopfinalfragmentpos
							- combinedExtension.startfinalfragmentpos + 1))) {
				StringBuilder finalprecigar = new StringBuilder();
				finalprecigar.append(leftExtension.precigar);

				finalprecigar.append(seed.getCigar(false).getPrecigar());

				finalprecigar.append(rightExtension.precigar);

				ExtensionResult tmpExtension = new ExtensionResult(ref.name, leftExtension.stopfinalrefpos, leftExtension.stopfinalfragmentpos, rightExtension.stopfinalrefpos,
						rightExtension.stopfinalfragmentpos, finalprecigar.toString(), finalscore, refinedratio);

				highestscore = finalscore;
				combinedExtension = tmpExtension;
				highestratio = ratio;
			}
			if (ear == 0)
				break;
		}
		if (times == 1)
			return combinedExtension;
		else
			return extensionLoop(data, seed, highestratio, ear / 2, times - 1);

	}

	/**
	 * Extends the <code>data</code> according to the given <code>seed</code>. Scaling range bound is set according to the seed information. Note that the boundary is not a strict boundary but an initial guess for the scaling factor
	 * 
	 * @param data
	 *            alignment of <code>data</code> to be extended according to <code>seed</code>
	 * @param seed
	 *            as the start point of extension
	 * @param rangeUBound
	 *            scaling range upper bound
	 * @param rangeLBound
	 *            scaling range lower bound
	 * @return the result of extension
	 */
	public ExtensionResult extension(DataNode data, Seed seed, double rangeUBound, double rangeLBound) {
		if (seed.rangeUBound == -1 || seed.rangeLBound == -1)
			return extensionLoop(data, seed, 1, ear, 3);
		else
			return extensionLoop(data, seed, (rangeUBound + rangeLBound) / 2, (rangeUBound - rangeLBound) / 2, 3);
	}

	/**
	 * Extends the <code>data</code> according to the given <code>seed</code>. Scaling range bound is set according to the seed information
	 * 
	 * @param data
	 *            <code>data</code> to be extended
	 * @param seed
	 *            <code>seed</code> as the start point of extension
	 * @return the result of extension
	 * @see #extension(DataNode, Seed, double, double)
	 */
	public ExtensionResult extension(DataNode data, Seed seed) {
		return extension(data, seed, seed.rangeUBound, seed.rangeLBound);
	}

	@Override
	public SeedExtension copy() {
		SeedExtension newse = new SeedExtension(optrefmap);
		newse.setParameters(measure, ear, matchscore, falseppenalty, falsenpenalty, falselimit, allowLocalAlignment);
		return newse;
	}
}
