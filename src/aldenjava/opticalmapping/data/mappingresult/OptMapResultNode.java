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


package aldenjava.opticalmapping.data.mappingresult;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;

import aldenjava.common.SimpleLocation;
import aldenjava.common.SimpleLongLocation;
import aldenjava.opticalmapping.Cigar;
import aldenjava.opticalmapping.GenomicPosNode;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.data.SimulationInfo;

public class OptMapResultNode {
	// Parameters
	public DataNode parentFrag;
	public GenomicPosNode mappedRegion;
	public int mappedstrand;
	public int subrefstart;
	public int subrefstop;
	public int subfragstart;
	public int subfragstop;
	public Cigar cigar;
	public double mappedscore;
	public double confidence;

	// refp indicator
	// Due to inefficient initialization, this is implemented but not used in constructor.
	public List<MatchingSignalPair> mspsQ;
	public List<MatchingSignalPair> mspsR;

	// Init
	public OptMapResultNode(DataNode f, GenomicPosNode mappedRegion, int mappedstrand, int subrefstart, int subrefstop, int subfragstart, int subfragstop, Cigar cigar, double mappedscore,
			double confidence) {
		parentFrag = f;
		this.mappedRegion = mappedRegion;
		this.mappedstrand = mappedstrand;
		this.subrefstart = subrefstart;
		this.subrefstop = subrefstop;
		this.subfragstart = subfragstart;
		this.subfragstop = subfragstop;
		this.cigar = cigar;
		this.mappedscore = mappedscore;
		this.confidence = confidence;
		// this.updateMSP();
	}

	public OptMapResultNode(OptMapResultNode r) {
		parentFrag = r.parentFrag;
		this.mappedRegion = r.mappedRegion;
		this.mappedstrand = r.mappedstrand;
		this.subrefstart = r.subrefstart;
		this.subrefstop = r.subrefstop;
		this.subfragstart = r.subfragstart;
		this.subfragstop = r.subfragstop;
		if (r.cigar != null)
			this.cigar = new Cigar(r.cigar);
		else
			this.cigar = null;
		this.mappedscore = r.mappedscore;
		this.confidence = r.confidence;
		// this.updateMSP();
	}

	public void setCigar(Cigar cigar) {
		this.cigar = cigar;
		// this.updateMSP();
	}

	@Deprecated
	public int[] getRelativeRefPos() {
		if (cigar == null)
			return null;
		if (parentFrag == null)
			return null;
		int[] relativeRefPos = new int[this.parentFrag.getTotalSegment()];
		for (int i = 0; i < relativeRefPos.length; i++)
			relativeRefPos[i] = -1;
		boolean lastMatched = false;
		int fragpos = subfragstart - mappedstrand;
		int refpos = subrefstart - 1;
		for (char c : cigar.getPrecigar().toCharArray()) {
			switch (c) {
				case 'M':
					if (lastMatched) {
						relativeRefPos[fragpos] = refpos;
					}
					refpos++;
					fragpos += mappedstrand;
					lastMatched = true;
					break;
				case 'I':
					// relative RefPos is still -1
					fragpos += mappedstrand;
					lastMatched = false;
					break;
				case 'D':
					refpos++;
					lastMatched = false;
					break;
				default:
					System.err.println("Warning!");
					break;
			}
		}
		return relativeRefPos;

	}

	public List<OptMapResultNode> getBreakResult(LinkedHashMap<String, DataNode> optrefmap, int meas, double ear) {
		if (parentFrag == null || optrefmap == null)
			return null;
		List<OptMapResultNode> resultList = new ArrayList<OptMapResultNode>();
		if (!isUsed()) {
			resultList.add(new OptMapResultNode(this));
			return resultList;
		}

		DataNode ref = optrefmap.get(mappedRegion.ref);
		String precigar = this.cigar.getPrecigar();
		int lastrefpos = -1;
		int lastfragpos = -1;
		int currentrefpos = subrefstart;
		int currentfragpos = subfragstart;
		int lastrefstart = -1;
		int lastfragstart = -1;
		for (char c : precigar.toCharArray())
			switch (c) {
				case 'M':
					if (lastrefpos != -1) {
						if (lastrefstart == -1) {
							lastrefstart = lastrefpos;
							lastfragstart = lastfragpos;
						}
						// start comparison
						long reflen = ref.length(lastrefpos, currentrefpos - 1);
						long fraglen;
						if (mappedstrand == 1)
							fraglen = parentFrag.length(lastfragpos, currentfragpos - mappedstrand);
						else if (mappedstrand == -1)
							fraglen = parentFrag.length(currentfragpos - mappedstrand, lastfragpos);
						else
							fraglen = 0;

						boolean pass = (reflen * (1 - ear) - meas < fraglen && fraglen < reflen * (1 + ear) + meas);
						if (!pass) {
							if (lastfragstart != lastfragpos) {
								resultList.add(this.getSubResult(ref, lastfragstart, lastfragpos - mappedstrand)); // the current one should not be included

							}
							lastrefstart = -1;
							lastfragstart = -1;
						} else {
							if (pass && currentrefpos == subrefstop + 1) {
								resultList.add(this.getSubResult(ref, lastfragstart, subfragstop)); // the current one should not be included
								lastrefstart = -1;
								lastfragstart = -1;
							}
						}
					}
					lastrefpos = currentrefpos;
					lastfragpos = currentfragpos;
					currentfragpos += mappedstrand;
					currentrefpos += 1;
					break;
				case 'I':
					currentfragpos += mappedstrand;
					break;
				case 'D':
					currentrefpos += 1;
					break;
				default:
					;
			}
		return resultList;
	}

	public List<OptMapResultNode> getBreakResult(LinkedHashMap<String, DataNode> optrefmap, GenomicPosNode unwantedRegion) {
		if (parentFrag == null || optrefmap == null)
			return null;
		List<OptMapResultNode> resultList = new ArrayList<OptMapResultNode>();
		if (!isUsed()) {
			resultList.add(new OptMapResultNode(this));
			return resultList;
		}

		DataNode ref = optrefmap.get(mappedRegion.ref);
		String precigar = this.cigar.getPrecigar();
		int lastrefpos = -1;
		int lastfragpos = -1;
		int currentrefpos = subrefstart;
		int currentfragpos = subfragstart;
		int lastrefstart = -1;
		int lastfragstart = -1;
		for (char c : precigar.toCharArray()) {
			switch (c) {
				case 'M':
					if (lastrefpos != -1) {
						if (lastrefstart == -1) {
							lastrefstart = lastrefpos;
							lastfragstart = lastfragpos;
						}
						GenomicPosNode currentRegion = new GenomicPosNode(ref.name, ref.refp[lastrefstart - 1], ref.refp[currentrefpos - 1]);

						boolean pass = !(unwantedRegion.overlapSize(currentRegion) > 0);
						if (!pass) {
							if (lastfragstart != lastfragpos)
								resultList.add(this.getSubResult(ref, lastfragstart, lastfragpos - mappedstrand));

							lastrefstart = -1;
							lastfragstart = -1;
						} else {
							if (pass && currentrefpos == subrefstop + 1) {
								resultList.add(this.getSubResult(ref, lastfragstart, subfragstop));
								lastrefstart = -1;
								lastfragstart = -1;
							}
						}
					}
					lastrefpos = currentrefpos;
					lastfragpos = currentfragpos;
					currentfragpos += mappedstrand;
					currentrefpos += 1;
					break;
				case 'I':
					currentfragpos += mappedstrand;
					break;
				case 'D':
					currentrefpos += 1;
					break;
				default:
					;
			}
		}
		return resultList;

	}

	public int getFP() {
		if (cigar != null)
			return cigar.getFP();
		else
			return -1;
	}

	public int getFN() {
		if (cigar != null)
			return cigar.getFN();
		else
			return -1;
	}

	public double getFPRate() {
		if (cigar == null)
			return -1;
		return this.getFP() / (double) (Math.abs(length(subfragstart, subfragstop)) + 1);
	}

	public double getFNRate() {
		if (cigar == null)
			return -1;
		return this.getFN() / (double) (Math.abs(subrefstop - subrefstart) + 1);
	}

	public long getMapLength() {
		if (subfragstart == -1 || subfragstop == -1)
			return -1;
		if (mappedstrand == 1)
			return parentFrag.length(subfragstart, subfragstop);
		else if (mappedstrand == -1)
			return parentFrag.length(subfragstop, subfragstart);
		else
			return 0;
	}

	public double getMapScale() {
		if (this.isUsed())
			return (getMapLength() / (double) (mappedRegion.length()));
		else
			return -1;
	}

	public double getMapSigRatio() {
		return (cigar.getMatch() / (double) (this.parentFrag.getTotalSegment() - 1));
	}

	public int[] getRefMapSignal() {
		if (cigar == null)
			return null;
		if (parentFrag == null)
			return null;

		int[] matchSignals = new int[cigar.getMatch()];
		int index = 0;
		int fragpos = subfragstart - mappedstrand;
		int refpos = subrefstart - 1;
		for (char c : cigar.getPrecigar().toCharArray()) {
			switch (c) {
				case 'M':
					matchSignals[index++] = refpos; // need to be checked
					refpos++;
					fragpos += mappedstrand;
					break;
				case 'I':
					// relative RefPos is still -1
					fragpos += mappedstrand;
					break;
				case 'D':
					refpos++;
					break;
				default:
					throw new RuntimeException(cigar + " character");
					// break;
			}
		}
		return matchSignals;
	}

	public int[] getMapSignal() {
		if (cigar == null)
			return null;
		if (parentFrag == null)
			return null;

		int[] matchSignals = new int[cigar.getMatch()];
		int index = 0;
		int fragpos = subfragstart - mappedstrand;
		int refpos = subrefstart - 1;
		for (char c : cigar.getPrecigar().toCharArray()) {
			switch (c) {
				case 'M':
					matchSignals[index++] = fragpos; // need to be checked
					refpos++;
					fragpos += mappedstrand;
					break;
				case 'I':
					// relative RefPos is still -1
					fragpos += mappedstrand;
					break;
				case 'D':
					refpos++;
					break;
				default:
					throw new RuntimeException(cigar + " character");
					// break;
			}
		}
		return matchSignals;

	}

	public int getMatch() {
		if (cigar != null)
			return cigar.getMatch();
		else
			return -1;
	}

	public GenomicPosNode getMoleMappedRegion() {
		// Not-reversed
		if (mappedstrand == 1)
			return new GenomicPosNode(parentFrag.name, parentFrag.refp[subfragstart - 1], parentFrag.refp[subfragstop]);
		else
			return new GenomicPosNode(parentFrag.name, parentFrag.refp[subfragstop - 1], parentFrag.refp[subfragstart]);
	}

	public int getSignal(String refname, int targetSig) {
		if (!mappedRegion.ref.equals(refname))
			return -1;
		int currRefSig = subrefstart - 1;
		int currFragSig = subfragstart;

		if (mappedstrand == 1)
			currFragSig -= 1;
		for (char c : cigar.getPrecigar().toCharArray()) {
			switch (c) {
				case 'M':
					if (targetSig == currRefSig) {
						return currFragSig;
					}
					currRefSig++;
					currFragSig += mappedstrand;
					break;
				case 'I':
					currFragSig += mappedstrand;
					break;
				case 'D':
					currRefSig++;
					break;
				default:
					break;
			}
		}
		return -1;
	}

	public long getSignalDisplace(String refname, int targetSig) {
		int currRefSig = subrefstart - 1;
		int currFragSig = subfragstart;

		if (mappedstrand == 1)
			currFragSig -= 1;
		for (char c : cigar.getPrecigar().toCharArray()) {
			switch (c) {
				case 'M':
					if (targetSig == currRefSig) {
						if (mappedstrand == 1)
							if (subfragstart > currFragSig)
								return 0;
							else
								return length(subfragstart, currFragSig);
						else if (currFragSig + 1 > subfragstart)
							return 0;
						else
							return length(currFragSig + 1, subfragstart);
					}
					currRefSig++;
					currFragSig += mappedstrand;
					break;
				case 'I':
					currFragSig += mappedstrand;
					break;
				case 'D':
					currRefSig++;
					break;
				default:
					break;
			}
		}
		return 0;
	}

	public double getSubFragRatio() {
		if (getMapLength() == -1)
			return -1;
		return (getMapLength() / (double) (length(1, this.getTotalSegment() - 2)));
	}

	public OptMapResultNode getSubResult(DataNode ref, int fragstart, int fragstop) {
		if ((mappedstrand == 1 && subfragstart <= fragstart && fragstart <= subfragstop && subfragstart <= fragstop && fragstop <= subfragstop)
				|| (mappedstrand == -1 && subfragstop <= fragstart && fragstart <= subfragstart && subfragstop <= fragstop && fragstop <= subfragstart)) {
			String precigar = this.cigar.getPrecigar();
			int currentrefpos = subrefstart;
			int currentfragpos = subfragstart;

			int newfragstart = -1;
			int newfragstop = -1;
			int newrefstart = -1;
			int newrefstop = -1;

			StringBuilder storedPrecigar = new StringBuilder();
			StringBuilder recentPrecigar = new StringBuilder();
			boolean start = false;
			boolean stop = false;
			for (char c : precigar.toCharArray()) {
				recentPrecigar.append(c);
				switch (c) {

					case 'M':
						if (fragstart * mappedstrand <= currentfragpos * mappedstrand) // we need direction
						{
							start = true;
							if (newrefstart == -1) {
								newrefstart = currentrefpos;
								newfragstart = currentfragpos;
							}
						}
						if ((fragstop + mappedstrand) * mappedstrand <= currentfragpos * mappedstrand) {
							stop = true;
							if (fragstop + mappedstrand == currentfragpos) {
								// if fragstop + mappedstrand < currentfragpos, some extra unmatches occur
								storedPrecigar.append(c);
								newrefstop = currentrefpos - 1;
								newfragstop = currentfragpos - mappedstrand;
							} else
								// Error occurs if you are using exact Match Position to draw sub result
								;

							break;
						} else {
							// continue to update stop position
							newrefstop = currentrefpos - 1; // correct?
							newfragstop = currentfragpos - mappedstrand;
						}
						currentfragpos += mappedstrand;
						currentrefpos += 1;
						break;
					case 'I':
						currentfragpos += mappedstrand;
						break;
					case 'D':
						currentrefpos += 1;
						break;
					default:
						;
				}
				if (stop)
					break;
				if (start)
					storedPrecigar.append(c);
			}
			Cigar newcigar;
			newcigar = new Cigar(storedPrecigar.toString());
			return (new OptMapResultNode(parentFrag, ref.getGenomicPos(newrefstart, newrefstop, true), mappedstrand, newrefstart, newrefstop, newfragstart, newfragstop, newcigar, -1, -1));
		} else {
			throw new IllegalArgumentException("Sub-result argument fails boundary checking. " + parentFrag == null ? "" : ("Result ID: " + parentFrag.name));

		}

	}

	public int getSubRefSigStart() {
		if (subrefstart < 0)
			return -1;
		return subrefstart - 1;
	}

	public int getSubRefSigStop() {
		if (subrefstop < 0)
			return -1;
		return subrefstop;
	}

	public int getSubMoleSigStart() {
		if (subfragstart < 0)
			return -1;
		if (mappedstrand == 1)
			return subfragstart - 1;
		else
			return subfragstart;

	}

	public int getSubMoleSigStop() {
		if (subfragstop < 0)
			return -1;
		if (mappedstrand == 1)
			return subfragstop;
		else
			return subfragstop - 1;
	}

	public void updateMSP() {
		if (cigar == null || subfragstart < 0 || subfragstop < 0 || subrefstart < 0 || subrefstop < 0)
			return;
		List<MatchingSignalPair> msps = new ArrayList<MatchingSignalPair>();

		int currRefSig = subrefstart - 1;
		int currQuerySig = subfragstart;

		if (mappedstrand == 1)
			currQuerySig -= 1;
		for (char c : cigar.getPrecigar().toCharArray()) {
			switch (c) {
				case 'M':
					msps.add(new MatchingSignalPair(currRefSig, currQuerySig));
					currRefSig++;
					currQuerySig += mappedstrand;
					break;
				case 'I':
					currQuerySig += mappedstrand;
					break;
				case 'D':
					currRefSig++;
					break;
				default:
					// Other characters are not supported
					break;
			}
		}
		List<MatchingSignalPair> mspsQ = new ArrayList<MatchingSignalPair>(msps);
		Collections.sort(msps, MatchingSignalPair.qposComparator);
		List<MatchingSignalPair> mspsR = new ArrayList<MatchingSignalPair>(msps);
		Collections.sort(msps, MatchingSignalPair.rposComparator);

		this.mspsQ = mspsQ;
		this.mspsR = mspsR;
	}

	public Set<MatchingSignalPair> getMSPOverlapped(OptMapResultNode result) {
		Set<MatchingSignalPair> mspSet = new HashSet<MatchingSignalPair>();

		int i1 = 0;
		int i2 = 0;
		while (i1 < this.mspsQ.size() && i2 < result.mspsQ.size())
			if (this.mspsQ.get(i1).qpos == result.mspsQ.get(i2).qpos) {
				mspSet.add(this.mspsQ.get(i1));
				i1++;
				i2++;
			} else if (this.mspsQ.get(i1).qpos > result.mspsQ.get(i2).qpos)
				i2++;
			else
				i1++;

		int j1 = 0;
		int j2 = 0;
		while (j1 < this.mspsR.size() && j2 < result.mspsR.size())
			if (this.mspsR.get(j1).rpos == result.mspsR.get(j2).rpos) {
				mspSet.add(this.mspsR.get(j1));
				j1++;
				j2++;
			} else if (this.mspsR.get(j1).rpos > result.mspsR.get(j2).rpos)
				j2++;
			else
				j1++;
		return mspSet;
	}

	public OptMapResultNode getSubResultbyRefPos(DataNode ref, int refstart, int refstop) {
		// input is signal start and stop
		// pter points at every signal
		int[] pter = new int[parentFrag.getTotalSegment() - 1];

		String precigar = this.cigar.getPrecigar();
		int currentrefpos = subrefstart;
		int currentfragpos = subfragstart;

		if (mappedstrand == 1)
			for (int i = 0; i < subfragstart - 1; i++)
				pter[i] = -1;
		else
			for (int i = subfragstart; i < parentFrag.getTotalSegment() - 1; i++)
				pter[i] = -1;
		for (char c : precigar.toCharArray()) {
			switch (c) {
				case 'M':
					if (mappedstrand == 1)
						pter[currentfragpos - 1] = currentrefpos - 1;
					else
						pter[currentfragpos] = currentrefpos - 1;
					currentfragpos += mappedstrand;
					currentrefpos += 1;
					break;
				case 'I':
					if (mappedstrand == 1)
						pter[currentfragpos - 1] = -1;
					else
						pter[currentfragpos] = -1;
					currentfragpos += mappedstrand;
					break;
				case 'D':
					currentrefpos += 1;
					break;
				default:
					;
			}
		}
		if (mappedstrand == 1)
			for (int i = subfragstop + 1; i < parentFrag.getTotalSegment() - 1; i++)
				pter[i] = -1;
		else
			for (int i = 0; i < subfragstop - 1; i++)
				pter[i] = -1;

		int targetfragstart = -1;
		int targetfragstop = -1;
		for (int i = 0; i < parentFrag.getTotalSegment() - 1; i++) {
			if (pter[i] == refstart)
				targetfragstart = i;
			if (pter[i] == refstop)
				targetfragstop = i;
		}
		if (targetfragstart == -1 || targetfragstop == -1)
			return null;
		else {
			if (mappedstrand == 1)
				return getSubResult(ref, targetfragstart + 1, targetfragstop);
			else
				return getSubResult(ref, targetfragstart, targetfragstop + 1);
		}

	}

	public int getTotalSegment() {
		return parentFrag.getTotalSegment();
	}

	public long length() {
		return parentFrag.length();
	}

	public long length(int start, int stop) {
		return parentFrag.length(start, stop);
	}

	@Override
	public String toString() {
		// Name, Aligned Region, Cigar
		if (this.isUsed())
			return String.format("%s\t%s\t%s", parentFrag.toString(), mappedRegion.toString(), cigar);
		else
			return parentFrag.toString();
	}

	// Modify
	public void reverse() {
		parentFrag = parentFrag.getReverse();
		this.mappedstrand *= -1;
		this.subfragstart = parentFrag.getTotalSegment() - subfragstart - 1;
		this.subfragstop = parentFrag.getTotalSegment() - subfragstop - 1;
	}

	public void scale(double ratio) {
		parentFrag = new DataNode(parentFrag);
		parentFrag.scale(ratio);
		// super.scale(ratio);
	}

	private void trimResult(int step) {
		String precigar = cigar.getPrecigar();
		if (step == 0) {
			// updateMSP();
			return;
		} else if (step < 0) // right to left trim, according to ref
			precigar = StringUtils.reverse(precigar);
		if (precigar.length() == 1) // will trim to nothing, cannot trim anymore
			return;
		String removed = precigar.substring(0, precigar.indexOf('M', 1));
		;
		int m = 0;
		int i = 0;
		int d = 0;
		for (char c : removed.toCharArray())
			switch (c) {
				case 'M':
					m++;
					break;
				case 'I':
					i++;
					break;
				case 'D':
					d++;
					break;
				default:
					break;
			}
		if (step > 0) {
			subrefstart += m + d;
			if (mappedstrand == -1)
				subfragstart -= m + i; // + or -
			else
				subfragstart += m + i;
		} else {
			subrefstop -= m + d; // + or -
			if (mappedstrand == -1)
				subfragstop += m + i;
			else
				subfragstop -= m + i; // + or -
		}
		String newprecigar = precigar.substring(precigar.indexOf('M', 1), precigar.length());
		if (step < 0)
			newprecigar = StringUtils.reverse(newprecigar);
		this.cigar = new Cigar(newprecigar);
		trimResult(step - Math.abs(step) / step);
	}

	public void trimResult(int step, LinkedHashMap<String, DataNode> optrefmap) {
		this.trimResult(step, optrefmap.get(mappedRegion.ref));
	}

	public void trimResult(int step, DataNode ref) {
		this.trimResult(step);
		this.updateMappedRegion(ref);
	}

	public void updateScore(LinkedHashMap<String, DataNode> optrefmap, int match, int fpp, int fnp) {
		this.mappedscore = (cigar.calcScore(match, fpp, fnp)) * (1 - Math.abs(1 - getMapScale()));
	}

	public void updateMappedRegion(DataNode ref) {
		long mappedstart;
		long mappedstop;
		if (subrefstart >= 1)
			mappedstart = ref.refp[subrefstart - 1];
		else
			mappedstart = 0;
		if (subrefstop < ref.refp.length)
			mappedstop = ref.refp[subrefstop];
		else
			mappedstop = ref.size;
		mappedRegion = new GenomicPosNode(ref.name, mappedstart, mappedstop);
	}

	// Test
	public boolean correctlyMapped() {
		// double defaultMaxMappedFactor = 0.50;
		// return correctlyMapped(defaultMaxMappedFactor);
		if (parentFrag == null)
			return false;
		if (!parentFrag.hasSimulationInfo())
			return false;
		return correctlyMapped(0, parentFrag.simuInfo);
	}

	public boolean correctlyMapped(SimulationInfo simuInfo) {
		return correctlyMapped(0, simuInfo);
	}

	public boolean correctlyMapped(long tolerance, SimulationInfo simuInfo) {
		return mappedRegion.isClose(simuInfo.simuRegion, tolerance) && mappedstrand == simuInfo.simuStrand;
	}

	public boolean isRefClose(OptMapResultNode map, long gapAllowed) {
		return this.mappedRegion.isClose(map.mappedRegion, gapAllowed);
	}

	public boolean isFragClose(OptMapResultNode map, long gapAllowed) {
		SimpleLongLocation s1;
		if (mappedstrand == 1)
			s1 = new SimpleLongLocation(length(0, subfragstart - 1), length(0, subfragstop));
		else if (mappedstrand == -1)
			s1 = new SimpleLongLocation(length(0, subfragstop - 1), length(0, subfragstart));
		else
			s1 = null;
		SimpleLongLocation s2;
		if (map.mappedstrand == 1)
			s2 = new SimpleLongLocation(length(0, map.subfragstart - 1), length(0, map.subfragstop));
		else if (map.mappedstrand == -1)
			s2 = new SimpleLongLocation(length(0, map.subfragstop - 1), length(0, map.subfragstart));
		else
			s2 = null;

		return (s1.overlap(s2, gapAllowed));

	}

	public boolean isClose(OptMapResultNode map, long gapAllowed) {
		return (this.isFragClose(map, gapAllowed) && this.isRefClose(map, gapAllowed));
	}

	public boolean isClose(GenomicPosNode region, long gapAllowed) {
		return this.mappedRegion.isClose(region, gapAllowed);
	}

	public boolean isSimilar(OptMapResultNode map, long gapAllowed) {
		return isClose(map, gapAllowed) && (this.mappedstrand == map.mappedstrand);
	}

	public boolean isUsed() {
		if (mappedRegion == null)
			return false;
		else if (mappedRegion.ref.equalsIgnoreCase("Discarded"))
			return false;
		else if (mappedRegion.ref.equalsIgnoreCase("Unmapped"))
			return false;
		return true;
	}

	public Boolean isSubFragInfoValid() {
		if (parentFrag.refp == null) {
			System.err.println("No fragment information for checking.");
			return null;
		}
		int max = parentFrag.getTotalSegment() - 1;
		int min = 0;
		return (subfragstart <= max && subfragstart >= min && subfragstop <= max && subfragstop >= min);
	}

	public boolean isSubRefInfoValid(DataNode ref) {
		int max = ref.refp.length;
		int min = 0;
		return (subrefstart <= max && subrefstart >= min && subrefstop <= max && subrefstop >= min) && (subrefstart <= subrefstop);

	}

	public boolean isSubRefInfoValid(LinkedHashMap<String, DataNode> optrefmap) {
		if (optrefmap.containsKey(mappedRegion.ref))
			return isSubRefInfoValid(optrefmap.get(mappedRegion.ref));
		else
			return false;
	}

	public boolean mapSignal(String refname, int targetSig) {
		if (!mappedRegion.ref.equalsIgnoreCase(refname))
			return false;
		int currRefSig = subrefstart - 1;
		for (char c : cigar.getPrecigar().toCharArray()) {
			switch (c) {
				case 'M':
					if (targetSig == currRefSig)
						return true;
					currRefSig++;
					break;
				case 'I':
					break;
				case 'D':
					currRefSig++;
					break;
				default:
					break;
			}
			if (currRefSig > targetSig)
				return false;
		}
		return false;
	}

	public boolean overlap(OptMapResultNode result) {
		return (overlapQuery(result) || overlapRef(result));
	}

	public boolean overlapQuery(OptMapResultNode result) {
		SimpleLocation thisQueryLoc = new SimpleLocation(this.subfragstart, this.subfragstop);
		SimpleLocation resultQueryloc = new SimpleLocation(result.subfragstart, result.subfragstop);

		return (thisQueryLoc.overlapsize(resultQueryloc) > -1); // even subfragstart not overlap is not enough; the signal is still shared.
	}

	public boolean overlapRef(OptMapResultNode map) {
		SimpleLocation thisrefloc = new SimpleLocation(this.subrefstart, this.subrefstop);
		SimpleLocation maprefloc = new SimpleLocation(map.subrefstart, map.subrefstop);
		return ((this.mappedRegion.ref.equals(map.mappedRegion.ref) && thisrefloc.overlapsize(maprefloc) > -1)); // even subrefstart not overlap is not enough; the signal is still shared.
	}

	// Copy
	public List<OptMapResultNode> getRealMap() {
		List<OptMapResultNode> mapList = new ArrayList<OptMapResultNode>();
		mapList.add(this);
		return mapList;
	}

	public OptMapResultNode newReverseFragment() {
		OptMapResultNode map = new OptMapResultNode(this);
		map.reverse();
		return map;
	}

	// Static
	public static OptMapResultNode newBlankMapNode(DataNode f) {
		return new OptMapResultNode(f, null, 0, -1, -1, -1, -1, null, Double.NEGATIVE_INFINITY, -1);
	}

	public static List<GenomicPosNode> getMappedRegion(List<OptMapResultNode> resultList) {
		if (resultList == null)
			throw new NullPointerException("resultList");
		List<GenomicPosNode> posList = new ArrayList<>();
		for (OptMapResultNode result : resultList)
			posList.add(result.mappedRegion);
		return posList;
	}

	public static List<GenomicPosNode> getMoleMappedRegion(List<OptMapResultNode> resultList) {
		if (resultList == null)
			throw new NullPointerException("resultList");
		List<GenomicPosNode> posList = new ArrayList<>();
		for (OptMapResultNode result : resultList)
			posList.add(result.getMoleMappedRegion());
		return posList;
	}

	public static List<GenomicPosNode> getPotentiallyMappedRegion(LinkedHashMap<String, DataNode> optrefmap, List<OptMapResultNode> resultList) {
		List<GenomicPosNode> targetRegionList = new ArrayList<GenomicPosNode>();
		if (resultList.isEmpty())
			return targetRegionList;
		if (!resultList.get(0).isUsed())
			return targetRegionList;

		for (OptMapResultNode result : resultList) {
			String ref = result.mappedRegion.ref;
			long start = result.mappedRegion.start - result.parentFrag.size;
			long stop = result.mappedRegion.stop + result.parentFrag.size;
			if (start < 1)
				start = 1;
			if (stop > optrefmap.get(ref).size)
				stop = optrefmap.get(ref).size;
			targetRegionList.add(new GenomicPosNode(ref, start, stop));
		}
		return GenomicPosNode.merge(targetRegionList);
	}

	public static LinkedHashMap<String, List<GenomicPosNode>> getPotentiallyMappedRegion(LinkedHashMap<String, DataNode> optrefmap, LinkedHashMap<String, List<OptMapResultNode>> resultListMap) {

		LinkedHashMap<String, List<GenomicPosNode>> targetRegionMap = new LinkedHashMap<>();
		for (List<OptMapResultNode> resultList : resultListMap.values()) {
			List<GenomicPosNode> targetRegionList = getPotentiallyMappedRegion(optrefmap, resultList);
			targetRegionMap.put(resultList.get(0).parentFrag.name, targetRegionList);
		}
		return targetRegionMap;
	}

	public static int[] getMapSignal(List<OptMapResultNode> resultList) {
		Set<Integer> set = new HashSet<Integer>();
		for (OptMapResultNode result : resultList) {
			for (int i : result.getMapSignal())
				set.add(i);
		}
		List<Integer> list = new ArrayList<Integer>(set);
		Collections.sort(list);
		return ArrayUtils.toPrimitive(list.toArray(new Integer[list.size()]));
	}

	public static List<SimpleLocation> getUnmapPortion(List<OptMapResultNode> fragmentmaplist) {
		if (fragmentmaplist == null || fragmentmaplist.isEmpty())
			return null;
		DataNode f = fragmentmaplist.get(0).parentFrag;
		// boolean[] usedSegment = new boolean[f.getTotalSegment()];
		List<SimpleLocation> subfraglist = new ArrayList<SimpleLocation>();
		for (OptMapResultNode map : fragmentmaplist)
			subfraglist.add(new SimpleLocation(map.subfragstart, map.subfragstop));
		SimpleLocation.mergeLocationList(subfraglist);
		List<SimpleLocation> remainfraglist = new ArrayList<SimpleLocation>();
		remainfraglist.add(new SimpleLocation(0, f.getTotalSegment() - 1));
		SimpleLocation.removeSimpleLocationList(remainfraglist, subfraglist);
		return remainfraglist;
	}

	public static List<OptMapResultNode> reconstruct(LinkedHashMap<String, DataNode> optrefmap, OptMapResultNode map) {
		// Only reconstruct the subrefstart,.... and subfragratio, and cigar string
		int direction = map.mappedstrand;
		List<OptMapResultNode> fragmentmaplist = new ArrayList<OptMapResultNode>();

		// String precigar = Cigar.convertpreCIGAR(map.cigar);
		String precigar = map.cigar.getPrecigar();
		// if (direction == -1)
		// precigar = StringUtils.reverse(precigar);
		String[] cigarlist = precigar.split("S");
		for (int i = 1; i <= cigarlist.length; i += 2) {
			int subrefstart = -1;
			int subrefstop = -1;
			int subfragstart = -1;
			int subfragstop = -1;
			// double scale = 0;
			if (i - 2 < 0) {
				subrefstart = map.subrefstart;
				subfragstart = map.subfragstart;
			} else {
				StringBuilder concat = new StringBuilder();
				for (int j = 0; j < i - 1; j++)
					concat.append(cigarlist[j]);
				String previousprecigar = concat.toString();
				int match = Cigar.getCertainNumberFromPrecigar(previousprecigar, 'M');
				int insert = Cigar.getCertainNumberFromPrecigar(previousprecigar, 'I');
				int delete = Cigar.getCertainNumberFromPrecigar(previousprecigar, 'D');
				subfragstart = map.subfragstart + (match + insert) * direction;
				subrefstart = map.subrefstart + match + delete;
			}
			String recentprecigar = cigarlist[i - 1];
			if (i == cigarlist.length) {
				subrefstop = map.subrefstop;
				subfragstop = map.subfragstop;
			} else {
				int match = Cigar.getCertainNumberFromPrecigar(recentprecigar, 'M');
				int insert = Cigar.getCertainNumberFromPrecigar(recentprecigar, 'I');
				int delete = Cigar.getCertainNumberFromPrecigar(recentprecigar, 'D');
				subfragstop = subfragstart + (match + insert - 1) * direction;
				subrefstop = subrefstart + match + delete - 1;
			}
			DataNode ref = optrefmap.get(map.mappedRegion.ref);
			long estimatestartpos = ref.refp[subrefstart - 1];
			long estimatestoppos = ref.refp[subrefstop];

			Cigar newcigar = new Cigar(recentprecigar);
			fragmentmaplist.add(new OptMapResultNode(map.parentFrag, new GenomicPosNode(map.mappedRegion.ref, estimatestartpos, estimatestoppos), map.mappedstrand, subrefstart, subrefstop,
					subfragstart, subfragstop, newcigar, map.mappedscore, -1));
		}
		return fragmentmaplist;
	}

	public static OptMapResultNode reverseRefAndFrag(OptMapResultNode result, DataNode originalRef) {
		DataNode newRef = new DataNode(result.parentFrag);
		OptMapResultNode newResult;
		int newsubfragstart;
		int newsubfragstop;
		int newsubrefstart;
		int newsubrefstop;
		Cigar newcigar;
		if (result.mappedstrand == 1) {
			newsubfragstart = result.subrefstart;
			newsubfragstop = result.subrefstop;
			newsubrefstart = result.subfragstart;
			newsubrefstop = result.subfragstop;
			newcigar = new Cigar(result.cigar);
		} else {
			newsubfragstart = result.subrefstop;
			newsubfragstop = result.subrefstart;
			newsubrefstart = result.subfragstop;
			newsubrefstop = result.subfragstart;
			newcigar = result.cigar.getReverseCigar();
		}
		newResult = new OptMapResultNode(new DataNode(originalRef), null, result.mappedstrand, newsubrefstart, newsubrefstop, newsubfragstart, newsubfragstop, newcigar, result.mappedscore,
				result.confidence);
		newResult.updateMappedRegion(newRef);
		return newResult;
	}

	public static OptMapResultNode reverseRefAndFrag(OptMapResultNode result, LinkedHashMap<String, DataNode> originalOptRefMap) {
		return OptMapResultNode.reverseRefAndFrag(result, originalOptRefMap.get(result.mappedRegion.ref));
	}

	public static List<OptMapResultNode> reverseRefAndFrag(List<OptMapResultNode> resultList, DataNode originalRef) {
		List<OptMapResultNode> newResultList = new ArrayList<OptMapResultNode>();
		for (OptMapResultNode result : resultList)
			newResultList.add(OptMapResultNode.reverseRefAndFrag(result, originalRef));
		return newResultList;
	}

	public static List<OptMapResultNode> reverseRefAndFrag(List<OptMapResultNode> resultList, LinkedHashMap<String, DataNode> originalOptRefMap) {
		List<OptMapResultNode> newResultList = new ArrayList<OptMapResultNode>();
		for (OptMapResultNode result : resultList)
			newResultList.add(OptMapResultNode.reverseRefAndFrag(result, originalOptRefMap));
		return newResultList;
	}

	public static Comparator<OptMapResultNode> mappedstartcomparator = new Comparator<OptMapResultNode>() {
		@Override
		public int compare(OptMapResultNode f1, OptMapResultNode f2) {
			return f1.mappedRegion.compareTo(f2.mappedRegion);
		}
	};
	public static Comparator<OptMapResultNode> subfragstartcomparator = new Comparator<OptMapResultNode>() {
		@Override
		public int compare(OptMapResultNode f1, OptMapResultNode f2) {
			int i = Integer.valueOf(f1.subfragstart).compareTo(Integer.valueOf(f2.subfragstart));
			if (i == 0)
				i = Integer.valueOf(f1.subfragstop).compareTo(Integer.valueOf(f2.subfragstop));
			return i;
		}
	};
	public static Comparator<OptMapResultNode> subfragstartstopcomparator = new Comparator<OptMapResultNode>() {
		// consider the comparison when reverse and forward strand is together
		@Override
		public int compare(OptMapResultNode f1, OptMapResultNode f2) {
			int i = Integer.valueOf(f1.mappedstrand == 1 ? f1.subfragstart : f1.subfragstop).compareTo(Integer.valueOf(f2.mappedstrand == 1 ? f2.subfragstart : f2.subfragstop));
			if (i == 0)
				i = Integer.valueOf(f1.mappedstrand == 1 ? f1.subfragstop : f1.subfragstart).compareTo(Integer.valueOf(f2.mappedstrand == 1 ? f2.subfragstop : f2.subfragstart));
			return i;
		}
	};
	public static Comparator<OptMapResultNode> mappedscorecomparator = new Comparator<OptMapResultNode>() {
		@Override
		public int compare(OptMapResultNode f1, OptMapResultNode f2) {
			return Double.valueOf(f1.mappedscore).compareTo(Double.valueOf(f2.mappedscore));
		}
	};

}

class MatchingSignalPair {
	public final int rpos;
	public final int qpos;

	public MatchingSignalPair(int rpos, int qpos) {
		this.rpos = rpos;
		this.qpos = qpos;
	}

	public boolean equals(MatchingSignalPair msp) {
		return this.rpos == msp.rpos && this.qpos == msp.qpos;
	}

	public boolean overlap(MatchingSignalPair msp) {
		return this.rpos == msp.rpos || this.qpos == msp.qpos;
	}

	// In a list of MSPs in an alignment, qpos and rpos should be unique
	public static Comparator<MatchingSignalPair> rposComparator = new Comparator<MatchingSignalPair>() {
		@Override
		public int compare(MatchingSignalPair msp1, MatchingSignalPair msp2) {
			return Integer.compare(msp1.rpos, msp2.rpos);
		}
	};
	public static Comparator<MatchingSignalPair> qposComparator = new Comparator<MatchingSignalPair>() {
		@Override
		public int compare(MatchingSignalPair msp1, MatchingSignalPair msp2) {
			return Integer.compare(msp1.qpos, msp2.qpos);
		}
	};
}
