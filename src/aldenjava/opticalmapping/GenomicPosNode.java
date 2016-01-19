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


package aldenjava.opticalmapping;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import joptsimple.OptionSet;
import aldenjava.common.SimpleLongLocation;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;

/**
 * GenomicPosNode objects are immutable
 * 
 * @author Alden
 *
 */
public final class GenomicPosNode implements Comparable<GenomicPosNode> {

	public final String ref;
	public final long start;
	public final long stop;
	private final static String separator = "[:\\-\\s]+"; // :, -, space

	public GenomicPosNode(String s) throws IllegalArgumentException {
		String[] l = s.split(GenomicPosNode.separator);
		if (s.isEmpty() || (l.length != 2 && l.length != 3))
			throw new IllegalArgumentException("Wrong genomic position format.");
		s = s.trim();

		String ref = l[0];
		long start;
		long stop;
		try {
			start = Long.parseLong(l[1]);
			if (l.length == 2)
				stop = start;
			else
				stop = Long.parseLong(l[2]);
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException("Wrong genomic position format.");
		}
		this.ref = ref;
		this.start = start;
		this.stop = stop;

	}

	public GenomicPosNode(String ref, long start) {
		this.ref = ref;
		this.start = start;
		this.stop = start;
	}

	public GenomicPosNode(String ref, long start, long stop) {
		this.ref = ref;
		this.start = start;
		this.stop = stop;
	}

	public GenomicPosNode(GenomicPosNode region) {
		this.ref = region.ref;
		this.start = region.start;
		this.stop = region.stop;
	}

	public GenomicPosNode(OptionSet options) {
		this((String) options.valueOf("region"));
	}

	public boolean isClose(GenomicPosNode region, long gapAllowed) {
		if (!this.ref.equals(region.ref))
			return false;
		SimpleLongLocation s1 = this.getLoc();
		SimpleLongLocation s2 = region.getLoc();
		return (s1.overlap(s2, gapAllowed));
	}

	public boolean isBPClose(GenomicPosNode region, long gapAllowed) {
		GenomicPosNode myregion1 = new GenomicPosNode(ref, start, start);
		GenomicPosNode myregion2 = new GenomicPosNode(ref, stop, stop);
		GenomicPosNode region1 = new GenomicPosNode(region.ref, region.start, region.start);
		GenomicPosNode region2 = new GenomicPosNode(region.ref, region.stop, region.stop);
		return myregion1.isClose(region1, gapAllowed) && myregion2.isClose(region2, gapAllowed);
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null)
			return false;
		if (!(obj instanceof GenomicPosNode))
			return false;
		GenomicPosNode region = (GenomicPosNode) obj;
		return region.ref.equals(ref) && region.start == start && region.stop == stop;
	}

	public long length() {
		return (stop - start + 1);
	}

	public SimpleLongLocation getLoc() {
		return new SimpleLongLocation(start, stop);
	}

	public long overlapSize(GenomicPosNode e) {
		if (this.ref.equalsIgnoreCase(e.ref))
			return this.getLoc().overlapsize(e.getLoc());
		else
			return 0;
	}

	@Override
	public String toString() {
		if (start != stop)
			return String.format("%s:%d-%d", ref, start, stop);
		else
			return String.format("%s:%d", ref, start);
	}

	public String toFileString() {
		if (start != stop)
			return String.format("%s-%d-%d", ref, start, stop);
		else
			return String.format("%s-%d", ref, start);
	}

	@Override
	public int compareTo(GenomicPosNode g2) {
		if (this.ref.compareToIgnoreCase(g2.ref) != 0)
			return this.ref.compareToIgnoreCase(g2.ref);
		else if (this.start != g2.start)
			return Long.compare(this.start, g2.start);
		else
			return Long.compare(this.stop, g2.stop);
	}

	public static String format() {
		return "chrX:NNNNN-NNNNN";
	}

	public static List<SimpleLongLocation> getLoc(List<GenomicPosNode> regionList) {
		if (regionList == null)
			throw new NullPointerException("regionList");
		List<SimpleLongLocation> locList = new ArrayList<>();
		for (GenomicPosNode region : regionList) {
			locList.add(region.getLoc());
		}
		return locList;
	}

	public static List<GenomicPosNode> merge(List<GenomicPosNode> regionList) {
		Collections.sort(regionList);
		String currRef = null;
		long currStart = -1;
		long currStop = -1;
		List<GenomicPosNode> mergedRegionList = new ArrayList<GenomicPosNode>();
		for (GenomicPosNode region : regionList) {
			if (region == null)
				continue;

			if (currRef != null)
				if (region.ref.equalsIgnoreCase(currRef) && region.start < currStop) {
					if (region.stop > currStop)
						currStop = region.stop;
				} else {
					mergedRegionList.add(new GenomicPosNode(currRef, currStart, currStop));
					currRef = region.ref;
					currStart = region.start;
					currStop = region.stop;
				}
			else {
				currRef = region.ref;
				currStart = region.start;
				currStop = region.stop;
			}
		}
		if (currRef != null)
			mergedRegionList.add(new GenomicPosNode(currRef, currStart, currStop));
		return mergedRegionList;
	}

	public static void assignOptions(ExtendOptionParser parser, int level) {
		parser.addHeader("Region", level);
		parser.accepts("region", "The region of interest.").withOptionalArg().ofType(String.class);
	}

	public boolean contains(GenomicPosNode region) {
		return (region.ref.equalsIgnoreCase(ref) && (start <= region.start && stop >= region.stop));
	}
}
