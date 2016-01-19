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


package aldenjava.common;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/** A class with two long integers describing a 1-dimension range
 *  
 * @author Alden
 *
 */
 
public class SimpleLongLocation implements Cloneable {
	public long min;
	public long max;

	public SimpleLongLocation(long x1, long x2) {
		if (x1 > x2) {
			this.max = x1;
			this.min = x2;
		} else {
			this.max = x2;
			this.min = x1;
		}
	}

	public SimpleLongLocation(SimpleLongLocation s) {
		this.min = s.min;
		this.max = s.max;
	}

	public long length() {
		return this.max - this.min + 1;
	}

	public boolean overlap(SimpleLongLocation s) {
		return overlap(s, 0);
	}

	public boolean overlap(SimpleLongLocation s, long allowedgap) {
		if (overlapsize(s) >= allowedgap * -1 + 1)
			return true;
		else
			return false;
	}

	public long overlapsize(SimpleLongLocation s) {
		long max = s.max > this.max ? this.max : s.max;
		long min = s.min < this.min ? this.min : s.min;
		return max - min + 1;
	}

	public boolean contain(long i) {
		if ((this.min <= i) && (this.max >= i))
			return true;
		else
			return false;
	}

	public boolean contain(SimpleLongLocation loc) {
		return ((this.min <= loc.min) && (this.max >= loc.max));
	}

	public long average() {
		return (min + max) / 2;
	}

	public boolean equalLoc(SimpleLongLocation loc) {
		return (this.min == loc.min) && (this.max == loc.max);
	}

	@Override
	public String 
	toString() {
		return (Long.toString(this.min) + "\t" + Long.toString(this.max));
	}

	public static List<SimpleLongLocation> createSimpleLongLocationList(List<Long> x1, List<Long> x2) {
		List<SimpleLongLocation> loc = new ArrayList<SimpleLongLocation>();
		if (x1.size() != x2.size())
			System.err.println("Size mismatch!");
		else {
			for (int i = 0; i < x1.size(); i++)
				loc.add(new SimpleLongLocation(x1.get(i), x2.get(i)));
		}
		return loc;
	}

	public static SimpleLongLocation[] createSimpleLongLocationArray(int[] x1, int[] x2) {
		SimpleLongLocation[] loc = new SimpleLongLocation[x1.length];
		if (x1.length != x2.length)
			System.err.println("Size mismatch!");
		else {
			for (int i = 0; i < x1.length; i++)
				loc[i] = new SimpleLongLocation(x1[i], x2[i]);
		}
		return loc;
	}

	public static List<SimpleLongLocation> mergeLocationList(List<SimpleLongLocation> locationlist) {
		List<SimpleLongLocation> newloclist = new ArrayList<SimpleLongLocation>();
		Collections.sort(locationlist, SimpleLongLocation.minc);
		SimpleLongLocation lastloc = null;
		for (SimpleLongLocation loc : locationlist) {
			if (lastloc == null)
				lastloc = new SimpleLongLocation(loc);
			if (lastloc.overlapsize(loc) >= 0) // no space between can be merged too
			{
				if (loc.max > lastloc.max)
					lastloc = new SimpleLongLocation(lastloc.min, loc.max);
			} else {
				newloclist.add(lastloc);
				lastloc = loc;
			}
		}
		if (lastloc != null)
			newloclist.add(lastloc);
		return newloclist;
	}

	public static List<SimpleLongLocation> getConsensusLocationList(List<SimpleLongLocation> loclist1, List<SimpleLongLocation> loclist2) {
		List<SimpleLongLocation> newloclist = new ArrayList<SimpleLongLocation>();
		int locsource = 0;
		int p1 = 0;
		int p2 = 0;
		SimpleLongLocation loc = null;
		while ((p1 < loclist1.size() || locsource == 1) && (p2 < loclist2.size() || locsource == 2)) {
			switch (locsource) {
				case 0: {
					if (loclist1.get(p1).min <= loclist2.get(p2).min) {
						locsource = 1;
						loc = loclist1.get(p1);
						p1++;
					} else {
						locsource = 2;
						loc = loclist2.get(p2);
						p2++;
					}
					break;
				}
				case 1: {
					if (loclist2.get(p2).overlap(loc)) {
						if (loclist2.get(p2).max > loc.max) {
							newloclist.add(new SimpleLongLocation(loclist2.get(p2).min, loc.max));
							locsource = 2;
							loc = loclist2.get(p2);
						} else
							newloclist.add(new SimpleLongLocation(loclist2.get(p2).min, loclist2.get(p2).max));
						p2++;
					} else {
						loc = null;
						locsource = 0;
					}
					break;
				}
				case 2: {
					if (loclist1.get(p1).overlap(loc)) {
						if (loclist1.get(p1).max > loc.max) {
							newloclist.add(new SimpleLongLocation(loclist1.get(p1).min, loc.max));
							locsource = 1;
							loc = loclist1.get(p1);
						} else
							newloclist.add(new SimpleLongLocation(loclist1.get(p1).min, loclist1.get(p1).max));
						p1++;
					} else {
						loc = null;
						locsource = 0;
					}
					break;
				}
				default:
					return null;
			}

		}
		return newloclist;
	}

	public static List<SimpleLongLocation> removeSimpleLongLocationList(List<SimpleLongLocation> targetlist, List<SimpleLongLocation> loclist) // merged + sorted list
	{
		List<SimpleLongLocation> newloclist = new ArrayList<SimpleLongLocation>();
		for (SimpleLongLocation targetloc : targetlist) {
			SimpleLongLocation newtargetloc = new SimpleLongLocation(targetloc);
			for (SimpleLongLocation loc : loclist)
				if (newtargetloc.overlap(loc)) {
					if (newtargetloc.min < loc.min)
						newloclist.add(new SimpleLongLocation(newtargetloc.min, loc.min - 1));
					if (newtargetloc.max > loc.max)
						newtargetloc.min = loc.max + 1;
					else {
						newtargetloc = null;
						break;
					}
				}
			if (newtargetloc != null)
				newloclist.add(newtargetloc);
		}
		return newloclist;
	}

	public static long getMin(List<SimpleLongLocation> locList) {
		long min = Long.MAX_VALUE;
		for (SimpleLongLocation loc : locList)
			if (min > loc.min)
				min = loc.min;
		return min;
	}

	public static long getMax(List<SimpleLongLocation> locList) {
		long max = Long.MIN_VALUE;
		for (SimpleLongLocation loc : locList)
			if (max < loc.max)
				max = loc.max;
		return max;
	}

	public static Comparator<SimpleLongLocation> maxc = new Comparator<SimpleLongLocation>() {
		@Override
		public int compare(SimpleLongLocation u1, SimpleLongLocation u2) {
			return Long.valueOf(u1.max).compareTo(Long.valueOf(u2.max));
		}
	};

	public static Comparator<SimpleLongLocation> minc = new Comparator<SimpleLongLocation>() {
		@Override
		public int compare(SimpleLongLocation u1, SimpleLongLocation u2) {
			return Long.valueOf(u1.min).compareTo(Long.valueOf(u2.min));
		}
	};

}