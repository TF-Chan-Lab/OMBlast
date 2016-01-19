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

/** A class with two integers describing a 1-dimension range
 *  
 * @author Alden
 *
 */
public class SimpleLocation {
	public int min;
	public int max;

	public SimpleLocation(int x1, int x2) {
		if (x1 > x2) {
			this.max = x1;
			this.min = x2;
		} else {
			this.max = x2;
			this.min = x1;
		}
	}

	public SimpleLocation(SimpleLocation s) {
		this.min = s.min;
		this.max = s.max;
	}

	public boolean overlap(SimpleLocation s) {
		return overlap(s, 0);
	}

	public boolean overlap(SimpleLocation s, long allowedgap) {
		if (overlapsize(s) >= allowedgap * -1 + 1)
			return true;
		else
			return false;
	}

	public int overlapsize(SimpleLocation s) {
		int max = s.max > this.max ? this.max : s.max;
		int min = s.min < this.min ? this.min : s.min;
		return max - min + 1;
	}

	public boolean contain(int i) {
		if ((this.min <= i) && (this.max >= i))
			return true;
		else
			return false;
	}

	public boolean contain(SimpleLocation loc) {
		return ((this.min <= loc.min) && (this.max >= loc.max));
	}

	public int average() {
		return (min + max) / 2;
	}

	public boolean equalLoc(SimpleLocation loc) {
		return (this.min == loc.min) && (this.max == loc.max);
	}

	@Override
	public String toString() {
		return (Integer.toString(this.min) + "\t" + Integer.toString(this.max));
	}

	public static List<SimpleLocation> createSimpleLocationList(List<Integer> x1, List<Integer> x2) {
		List<SimpleLocation> loc = new ArrayList<SimpleLocation>();
		if (x1.size() != x2.size())
			System.err.println("Size mismatch!");
		else {
			for (int i = 0; i < x1.size(); i++)
				loc.add(new SimpleLocation(x1.get(i), x2.get(i)));
		}
		return loc;
	}

	public static SimpleLocation[] createSimpleLocationArray(int[] x1, int[] x2) {
		SimpleLocation[] loc = new SimpleLocation[x1.length];
		if (x1.length != x2.length)
			System.err.println("Size mismatch!");
		else {
			for (int i = 0; i < x1.length; i++)
				loc[i] = new SimpleLocation(x1[i], x2[i]);
		}
		return loc;
	}

	public static List<SimpleLocation> mergeLocationList(List<SimpleLocation> locationlist) {
		List<SimpleLocation> newloclist = new ArrayList<SimpleLocation>();
		Collections.sort(locationlist, SimpleLocation.minc);
		SimpleLocation lastloc = null;
		for (SimpleLocation loc : locationlist) {
			if (lastloc == null)
				lastloc = new SimpleLocation(loc);
			if (lastloc.overlapsize(loc) >= 0) // no space between can be merged too
			{
				if (loc.max > lastloc.max)
					lastloc = new SimpleLocation(lastloc.min, loc.max);
			} else {
				newloclist.add(lastloc);
				lastloc = loc;
			}
		}
		if (lastloc != null)
			newloclist.add(lastloc);
		return newloclist;
	}

	public static List<SimpleLocation> getConsensusLocationList(List<SimpleLocation> loclist1, List<SimpleLocation> loclist2) {
		List<SimpleLocation> newloclist = new ArrayList<SimpleLocation>();
		int locsource = 0;
		int p1 = 0;
		int p2 = 0;
		SimpleLocation loc = null;
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
							newloclist.add(new SimpleLocation(loclist2.get(p2).min, loc.max));
							locsource = 2;
							loc = loclist2.get(p2);
						} else
							newloclist.add(new SimpleLocation(loclist2.get(p2).min, loclist2.get(p2).max));
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
							newloclist.add(new SimpleLocation(loclist1.get(p1).min, loc.max));
							locsource = 1;
							loc = loclist1.get(p1);
						} else
							newloclist.add(new SimpleLocation(loclist1.get(p1).min, loclist1.get(p1).max));
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

	public static List<SimpleLocation> removeSimpleLocationList(List<SimpleLocation> targetlist, List<SimpleLocation> loclist) // merged + sorted list
	{
		List<SimpleLocation> newloclist = new ArrayList<SimpleLocation>();
		for (SimpleLocation targetloc : targetlist) {
			SimpleLocation newtargetloc = new SimpleLocation(targetloc);
			for (SimpleLocation loc : loclist)
				if (newtargetloc.overlap(loc)) {
					if (newtargetloc.min < loc.min)
						newloclist.add(new SimpleLocation(newtargetloc.min, loc.min - 1));
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

	public static Comparator<SimpleLocation> maxc = new Comparator<SimpleLocation>() {
		@Override
		public int compare(SimpleLocation u1, SimpleLocation u2) {
			return Integer.valueOf(u1.max).compareTo(Integer.valueOf(u2.max));
		}
	};

	public static Comparator<SimpleLocation> minc = new Comparator<SimpleLocation>() {
		@Override
		public int compare(SimpleLocation u1, SimpleLocation u2) {
			return Integer.valueOf(u1.min).compareTo(Integer.valueOf(u2.min));
		}
	};

}