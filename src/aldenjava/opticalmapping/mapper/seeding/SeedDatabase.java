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


package aldenjava.opticalmapping.mapper.seeding;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import joptsimple.OptionSet;
import aldenjava.opticalmapping.GenomicPosNode;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;
import aldenjava.opticalmapping.miscellaneous.SelectableMode;

/**
 * The SeedDatabase stores reference kmers. It offers a list of reference kmers as seeds upon a query kmer request
 * 
 * @author Alden
 *
 */
public class SeedDatabase implements SelectableMode {

	private int seedingmode;
	private final LinkedHashMap<String, DataNode> optrefmap;
	private List<List<Kmer>> databaseSeedList;
	private List<GenomicPosNode> restrictedRegions;
	private FastConversionTable table = FastConversionTable.standardTable();
	private LinkedHashMap<String, List<Kmer>> fastDatabaseSeedMap;

	private int k;
	private int maxnosignalregion;
	private boolean referenceChanged;

	public SeedDatabase(LinkedHashMap<String, DataNode> optrefmap) {
		this.optrefmap = optrefmap;
		this.referenceChanged = true;
	}

	public SeedDatabase(LinkedHashMap<String, DataNode> optrefmap, List<List<Kmer>> databaseSeedList) {
		this.optrefmap = optrefmap;
		this.referenceChanged = false;
		this.databaseSeedList = databaseSeedList;
	}

	public SeedDatabase(List<List<Kmer>> databaseSeedList) {
		this.optrefmap = null;
		this.referenceChanged = false;
		this.databaseSeedList = databaseSeedList;
	}

	public SeedDatabase(List<Kmer> kmerList, int kmerlen) {
		this.optrefmap = null;
		this.databaseSeedList = convertKmerList(kmerList, kmerlen);
	}

	@Override
	public void setMode(OptionSet options) {
		setMode((int) options.valueOf("seedingmode"));
	}

	@Override
	public void setMode(int mode) {
		this.seedingmode = mode;
	}

	@Override
	public int getMode() {
		return seedingmode;
	}

	public void setParameters(OptionSet options) {
		setParameters((int) options.valueOf("k"), (int) options.valueOf("maxnosignal"));
	}

	public void setParameters(int k, int maxnosignalregion) {
		this.k = k;
		this.maxnosignalregion = maxnosignalregion;
	}

	public void buildDatabase() {
		if (referenceChanged) // To prevent rebuilding reference frequently
			switch (seedingmode) {
				case -1:
					if (k > 10)
						seedingmode = 1;
					else
						seedingmode = 2;
					buildDatabase();
					break;
				case 1:
					buildSortListDatabase();
					break;
				case 2:
					buildBinningDatabase();
					break;
				default:
					System.err.println("Warning! Unknown mode " + Integer.toString(seedingmode));
					buildSortListDatabase();
					break;
			}
	}

	private void buildSortListDatabase() {
		// this.databaseSeedList = convertReference(optrefmap, k, maxnosignalregion);
		if (restrictedRegions != null)
			this.databaseSeedList = convertKmerList(DataNode.getKmerWord(optrefmap, k, maxnosignalregion, restrictedRegions), k);
		else
			this.databaseSeedList = convertKmerList(DataNode.getKmerWord(optrefmap, k, maxnosignalregion), k);
	}

	private void buildBinningDatabase() {
		buildSortListDatabase();
		processFastAccessDatabase(k);
	}

	public void restrictRegion(List<GenomicPosNode> regionList) {
		if (regionList == null && restrictedRegions == null)
			referenceChanged = false;
		else {
			this.restrictedRegions = regionList;
			referenceChanged = true;
		}
	}

	private void processFastAccessDatabase(int k) {
		// generate words instead of
		fastDatabaseSeedMap = new LinkedHashMap<String, List<Kmer>>();
		List<Kmer> kmerlist = databaseSeedList.get(0);
		for (Kmer kmer : kmerlist) {
			String key = table.convertKmer(kmer);
			if (!fastDatabaseSeedMap.containsKey(key))
				fastDatabaseSeedMap.put(key, new ArrayList<Kmer>());
			fastDatabaseSeedMap.get(key).add(kmer);
		}

	}

	private List<List<Kmer>> convertKmerList(List<Kmer> kmerList, int kmerlen) {
		if (kmerList == null)
			return null;
		List<List<Kmer>> kmerlistlist = new ArrayList<List<Kmer>>(); // each sort by kmer pos 1, pos 2, pos 3...

		for (int i = 0; i < kmerlen; i++) {
			List<Kmer> dummyrefkmerlist = new ArrayList<Kmer>(kmerList);
			Collections.sort(dummyrefkmerlist, Kmer.comparator(i));
			kmerlistlist.add(dummyrefkmerlist);
		}
		return kmerlistlist;
	}

	public List<Kmer> getKmerListFromBinning(Kmer kmer, double ear, int measure) {
		List<Kmer> kmerList = new ArrayList<Kmer>();
		List<String> keylist = table.convertKmer(kmer, ear, measure);
		Kmer smallkmer = kmer.newKmer(1 - ear, measure * -1);
		Kmer largekmer = kmer.newKmer(1 + ear, measure);
		for (String key : keylist) {
			List<Kmer> kmerlist = fastDatabaseSeedMap.get(key);
			if (kmerlist != null)
				for (Kmer refkmer : kmerlist) {
					boolean wrong = false;
					for (int i = 0; i < refkmer.k(); i++)
						if (refkmer.get(i) < smallkmer.get(i) || refkmer.get(i) > largekmer.get(i)) {
							wrong = true;
							break;
						}
					if (!wrong) {
						if (refkmer.limitRange(kmer, measure, ear))
							kmerList.add(refkmer);
					}
				}
		}
		return kmerList;

	}

	public List<Kmer> getKmerListFromSortList(Kmer kmer, double ear, int measure) {
		List<Kmer> matchedkmerlist = null;
		Kmer smallkmer = kmer.newKmer(1 - ear, measure * -1);
		Kmer largekmer = kmer.newKmer(1 + ear, measure);
		// long smalllen = smallkmer.length();
		// long largelen = largekmer.length();
		// HashSet<Kmer> kmerHash = new HashSet<Kmer>();
		for (int i = 0; i < kmer.k(); i++) {
			List<Kmer> refkmerlist = databaseSeedList.get(i);
			int startpos = Collections.binarySearch(refkmerlist, smallkmer, Kmer.comparator(i));
			if (startpos < 0)
				startpos = (startpos + 1) * -1;
			else {
				while (startpos >= 0) {
					if (smallkmer.compare(refkmerlist.get(startpos), i) > 0)
						break;
					startpos--;
				}
				startpos++;
			}
			int stoppos = Collections.binarySearch(refkmerlist, largekmer, Kmer.comparator(i));
			if (stoppos < 0)
				stoppos = (stoppos + 1) * -1;
			else {
				while (stoppos < refkmerlist.size()) {
					if (largekmer.compare(refkmerlist.get(stoppos), i) < 0)
						break;
					stoppos++;
				}
			}
			stoppos--;

			List<Kmer> sublist = refkmerlist.subList(startpos, stoppos + 1);
			if (matchedkmerlist == null) {
				matchedkmerlist = new ArrayList<Kmer>(sublist);
				// matchedkmerlist = new ArrayList<Kmer>();
				// for (Kmer km : sublist)
				// {
				// if (km.len <= largelen && km.len >= smalllen)
				// matchedkmerlist.add(km);
				// }
			} else {
				List<Kmer> newlist = new ArrayList<Kmer>();
				if (matchedkmerlist.size() > sublist.size()) {
					Set<Kmer> subset = new HashSet<Kmer>(sublist);
					for (Kmer candidatekmer : matchedkmerlist)
						if (subset.contains(candidatekmer))
							newlist.add(candidatekmer);
				} else {
					Set<Kmer> subset = new HashSet<Kmer>(matchedkmerlist);
					for (Kmer candidatekmer : sublist) {
						if (subset.contains(candidatekmer))
							newlist.add(candidatekmer);
					}
				}
				matchedkmerlist = newlist;
			}

		}
		List<Kmer> kmerList = new ArrayList<Kmer>();
		for (Kmer matchedkmer : matchedkmerlist) {
			if (matchedkmer.limitRange(kmer, measure, ear))
				kmerList.add(matchedkmer);
		}
		return kmerList;

	}

	public List<Kmer> getKmerList(Kmer kmer, double ear, int measure) {
		switch (seedingmode) {
			case -1:
				if (k > 10)
					seedingmode = 1;
				else
					seedingmode = 2;
				return getKmerList(kmer, ear, measure);
			case 1:
				return getKmerListFromSortList(kmer, ear, measure);
			case 2:
				return getKmerListFromBinning(kmer, ear, measure);
			default:
				System.err.println("Warning! Unknown mode " + Integer.toString(seedingmode));
				return getKmerListFromSortList(kmer, ear, measure);
		}

	}

	public List<Seed> getSeed(Kmer kmer, double ear, int measure) {

		List<Kmer> matchedkmerlist = getKmerList(kmer, ear, measure);
		List<Seed> seedlist = new ArrayList<Seed>();
		for (Kmer matchedkmer : matchedkmerlist) {

			Seed s = new Seed(matchedkmer, new Kmer(kmer));
			if (s.limitRange(measure, ear))
				seedlist.add(s);
		}
		return seedlist;

	}

	public List<Seed> getJoinedSeed(Kmer kmer, double ear, int measure) {
		List<Seed> seedlist = getSeed(kmer, ear, measure);
		return seedJoin(seedlist);
	}

	private List<Seed> seedJoin(List<Seed> pooledseedlist) {
		Collections.sort(pooledseedlist, Kmer.comparatorSourcePos());
		List<Seed> joinedseedlist = new ArrayList<Seed>();
		for (int i = pooledseedlist.size() - 1; i >= 0; i--) {
			int j = i;
			while ((pooledseedlist.get(j).pos - pooledseedlist.get(i).pos <= 2) && (pooledseedlist.get(i).source.equalsIgnoreCase(pooledseedlist.get(j).source))) {
				Seed previousseed = pooledseedlist.get(j);
				Seed recentseed = pooledseedlist.get(i);
				if (previousseed.pos - recentseed.pos == 1) {
					if (previousseed.kmerpointer.pos - recentseed.kmerpointer.pos == 1) {
						// new
						double ubound = recentseed.rangeUBound;
						double lbound = recentseed.rangeLBound;
						if (ubound > previousseed.rangeUBound)
							ubound = previousseed.rangeUBound;
						if (lbound < previousseed.rangeLBound)
							lbound = previousseed.rangeLBound;
						if (ubound >= lbound) {
							recentseed.sizelist.addAll(previousseed.sizelist.subList(recentseed.k() - 1, previousseed.k()));
							recentseed.kmerpointer.sizelist.addAll(previousseed.kmerpointer.sizelist.subList(recentseed.kmerpointer.k() - 1, previousseed.kmerpointer.k()));
							recentseed.rangeLBound = lbound;
							recentseed.rangeUBound = ubound;
							pooledseedlist.set(j, null);

							break; // only one can be joined, must be only one
							// joined;
						}
					}
				}
				do {
					j++;
					if (j == pooledseedlist.size())
						break;
				} while (pooledseedlist.get(j) == null);
				if (j == pooledseedlist.size())
					break;

			}
		}
		for (Seed seed : pooledseedlist)
			if (seed != null)
				joinedseedlist.add(seed);
		return joinedseedlist;

	}

	// filter those similar kmers
	public List<Kmer> filter(List<Kmer> fragmentkmerlist, double ear, int measure, int maxSeedNumber, int maxSignalConsidered) {
		List<Kmer> filteredKmerList = new ArrayList<Kmer>();
		if (maxSignalConsidered == -1)
			maxSignalConsidered = Integer.MAX_VALUE;
		for (Kmer kmer : fragmentkmerlist) {
			List<Seed> seedList = this.getSeed(kmer, ear, measure);
			int r = 0;
			for (Seed seed : seedList)
				if (kmer.source.equalsIgnoreCase(seed.source))
					if (kmer.pos > seed.pos) {
						if (kmer.pos - seed.pos - seed.k() <= maxSignalConsidered)
							r++;
					} else if (seed.pos - kmer.pos - kmer.k() <= maxSignalConsidered)
						r++;

			r++;
			if (r <= maxSeedNumber)
				filteredKmerList.add(kmer);
		}
		return filteredKmerList;
	}

	public SeedDatabase copy() {
		SeedDatabase seedDatabase = new SeedDatabase(optrefmap, databaseSeedList);
		seedDatabase.setMode(seedingmode);
		seedDatabase.setParameters(k, maxnosignalregion);
		seedDatabase.fastDatabaseSeedMap = this.fastDatabaseSeedMap;
		seedDatabase.table = this.table;
		seedDatabase.referenceChanged = false;
		return seedDatabase;
	}

	public static void assignOptions(ExtendOptionParser parser, int level) {
		parser.addHeader("Seeding Options", level);
		parser.accepts("seedingmode", "Seeding mode: 1: Opt for long k-mer; 2: Opt for short k-mer; -1: Auto-selection").withRequiredArg().ofType(Integer.class).defaultsTo(-1);
		parser.accepts("k", "Kmer length.").withOptionalArg().ofType(Integer.class).defaultsTo(3);
		parser.accepts("maxnosignal", "Maximum no signal region between signals for seeding.").withOptionalArg().ofType(Integer.class).defaultsTo(10000000);

	}
}

class FastConversionTable {
	private List<Integer> sizelist; // final size must be Integer.max
	private List<Character> rclist;

	public FastConversionTable(List<Integer> sizelist, List<Character> rclist) {
		this.sizelist = sizelist;
		this.rclist = rclist;
		if (rclist.size() - sizelist.size() != 0)
			System.err.println("Table corrupted. List Size Inappropriate. Proceed anyway.");
	}

	public String convertKmer(Kmer kmer) {
		StringBuilder s = new StringBuilder();
		for (long possize : kmer.sizelist) {
			int pt = 0;
			while (possize > sizelist.get(pt))
				pt++;
			s.append(rclist.get(pt));
		}
		return s.toString();
	}

	private List<String> join(List<String> slist, List<Character> clist) {
		// if (slist.size() != clist.size())
		// return null;

		List<String> newslist = new ArrayList<String>();
		for (String s : slist)
			for (char c : clist)
				newslist.add(s + c);
		return newslist;
	}

	public List<String> convertKmer(Kmer kmer, double ear, int measure) {
		Kmer smallkmer = kmer.newKmer(1 - ear, measure * -1);
		Kmer largekmer = kmer.newKmer(1 + ear, measure);
		String smallkey = convertKmer(smallkmer);
		String largekey = convertKmer(largekmer);
		List<String> recentList = new ArrayList<String>();
		recentList.add("");
		for (int i = 0; i < smallkey.length(); i++) {
			List<Character> clist = new ArrayList<Character>();
			for (char c = smallkey.charAt(i); c <= largekey.charAt(i); c++)
				clist.add(c);
			recentList = join(recentList, clist);
		}
		return recentList;
	}

	public static FastConversionTable standardTable() {
		int gap = 5000;
		List<Character> clist = new ArrayList<Character>();
		List<Integer> sizelist = new ArrayList<Integer>();
		int recent = 0;
		char next = 1;
		for (int i = 1; i < 255; i++) {
			recent += gap;
			sizelist.add(recent);
			clist.add(next);
			next++;
		}
		sizelist.add(Integer.MAX_VALUE);
		clist.add(next);
		return new FastConversionTable(sizelist, clist);
	}
}

// class FastKmerList
// {
// public List<Object> kmerlist;
// public FastKmerList(List<Object> kmerlist) {
// this.kmerlist = kmerlist;
// }
// public String getKmer(String s)
// {
// s.substring(1, s.length());
// }
// private List<Kmer> get(String s, List<Object> next)
// {
// if (s.length() == 1)
// if (next.get(s.charAt(0))
// }

// }
// private List<Seed> limitSeed(List<Seed> seedlist, double ear, int measure)
// {
// List<Seed> newseedlist = new ArrayList<Seed>();
//
// for (Seed seed : seedlist)
// {
// double ubound = 1 + ear;
// double lbound = 1 - ear;
// for (int pos = 0; pos < seed.k(); pos++)
// {
// double newubound = (seed.kmerpointer.get(pos) + measure) / (double) seed.get(pos);
// double newlbound = (seed.kmerpointer.get(pos) - measure) / (double) seed.get(pos);
// if (newubound < ubound)
// ubound = newubound;
// if (newlbound > lbound)
// lbound = newlbound;
// }
//
// }
// return newseedlist;
// }
