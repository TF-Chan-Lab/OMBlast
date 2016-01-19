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


package aldenjava.opticalmapping.data.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.lang.ArrayUtils;

import aldenjava.opticalmapping.GenomicPosNode;
import aldenjava.opticalmapping.data.Identifiable;
import aldenjava.opticalmapping.mapper.seeding.Kmer;

/**
 * Basic Data Node. Two schemes have been used for an optical map record. (1) Size (bp) of molecule and positions of labeling/segment break; and (2) ordered sizes of segments Here (1) is selected for easy segment size deduction, as represented by <code>refp</code>. To allow convenient segment size as in scheme (2), one can use the function <code>getRefl</code>. It is a 1-based coordinate system. Each position occupies 1bp.
 * <p>
 * Example: A molecule ("|": signal) ----|------|--|-----
 * <li>(1) Size: 20; refp: [5, 12, 15]</li>
 * <li>(2) refl: [4, 6, 2, 5]</li>
 * 
 * @author Alden
 *
 */

// Require better encapsulation
public class DataNode implements Identifiable<String> {

	public String name = null;
	public long size = -1;
	public long refp[] = null;
	public boolean circular = false; // Temporarily not used
	public SimulationInfo simuInfo = null;

	/**
	 * Creates a blank <code>DataNode</code> without any information.
	 */
	public DataNode() {
	}

	/**
	 * Constructs a new <code>DataNode</code>, initialized to match the values of the specified <code>DataNode</code>.
	 * 
	 * @param data
	 *            the <code>DataNode</code> from which to copy initial values to a newly constructed <code>DataNode</code>
	 */
	public DataNode(DataNode data) {
		this.name = data.name;
		this.size = data.size;
		this.refp = Arrays.copyOf(data.refp, data.refp.length);
		;
		this.importSimulationInfo(data.simuInfo);
	}

	/**
	 * Constructs a new <code>DataNode</code> whose name, size, label positions are specified by the arguments
	 * 
	 * @param name
	 *            the specified name for this Data
	 * @param size
	 *            the specified size for this Data
	 * @param refp
	 *            the specified label positions
	 */
	public DataNode(String name, long size, long[] refp) {
		this.name = name;
		this.size = size;
		this.refp = refp;
	}

	/**
	 * 
	 * @param name
	 *            the specified name for this Data
	 * @param refl
	 *            the size of ordered segments
	 */
	public DataNode(String name, long[] refl) {
		this.name = name;
		long[] refp = new long[refl.length - 1];
		long size = 0;
		for (int i = 0; i < refl.length - 1; i++) {
			size += refl[i] + 1; // need +1 here? // Yes
			refp[i] = size; // Need -1 here?
		}
		size += refl[refl.length - 1]; // need +1 here?
		this.refp = refp;
		this.size = size;
	}

	/**
	 * Constructs a new <code>DataNode</code> based on a boolean array, where {@code true} represents existence of a signal
	 * 
	 * @param name
	 *            the specified name for this Data
	 * @param zeroOneSeq
	 *            the specified boolean array to indicate the label positions
	 */
	public DataNode(String name, boolean[] zeroOneSeq) {
		long size = zeroOneSeq.length;
		List<Long> refp = new ArrayList<Long>();
		for (int i = 0; i < zeroOneSeq.length; i++) {
			if (zeroOneSeq[i])
				refp.add((long) i + 1);
		}
		this.name = name;
		this.size = size;
		this.refp = ArrayUtils.toPrimitive(refp.toArray(new Long[refp.size()]));
		this.simuInfo = null;
	}

	/**
	 * Constructs a new <code>DataNode</code> without any labels
	 * 
	 * @param name
	 *            the specified name for this Data
	 * @param size
	 *            the specified size for this Data
	 */
	public DataNode(String name, long size) {
		this(name, size, new long[0]);
	}

	public void importSimulationInfo(GenomicPosNode simuRegion, int simuStrand) {
		this.simuInfo = new SimulationInfo(simuRegion, simuStrand);
	}

	public void importSimulationInfo(SimulationInfo simuInfo) {
		if (simuInfo == null)
			this.simuInfo = null;
		else
			this.simuInfo = new SimulationInfo(simuInfo);
	}

	// Information
	@Override
	public String getIdentifier() {
		return name;
	}

	public GenomicPosNode getGenomicPos(int subrefstart, int subrefstop, boolean flankingSig) {
		return new GenomicPosNode(this.name, subrefstart > 0 ? refp[subrefstart - 1] + (flankingSig ? 0 : 1) : 1, subrefstop < refp.length ? refp[subrefstop] - (flankingSig ? 0 : 1) : size);
	}

	/**
	 * Returns the segment length excluding the length of signals
	 * 
	 * @param index
	 * @return The segment length
	 */
	public long getRefl(int index) {
		if (index > refp.length || index < 0)
			throw new IndexOutOfBoundsException("Segment index is out of bound");
		if (index == 0 && refp.length == 0)
			return size;
		if (index == refp.length)
			return size - refp[refp.length - 1];
		if (index == 0)
			return refp[index] - 1;
		return refp[index] - refp[index - 1] - 1;

	}

	/**
	 * Returns the signal density of this <code>DataNode</code>
	 * 
	 * @return the signal density along the <code>DataNode</code>
	 */
	public double getSignalDensity() {
		return ((double) refp.length / (double) size);
	}

	/**
	 * Returns the total signals of this <code>DataNode</code>. The results should be equal to <code>getTotalSegment</code> - 1.
	 * 
	 * @return Number of elements in <code>refp</code>
	 * @see #getTotalSegment
	 */
	public int getTotalSignal() {
		return refp.length;
	}

	/**
	 * Returns the total segments of this <code>DataNode</code>. The results should be equal to <code>getTotalSignal</code> + 1.
	 * 
	 * @return Number of elements in <code>refp</code>
	 * @see #getTotalSignal
	 */
	public int getTotalSegment() {
		return refp.length + 1;
	}

	/**
	 * Returns the signal density along <code>DataNode</code> using a sliding window of size {@code regionSize}.
	 * 
	 * @param regionSize
	 *            size of the sliding window
	 * @return An array of signal density
	 */
	public double[] getSignalDensityAcrossRegion(long regionSize) {
		int startSig = 0;
		int stopSig = 0;
		List<Double> densityList = new ArrayList<Double>();
		while (startSig < refp.length - 1 && stopSig < refp.length - 1) {
			if (getSignalLength(startSig, stopSig) + 2 > regionSize) {
				startSig++;
				continue;
			}

			if (getSignalLength(startSig, stopSig) + 2 + getRefl(startSig) + getRefl(stopSig + 1) < regionSize) {
				stopSig++;
				continue;
			}

			// Must be in this way
			double density = (stopSig - startSig + 1) / (double) regionSize;
			densityList.add(density);
			startSig++;
			stopSig++;
		}
		return ArrayUtils.toPrimitive(densityList.toArray(new Double[densityList.size()]));
	}

	/**
	 * Indicates if this data is in circular form. This function is not yet implemented
	 * 
	 * @return Always return <code>false</code>
	 */
	public boolean isCircular() {
		return circular;
	}

	/**
	 * Fixes a molecule using 4 steps <li>If data with zero size has any label, all labels are erased <li>Label positions are changed to 1 or <code>size</code> if they locate at position smaller than 1 or larger than <code>size</code> respectively. This fix can lead to overlapped signals <li>Label positions are sorted in ascending order <li>Overlapped labels are separated to consecutive signals. This shifts all followed signals and augments the size
	 * 
	 * @return <code>true</code> if the molecule has undergone modification
	 */
	public boolean fix() {
		boolean fixed = false;
		// 4 fixing procedures

		// 1. Fix zero size case
		if (size == 0 && refp.length > 0) {
			refp = new long[0];
			fixed = true;
		}

		// 2. Correct invalid refp position
		for (int i = 0; i < refp.length; i++) {
			if (refp[i] <= 0) {
				refp[i] = 1;
				fixed = true;
			}
			if (refp[i] > size) {
				refp[i] = size;
				fixed = true;
			}
		}

		// 3. Sort signals in wrong order
		long prevP = 0;
		boolean requiredSorting = false;
		for (long p : refp) {
			if (p < prevP) {
				requiredSorting = true;
				break;
			} else
				prevP = p;
		}
		if (requiredSorting) {
			Arrays.sort(refp);
			fixed = true;
		}

		// 4. Allocate space for overlapped signals, size will be augmented by 1 for each allocated space
		long shift = 0;
		long prev = -1;
		for (int i = 0; i < refp.length; i++) {
			if (refp[i] <= prev)
				shift++;
			prev = refp[i];
			refp[i] += shift;
		}
		if (shift > 0) {
			size += shift;
			fixed = true;
		}
		return fixed;
	}

	/**
	 * Returns the length of the molecule using the start and stop signal as index, excluding the length of flanking signals
	 * 
	 * @param sigStart
	 * @param sigStop
	 * @return the specified length
	 * @see #length(int, int)
	 */
	public long getSignalLength(int sigStart, int sigStop) {
		return length(sigStart + 1, sigStop);
	}

	/**
	 * Return the size of the molecule
	 * 
	 * @return size of the molecule
	 */
	public long length() {
		return size;
	}

	/**
	 * Returns the length of the molecule within start and stop segment, excluding the flanking signals
	 * 
	 * @param start
	 *            The start segment, inclusive
	 * @param end
	 *            The stop segment, inclusive
	 * @return the specified length
	 */
	public long length(int start, int end) {
		if (start < 0 || start > refp.length || end < 0 || end > refp.length)
			throw new IndexOutOfBoundsException("Index range is out of bound");
		return ((end == refp.length ? (size + 1) : refp[end]) - (start == 0 ? 0 : refp[start - 1]) - 1);
	}

	public String getReflString() {
		return getReflString(';');
	}

	public String getReflString(char separator) {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < refp.length; i++) {
			s.append(Long.toString(getRefl(i)));
			s.append(separator);
		}
		s.append(Long.toString(getRefl(refp.length)));
		return s.toString();
	}

	public String getRefpString() {
		return getRefpString(';');
	}

	public String getRefpString(char separator) {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < refp.length - 1; i++) {
			s.append(Long.toString(refp[i]));
			s.append(separator);
		}
		if (refp.length > 0)
			s.append(Long.toString(refp[refp.length - 1]));
		return s.toString();
	}

	@Override
	public String toString() {
		return (name);
	}

	/**
	 * Check if simulation information is available
	 * 
	 * @return <code>true</code> if simuInfo is not <code>null</code>
	 */
	public boolean hasSimulationInfo() {
		return simuInfo != null;
	}

	/**
	 * Check if there is a signal at the specified position
	 * 
	 * @param pos
	 *            Specified position
	 * @return <code>true</code> if the signal exists
	 */
	public boolean hasSignal(long pos) {
		return (Arrays.binarySearch(refp, pos) >= 0);
	}

	/**
	 * Search the refp index of the signal at the specified position using the binary search algorithm
	 * 
	 * @param pos
	 *            Specified position
	 * @return the refp index of the signal at the specified position, or <code>-1</code> if the signal does not exist
	 */
	public int findExactRefpIndex(long pos) {
		int index = findRefpIndex(pos);
		if (index >= refp.length)
			return -1;
		if (refp[index] != pos)
			return -1;
		return index;
	}

	/**
	 * Search the refp index signals using the binary search algorithm
	 * 
	 * @param pos
	 *            Specified position
	 * @return the refp index of the signal at the specified position, or the refp index of the <code>insertion point</code> if the signal does not exist
	 * @see findExactRefpIndex(long)
	 */
	public int findRefpIndex(long pos) {
		int index = Arrays.binarySearch(refp, pos); // index points to the last refp just larger than or equal to start
		if (index < 0) {
			index = (index + 1) * -1;
		} else {
			while (index >= 0 && refp[index] == pos)
				index--;
			if (index < 0)
				index = 0;
			else
				index++;
		}
		return index;
	}

	/**
	 * Get a new reversed data. Simulation information is not reversed.
	 * 
	 * @return new reversed <code>DataNode</code>
	 */
	public DataNode getReverse() {
		long[] newrefp = new long[refp.length];
		for (int i = 0; i < refp.length; i++)
			newrefp[refp.length - i - 1] = size - refp[i] + 1;
		DataNode data = new DataNode(name, size, newrefp);
		if (hasSimulationInfo())
			data.importSimulationInfo(simuInfo);
		return data;
	}

	public DataNode subRefNode(String newname, long start, long stop) {
		long length = stop - start + 1;
		int index = findRefpIndex(start);
		int index2 = findRefpIndex(stop);
		if (!hasSignal(stop))
			index2--;
		long[] subrefp = new long[index2 - index + 1];
		System.arraycopy(refp, index, subrefp, 0, index2 - index + 1);
		for (int i = 0; i < subrefp.length; i++)
			subrefp[i] = subrefp[i] - start + 1;
		return new DataNode(newname, length, subrefp);
	}

	public DataNode subRefNode(long start, long stop) {
		return subRefNode(name, start, stop);
	}

	public DataNode subRefNode(GenomicPosNode region) {
		if (!region.ref.equalsIgnoreCase(name)) {
			System.err.println("Region doesn't match!");
			return null;
		}
		return subRefNode(region.start, region.stop);
	}

	public DataNode subRefNode(String newname, int subrefstart, int subrefstop, boolean flankingSignal) {
		// long start;
		// long stop;
		// if (subrefstart == 0)
		// start = 1;
		// else
		// {
		// start = refp[subrefstart - 1];
		// if (!flankingSignal)
		// start++;
		// }
		// if (subrefstop == refp.length)
		// stop = size;
		// else
		// {
		// stop = refp[subrefstop];
		// if (!flankingSignal)
		// stop--;
		// }
		// return subRefNode(newname, start, stop);
		GenomicPosNode pos = getGenomicPos(subrefstart, subrefstop, flankingSignal);
		return subRefNode(newname, pos.start, pos.stop);
	}

	public DataNode subRefNode(int subrefstart, int subrefstop, boolean flankingSignal) {
		return subRefNode(name, subrefstart, subrefstop, flankingSignal);
	}

	public List<Kmer> getKmerWord(int kmerlen, long maxnosignalregion) {
		List<Kmer> kmerlist = new ArrayList<Kmer>();
		for (int i = 1; i < getTotalSegment() - kmerlen; i++) // +1 is no need to be added: 2014/05/27
		{
			List<Long> sizelist = new ArrayList<Long>();
			boolean discarded = false;
			for (int j = i; j < i + kmerlen; j++)
				if (getRefl(j) <= maxnosignalregion)
					sizelist.add(getRefl(j));
				else
					discarded = true;
			if (!discarded)
				kmerlist.add(new Kmer(name, i, sizelist));
		}

		return kmerlist;
	}

	public List<Kmer> getKmerWord(int kmerlen, long maxnosignalregion, List<GenomicPosNode> restrictedRegions) {
		List<GenomicPosNode> regionList = new ArrayList<GenomicPosNode>();
		for (GenomicPosNode region : restrictedRegions) {
			if (region.ref.equals(this.name))
				regionList.add(region);
		}

		List<Kmer> kmerlist = new ArrayList<Kmer>();

		for (GenomicPosNode region : regionList) {
			int start = this.findRefpIndex(region.start);
			int stop = this.findRefpIndex(region.stop);

			if (stop >= refp.length)
				stop--;
			else if (refp[stop] > region.stop)
				stop--;
			NEXTKmer: for (int i = start + 1; i <= stop - kmerlen + 1; i++) {
				List<Long> sizelist = new ArrayList<Long>();
				// boolean discarded = false;
				// for (int j = i - 1; j < i + kmerlen; j++) {
				// boolean pass = false;
				// for (GenomicPosNode region : regionList)
				// if (region.getLoc().contain(refp[j]))
				// pass = true;
				// if (!pass)
				// // discarded = true;
				// continue NEXTKmer;
				// }
				for (int j = i; j < i + kmerlen; j++) {
					if (getRefl(j) <= maxnosignalregion)
						sizelist.add(getRefl(j));
					else
						// discarded = true;
						continue NEXTKmer;
				}
				// if (!discarded)
				kmerlist.add(new Kmer(name, i, sizelist));
			}
		}
		return kmerlist;
	}

	public List<Kmer> getErrorKmerWord(int kmerlen, int maxnosignalregion, int errorno) {
		List<Kmer> kmerlist = new ArrayList<Kmer>();
		// Generate error-free kmer list
		int finalkmerlen = kmerlen + errorno;
		for (int i = 1; i < getTotalSegment() - finalkmerlen; i++) {
			List<Long> sizelist = new ArrayList<Long>();
			boolean discarded = false;
			for (int j = i; j < i + finalkmerlen; j++)
				if (getRefl(j) <= maxnosignalregion)
					sizelist.add(getRefl(j));
				else
					discarded = true;
			if (!discarded)
				kmerlist.add(new Kmer(name, i, sizelist));
		}
		// Induce error
		if (errorno == 0)
			return kmerlist;
		else {
			if (errorno > 1) {
				System.err.println("Only support error number <= 1.");
				return null;
			}
			List<Kmer> errorkmerlist = new ArrayList<Kmer>();
			// int pos = 1;
			// int error = 0;
			// List<Integer> errorposlist = new ArrayList<Integer>();
			// for (Kmer kmer : kmerlist)
			// {
			// while (true)
			// {
			// if (error == errorno)
			// {
			// Kmer k = new Kmer(kmer.source, kmer.pos, kmer.sizelist, errorposlist);
			// errorkmerlist.add();
			// errorposlist.add(pos);
			// }
			//
			// }
			// }
			// temporarily used for errorno = 1
			for (Kmer kmer : kmerlist) {
				for (int i = 0; i < finalkmerlen - 1; i++) {
					List<Integer> errorposlist = new ArrayList<Integer>();
					errorposlist.add(i);
					errorkmerlist.add(new Kmer(kmer.source, kmer.pos, kmer.sizelist, errorposlist));
				}
			}
			return errorkmerlist;
		}
	}

	// Modify
	public void insertSignal(long pos) {
		int index = Arrays.binarySearch(refp, pos);
		if (index < 0) {
			long[] newrefp = new long[refp.length + 1];
			int insertionpoint = (index + 1) * -1;
			System.arraycopy(refp, 0, newrefp, 0, insertionpoint);
			newrefp[insertionpoint] = pos;
			System.arraycopy(refp, insertionpoint, newrefp, insertionpoint + 1, refp.length - insertionpoint);
			refp = newrefp;
			// this.rebuildRefl();
		} else
			System.err.println("Warning: Signal already exists at " + Long.toString(pos));
	}

	public void removeSignal(long pos) {
		int index = Arrays.binarySearch(refp, pos);
		if (index < 0)
			System.err.println("Warning: Signal does not exist at " + Long.toString(pos));
		else {
			long[] newrefp = new long[refp.length - 1];
			System.arraycopy(refp, 0, newrefp, 0, index);
			System.arraycopy(refp, index + 1, newrefp, index, refp.length - index - 1);
			refp = newrefp;
			// this.rebuildRefl();
		}
	}

	public void removeRefpSig(int refpPos) {
		long[] newrefp = new long[refp.length - 1];
		System.arraycopy(refp, 0, newrefp, 0, refpPos);
		System.arraycopy(refp, refpPos + 1, newrefp, refpPos, refp.length - refpPos - 1);
		refp = newrefp;
	}

	public void insert(long start, DataNode ref) {
		int index = findRefpIndex(start);
		long[] newrefp = new long[refp.length + ref.refp.length];
		System.arraycopy(refp, 0, newrefp, 0, index);
		System.arraycopy(ref.refp, 0, newrefp, index, ref.refp.length);
		System.arraycopy(refp, index, newrefp, index + ref.refp.length, refp.length - index);

		for (int i = index; i < index + ref.refp.length; i++)
			newrefp[i] += start;
		for (int i = index + ref.refp.length; i < newrefp.length; i++)
			newrefp[i] += ref.length();
		refp = newrefp;
		size += ref.size;
	}

	public DataNode remove(long start, long length) {
		if (start < 1)
			throw new IllegalArgumentException("Illegal start: " + start);
		if (start > size)
			throw new IllegalArgumentException("Start " + start + " is larger than data size " + size);
		if (length < 0)
			throw new IllegalArgumentException("Illegal length: " + length);
		if (start + length - 1 > size)
			throw new IllegalArgumentException("The target region " + (start + length - 1) + " is out of range to " + size);

		if (length == 0)
			return new DataNode(name, 0);

		int index = Arrays.binarySearch(refp, start); // refp[index] >= start
		if (index < 0)
			index = (index + 1) * -1;
		else {
			while (index >= 0 && refp[index] == start)
				index--;
			if (index < 0)
				index = 0;
			else
				index++;
		}

		long stop = start + length - 1;
		int index2 = Arrays.binarySearch(refp, stop); // refp[index2] >= stop
		if (index2 < 0)
			index2 = (index2 + 1) * -1;
		else {
			while (index2 >= 0 && refp[index2] == stop)
				index2--;
			if (index2 < 0)
				index2 = 0;
			else
				index2++;
		}

		long[] removedrefp = new long[index2 - index];
		long[] newrefp = new long[refp.length - removedrefp.length];
		System.arraycopy(refp, 0, newrefp, 0, index);
		System.arraycopy(refp, index2, newrefp, index, refp.length - index2);
		for (int i = index; i < newrefp.length; i++)
			newrefp[i] -= length;
		System.arraycopy(refp, index, removedrefp, 0, index2 - index);
		for (int i = 0; i < removedrefp.length; i++)
			removedrefp[i] -= start;
		refp = newrefp;
		size -= length;
		return new DataNode(name, length, removedrefp);
	}

	public DataNode remove(long start) {
		// Remove till the end
		return remove(start, size - start + 1);
	}

	public void trim(long left, long right) {
		if (left + right > size) {
			size = 0;
			refp = new long[0];
		}

		remove(size - right + 1, right);
		remove(1, left);
	}

	public void replace(long start, long orgLength, DataNode ref) {
		this.remove(start, orgLength);
		this.insert(start, ref);
	}

	/**
	 * Join the refs in order into current data
	 * 
	 * @param refs
	 *            the data to be joined to current data
	 */
	public void join(DataNode... refs) {
		int newrefpLength = refp.length;
		for (int i = 0; i < refs.length; i++)
			newrefpLength += refs[i].refp.length;

		long[] newrefp = new long[newrefpLength];
		System.arraycopy(refp, 0, newrefp, 0, refp.length);
		long newsize = size;
		int startpos = refp.length;

		for (int i = 0; i < refs.length; i++) {
			DataNode ref = refs[i];
			System.arraycopy(ref.refp, 0, newrefp, startpos, ref.refp.length);
			for (int j = startpos; j < ref.refp.length + startpos; j++)
				newrefp[j] += newsize;
			newsize += ref.size;
			startpos += ref.refp.length;

		}
		refp = newrefp;
		size = newsize;
	}

	/**
	 * Reverses a portion of data as indicated by the start position and length
	 * 
	 * @param start
	 *            the starting position of target region
	 * @param length
	 *            the length of target region
	 */
	public void reverse(long start, long length) {
		if (start < 1)
			throw new IllegalArgumentException("Illegal start: " + start);
		if (start > size)
			throw new IllegalArgumentException("Start " + start + " is larger than data size " + size);
		if (length < 0)
			throw new IllegalArgumentException("Illegal length: " + length);
		if (start + length - 1 > size)
			throw new IllegalArgumentException("The target region " + (start + length - 1) + " is out of range to " + size);
		DataNode invertRef = this.remove(start, length);
		invertRef.reverse();
		this.insert(start, invertRef);
	}

	/**
	 * Reverses the current data
	 */
	public void reverse() {
		long[] newrefp = new long[refp.length];
		for (int i = 0; i < refp.length; i++)
			newrefp[i] = size - refp[refp.length - i - 1] + 1;
		refp = newrefp;
	}

	/**
	 * Scales the data according to <code>ratio</code>.
	 * <p>
	 * This method does not guarantee checking on overlapped signals;
	 * <p>
	 * The following is not guaranteed to result in no change in data
	 * 
	 * <pre>
	 * data.scale(ratio).scale(1 / ratio)
	 * </pre>
	 * 
	 * @param ratio
	 *            the scaling factor
	 */
	public void scale(double ratio) {
		for (int i = 0; i < refp.length; i++)
			refp[i] = (long) (refp[i] * ratio + 0.5);
		size = (long) (size * ratio + 0.5);
	}

	/**
	 * 
	 * @param degeneracy
	 */
	public void degenerate(int degeneracy) {
		if (degeneracy < 0)
			throw new IllegalArgumentException("Negative degeneracy is not accepted " + degeneracy);
		if (degeneracy == 0) // Nothing needed to be done
			return;
		if (getTotalSignal() <= 1) // Nothing needed to be done
			return;
		List<Long> newRefpList = new ArrayList<Long>();
		long sum = 0;
		int element = 0;
		long prev = -1;
		for (int i = 0; i < refp.length; i++) {
			if (prev != -1) {
				if (refp[i] - prev > degeneracy) {
					newRefpList.add((long) (sum / element));
					sum = 0;
					element = 0;
				}
			}
			sum += refp[i];
			element++;
			prev = refp[i];
		}
		newRefpList.add((long) (sum / element));
		refp = ArrayUtils.toPrimitive(newRefpList.toArray(new Long[newRefpList.size()]));
	}

	public static int getTotalSignal(LinkedHashMap<String, DataNode> optrefmap) {
		int total = 0;
		for (DataNode data : optrefmap.values())
			total += data.getTotalSignal();
		return total;
	}

	public static long getTotalSize(LinkedHashMap<String, DataNode> optrefmap) {
		long total = 0;
		for (DataNode item : optrefmap.values()) {
			total += item.size;
		}
		return total;
	}

	public static double getDensity(LinkedHashMap<String, DataNode> optrefmap) {
		return getTotalSignal(optrefmap) / (double) getTotalSize(optrefmap);
	}

	public static long[] parseReflInString(String s, String separator) {
		if (s.trim().isEmpty())
			return new long[0];
		if (s.endsWith(separator)) // Dealing with wrong output with last ";"
			s = s.substring(0, s.length() - 1);
		String[] f = s.trim().split(separator);
		long[] fragment = new long[f.length];
		for (int i = 0; i < f.length; i++)
			fragment[i] = Long.parseLong(f[i]);
		return fragment;
	}

	public static List<Kmer> getKmerWord(LinkedHashMap<String, DataNode> optrefmap, int kmerlen, long maxnosignalregion) {
		if (optrefmap == null)
			throw new NullPointerException("optrefmap");
		List<Kmer> kmerList = new ArrayList<Kmer>();
		for (DataNode ref : optrefmap.values()) {
			kmerList.addAll(ref.getKmerWord(kmerlen, maxnosignalregion));
		}
		return kmerList;
	}

	public static List<Kmer> getKmerWord(LinkedHashMap<String, DataNode> optrefmap, int kmerlen, long maxnosignalregion, List<GenomicPosNode> restrictedRegions) {
		List<Kmer> kmerList = new ArrayList<Kmer>();
		for (DataNode ref : optrefmap.values()) {
			kmerList.addAll(ref.getKmerWord(kmerlen, maxnosignalregion, restrictedRegions));
		}
		return kmerList;
	}

	public static Comparator<DataNode> sizecomparator = new Comparator<DataNode>() {
		@Override
		public int compare(DataNode d1, DataNode d2) {
			return Long.valueOf(d1.size).compareTo(Long.valueOf(d2.size));
		}
	};
	public static Comparator<DataNode> signalcomparator = new Comparator<DataNode>() {
		@Override
		public int compare(DataNode d1, DataNode d2) {
			return Long.valueOf(d1.getTotalSignal()).compareTo(Long.valueOf(d2.getTotalSignal()));
		}
	};

}