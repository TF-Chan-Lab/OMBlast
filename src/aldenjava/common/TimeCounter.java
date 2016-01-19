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
import java.util.List;

/**
 * A class to count the running time by System.currentTimeMillis()
 * 
 * @author Alden
 *
 */
public class TimeCounter {

	class TimeItem {
		long totaltime;
		List<Long> starttime;
		String header;

		public TimeItem() {
			totaltime = 0;
			starttime = new ArrayList<Long>();
		}

		public TimeItem(String header) {
			this();
			this.header = header;
		}

		public TimeItem(TimeItem timeItem) {
			this.totaltime = timeItem.totaltime;
			this.starttime = new ArrayList<Long>(timeItem.starttime);
			this.header = timeItem.header;
		}

		public synchronized void reset() {
			totaltime = 0;
			starttime = new ArrayList<Long>();
		}

		public synchronized void start() {
			starttime.add(System.currentTimeMillis());
		}

		public synchronized void stop() {
			if (starttime.size() == 0) {
				System.err.println("Time counter is not started.");
				return;
			}
			totaltime += (System.currentTimeMillis() - starttime.remove(starttime.size() - 1));
		}

		public synchronized void addTime(long time) {
			totaltime += time;
		}

		public synchronized void setTime(long time) {
			totaltime = time;
		}

		public synchronized long getTime() {
			return totaltime;
		}
	}

	private List<TimeItem> timeItemList;

	public TimeCounter(int size) {
		timeItemList = new ArrayList<TimeItem>();
		for (int i = 0; i < size; i++)
			timeItemList.add(new TimeItem());
	}

	public TimeCounter(int size, String... header) {
		timeItemList = new ArrayList<TimeItem>();
		if (header.length == size)
			for (int i = 0; i < size; i++)
				timeItemList.add(new TimeItem(header[i]));
		else {
			if (header.length > 0)
				System.err.println("Wrong size of headers! Proceed with no headers.");
			for (int i = 0; i < size; i++)
				timeItemList.add(new TimeItem());
		}
	}

	public TimeCounter(TimeCounter tc) {
		timeItemList = new ArrayList<TimeItem>();
		for (int i = 0; i < tc.size(); i++)
			timeItemList.add(new TimeItem(tc.timeItemList.get(i)));
	}

	public void reset(int i) {
		timeItemList.get(i).reset();
	}

	public void resetAll() {
		for (TimeItem timeItem : timeItemList)
			timeItem.reset();
	}

	public void add(int i, long time) {
		timeItemList.get(i).addTime(time);
	}

	public void set(int i, long time) {
		timeItemList.get(i).setTime(time);
	}

	public void start(int i) {
		timeItemList.get(i).start();
	}

	public void end(int i) {
		timeItemList.get(i).stop();
	}

	public void stop(int i) {
		timeItemList.get(i).stop();
	}

	public long get(int i) {
		return timeItemList.get(i).getTime();
	}

	public int size() {
		return timeItemList.size();
	}

	public void outputtime() {
		System.out.println("Time Counter");
		for (int i = 0; i < timeItemList.size(); i++) {
			String description = timeItemList.get(i).header;
			if (description == null)
				description = Integer.toString(i);
			System.out.print(description + ": ");
			System.out.println((double) timeItemList.get(i).getTime() / 1000);
		}
	}

	public static TimeCounter mergeTimeCounter(List<TimeCounter> tclist) {
		if (tclist.isEmpty())
			return null;
		int size = tclist.get(0).size();
		TimeCounter mergetc = new TimeCounter(tclist.get(0));
		mergetc.resetAll();
		for (TimeCounter tc : tclist)
			if (tc.size() != size)
				return null;
			else
				for (int i = 0; i < size; i++)
					mergetc.add(i, tc.get(i));
		return mergetc;
	}
}
