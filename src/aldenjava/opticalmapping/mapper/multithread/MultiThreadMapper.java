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


package aldenjava.opticalmapping.mapper.multithread;

import java.io.Closeable;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import joptsimple.OptionSet;
import aldenjava.common.TimeCounter;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;
import aldenjava.opticalmapping.mapper.Mapper;
import aldenjava.opticalmapping.mapper.MapperConstructionException;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;

/**
 * <code>MultiThreadMapper</code> provides a standard multi-threading procedures on any <code>Mapper</code>.
 * 
 * @author Alden
 *
 */
public class MultiThreadMapper implements Closeable {

	// Parameters
	private Mapper[] ommapper;
	private Mapper targetmapper;
	private boolean mapperNeedToSetParameters;
	private List<Future<List<OptMapResultNode>>> futureList;
	private int nrOfProcessors = 1;
	private ExecutorService es;
	private CompletionService<List<OptMapResultNode>> ecs;
	private Constructor<? extends Mapper> ctor;

	/**
	 * Constructs a new <code>MultiThreadMapper</code> based on any class extending <code>Mapper</code> using reflection
	 * 
	 * @param mapperClass
	 * @param optrefmap
	 * @throws MapperConstructionException 
	 */
	public MultiThreadMapper(Class<? extends Mapper> mapperClass, LinkedHashMap<String, DataNode> optrefmap) throws MapperConstructionException {
		try {
			this.ctor = mapperClass.getConstructor(LinkedHashMap.class);
			targetmapper = (Mapper) ctor.newInstance(optrefmap);
		} catch (NoSuchMethodException | SecurityException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
			throw new MapperConstructionException(e);
		}
		mapperNeedToSetParameters = true;
	}

	/**
	 * Constructs a new <code>MultiThreadMapper</code> based on an existing <code>Mapper</code>
	 * 
	 * @param mapper
	 */
	public MultiThreadMapper(Mapper mapper) {
		targetmapper = mapper;
		mapperNeedToSetParameters = false;
	}

	/**
	 * Sets parameters for <code>MultiThreadMapper</code>. If the parameters of <code>targetmapper</code> is not set, also set the parameters for the it as well
	 * 
	 * @param options
	 * @throws IOException
	 */
	public void setParameters(OptionSet options) throws IOException {
		if (mapperNeedToSetParameters)
			targetmapper.setParameters(options);
		setParameters((int) options.valueOf("thread"));
	}

	/**
	 * Sets parameters for <code>MultiThreadMapper</code>. Once the <code>nrOfProcessors</code> is set, a new thread pool is initialized. Each thread is corresponding to a copy of the <code>targetmapper</code>.
	 * 
	 * @param thread
	 *            number of copies of <code>ommapper</code> to be used
	 */
	public void setParameters(int thread) {
		this.nrOfProcessors = thread;
		es = Executors.newFixedThreadPool(nrOfProcessors);
		ecs = new ExecutorCompletionService<List<OptMapResultNode>>(es);

		futureList = new ArrayList<Future<List<OptMapResultNode>>>();
		for (int i = 0; i < nrOfProcessors; i++)
			futureList.add(null);
		ommapper = new Mapper[nrOfProcessors];
		ommapper[0] = targetmapper;

		for (int i = 1; i < ommapper.length; i++) {
			ommapper[i] = ommapper[0].copy();
		}
	}

	// Core
	/**
	 * Feeds this instance of <code>MultiThreadMapper</code> a new data to be processed. Returns <code>false</code> If all <code>ommapper</code> are busy and this data is not processed.
	 * 
	 * @param data
	 * @return <code>true</code> if the data is successfully fed to one of the <code>ommapper</code>; <code>false</code> if the data is not processed and need to be fed again.
	 */
	public boolean startNext(DataNode data) {
		for (int i = 0; i < nrOfProcessors; i++)
			if (futureList.get(i) == null) {
				ommapper[i].setData(data);
				futureList.set(i, ecs.submit(ommapper[i]));
				return true;
			}
		return false;
	}

	/**
	 * Returns the status of this instance of <code>MultiThreadMapper</code> according to running status of <code>ommapper</code> and availability of results
	 * 
	 * @return <code>-1</code> if all <code>ommapper</code> are idle and no results are available; <code>0</code> if some of the <code>ommapper</code> are running but no results are available yet; <code>1</code> if some of the <code>ommapper</code> are running and some results are available to be taken
	 */
	public int getStatus() {
		// -1 nothing running
		// 0 running
		// 1 running and have results
		boolean nothingRun = true;
		for (int i = 0; i < nrOfProcessors; i++)
			if (futureList.get(i) != null) {
				nothingRun = false;
				if (futureList.get(i).isDone())
					return 1;
			}
		if (nothingRun)
			return -1;
		else
			return 0;
	}

	/**
	 * Returns the alignment results from next available <code>ommapper</code>. If all <code>ommapper</code> are running, wait until the next result is available.
	 * 
	 * @return alignment result wrapped in <code>MultiThreadResultNode</code>
	 * @throws IllegalStateException
	 *             if this instance of <code>MultiThreadMapper</code> has status of <code>-1</code>, i.e. no results are available
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
	public MultiThreadResultNode getNextResult() throws IllegalStateException, InterruptedException, ExecutionException {
		if (getStatus() == -1) // all tasks finish, nothing to wait
			throw new IllegalStateException("All results are taken");
		Future<List<OptMapResultNode>> future = ecs.take();

		int i = this.futureList.indexOf(future);
		MultiThreadResultNode multinode = new MultiThreadResultNode(ommapper[i].getData(), futureList.get(i).get());
		futureList.set(i, null);
		return multinode;
	}

	// Call to activate
	/**
	 * Perform alignments on all the data in <code>fragmentmap</code>
	 * 
	 * @param fragmentmap
	 * @return alignment results
	 */
	public LinkedHashMap<String, List<OptMapResultNode>> mapAll(LinkedHashMap<String, DataNode> fragmentmap) {
		if (getStatus() != -1) {
			throw new IllegalStateException("Mapper is already running.");
		} else {
			LinkedHashMap<String, List<OptMapResultNode>> fragmentmaplistmap = new LinkedHashMap<String, List<OptMapResultNode>>();
			try {
				for (DataNode fragment : fragmentmap.values()) {
					while (!startNext(fragment)) {
						MultiThreadResultNode multinode = getNextResult();
						fragmentmaplistmap.put(multinode.data.name, multinode.alignmentResults);
					}

				}
				while (getStatus() != -1) {
					// while (getStatus() == 0);
					MultiThreadResultNode multinode = getNextResult();
					fragmentmaplistmap.put(multinode.data.name, multinode.alignmentResults);
				}
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
				// Unknown reason for interruption, but should continue to handle the result.
			}
			return fragmentmaplistmap;
		}
	}

	/**
	 * Returns the merged alignment time of all <code>ommapper</code>
	 * 
	 * @return merged alignment time
	 */
	public TimeCounter getMappingTime() {
		List<TimeCounter> tclist = new ArrayList<TimeCounter>();
		for (Mapper mapper : ommapper)
			tclist.add(mapper.tc);
		return TimeCounter.mergeTimeCounter(tclist);
	}

	@Override
	public void close() {
		es.shutdown();
	}

	public static void assignOptions(ExtendOptionParser parser, int level) {
		parser.addHeader("Multi-thread Options", level);
		parser.accepts("thread", "Number of threads").withOptionalArg().ofType(Integer.class).defaultsTo(1);
	}
}
