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

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

import joptsimple.OptionSet;
import aldenjava.opticalmapping.GenomicPosNode;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;
import aldenjava.opticalmapping.mapper.Mapper;
import aldenjava.opticalmapping.mapper.MapperConstructionException;
import aldenjava.opticalmapping.mapper.seeding.SeedDatabase;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;

/**
 * The OMBlastMapper class
 * 
 * @author Alden
 *
 */
public class OMBlastMapper extends Mapper {

	private OMBlastCore blastcore = null;

	/**
	 * Constructs an <code>OMBlastMapper</code> based on the reference information
	 * 
	 * @param optrefmap
	 *            the reference information
	 */
	public OMBlastMapper(LinkedHashMap<String, DataNode> optrefmap) {
		super(optrefmap);
	}

	@Override
	public void setParameters(OptionSet options) throws IOException {
		super.setParameters(options);
		this.setParameters((int) options.valueOf("seedingmode"), (boolean) options.valueOf("local"), (int) options.valueOf("falselimit"), (int) options.valueOf("k"),
				(int) options.valueOf("maxnosignal"), (int) options.valueOf("meas"), (double) options.valueOf("ear"), (int) options.valueOf("match"), (int) options.valueOf("fpp"),
				(int) options.valueOf("fnp"), (int) options.valueOf("maxseedno"));
	}

	public void setParameters(int seedingmode, boolean allowLocalAlignment, int falselimit, int kmerlen, int maxnosignalregion, int measure, double ear, int matchscore, int falseppenalty,
			int falsenpenalty, int maxSeedNumber) {
		blastcore = new OMBlastCore(optrefmap);
		blastcore.setParameters(seedingmode, kmerlen, maxnosignalregion, allowLocalAlignment, measure, ear, matchscore, falseppenalty, falsenpenalty, falselimit, maxSeedNumber);
	}

	@Override
	public List<OptMapResultNode> getResult(DataNode data, List<GenomicPosNode> regionList) {
		if (data == null)
			return null;
		else {
			blastcore.restrictRegion(regionList);
			return blastcore.getResult(data);
		}
	}

	@Override
	public OMBlastMapper copy() {
		OMBlastMapper mapper = new OMBlastMapper(optrefmap);
		mapper.blastcore = this.blastcore.copy();
		super.setCopyMapperParameters(mapper);
		return mapper;
	}

	public static void assignOptions(ExtendOptionParser parser, int level) {

		Mapper.assignOptions(parser, level);

		parser.addHeader("OMBlastMapper Options", level);
		parser.accepts("local", "Enable Local Alignment.").withOptionalArg().ofType(Boolean.class).defaultsTo(true);
		parser.accepts("meas", "Measurement error.").withOptionalArg().ofType(Integer.class).defaultsTo(500);
		parser.accepts("ear", "Error acceptable range.").withOptionalArg().ofType(Double.class).defaultsTo(0.1);
		parser.accepts("match", "Match score.").withOptionalArg().ofType(Integer.class).defaultsTo(5);
		parser.accepts("fpp", "False positive penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
		parser.accepts("fnp", "False negative penalty").withOptionalArg().ofType(Integer.class).defaultsTo(2);
		parser.accepts("falselimit", "Max consecutive false signals").withOptionalArg().ofType(Integer.class).defaultsTo(5);
		parser.accepts("maxseedno", "Maximum similar seed number on query").withOptionalArg().ofType(Integer.class).defaultsTo(10);
		SeedDatabase.assignOptions(parser, level + 1);
	}

	public static void main(String[] args) throws IOException, MapperConstructionException {
		Mapper.standardMapperProcedure(args, OMBlastMapper.class);
	}
}