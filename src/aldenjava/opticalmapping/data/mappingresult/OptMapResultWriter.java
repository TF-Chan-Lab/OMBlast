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

import java.io.IOException;
import java.util.List;

import joptsimple.OptionSet;
import aldenjava.opticalmapping.data.OMWriter;
import aldenjava.opticalmapping.data.data.DataNode;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;

/**
 * Alignment Result Writer
 * 
 * @author Alden
 *
 */
public class OptMapResultWriter extends OMWriter<OptMapResultNode> {

	private ResultFormat rformat;
	private boolean writeunmap;
	private boolean multiple;
	private boolean writeinfo;
	private int xmapDummyID;

	public OptMapResultWriter(OptionSet options) throws IOException {
		this((String) options.valueOf("optresout"), (int) options.valueOf("optresoutformat"), (boolean) options.valueOf("writeunmap"), (boolean) options.valueOf("multiple"), (boolean) options
				.valueOf("writeinfo"));
	}

	public OptMapResultWriter(String filename) throws IOException {
		this(filename, -1);
	}

	public OptMapResultWriter(String filename, int format) throws IOException {
		this(filename, format, true, true, true);
	}

	public OptMapResultWriter(String filename, int format, boolean writeunmap, boolean multiple, boolean writeinfo) throws IOException {
		this(filename, ResultFormat.lookup(filename, format), writeunmap, multiple, writeinfo);
	}

	public OptMapResultWriter(String filename, ResultFormat rformat) throws IOException {
		this(filename, rformat, true, true, true);
	}

	public OptMapResultWriter(String filename, ResultFormat rformat, boolean writeunmap, boolean multiple, boolean writeinfo) throws IOException {
		super(filename, false);
		this.rformat = rformat;
		this.initializeHeader();
		this.writeunmap = writeunmap;
		this.multiple = multiple;
		this.writeinfo = writeinfo;
		if (rformat == ResultFormat.XMAP && this.writeunmap) {
			System.out.println("XMAP does not support unmap entry. Forced writeunmap to false");
			this.writeunmap = false;
		}
		if (rformat == ResultFormat.XMAP && this.writeunmap) {
			System.out.println("XMAP does not support molecule information. Forced writeinfo to false");
			this.writeinfo = false;
		}
	}

	@Override
	protected void initializeHeader() throws IOException {
		super.initializeHeader();
		if (rformat != null)
			switch (rformat) {
				case OMA:
					bw.write("#OMA File format version v1.1\n");
					bw.write("#MoleID\tNoOfSeg\tMoleInfo\tRName\tOrient\tScore\tConfidence\tRegSegStart\tRefSegStop\tMoleSegStart\tMoleSegStop\tRefStartCoord\tRefStopCoord\tCigar\n");
					break;
				case OMD:
					bw.write("#MoleID\tfromRef\tgenomeStrand\tgenomestart\tgenomestop\tFragmentSize\tFragmentLabels\tFragment\tRefName\tStartPos\tEndPos\tStrand\tSubRefStart\tSubRefStop\tSubFragmentStart\tSubFragmentStop\tSubFragmentRatio\tScore\tCIGAR\tConfidence\tFalseP\tFalseN\tScale\tFPRate\tFNRate\tMapped\n");
					break;
				case XMAP:
					bw.write("# XMAP File Version:\t0.1\n");
					bw.write("#h XmapEntryID\tQryContigID\tRefcontigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\n");
					bw.write("#f int\tstring\tstring\tfloat\tfloat\tfloat\tfloat\tstring\tfloat\tstring\n");
					xmapDummyID = 0;
					break;
				case VAL:
				case SOMA:
				case PSL:
					System.err.println("The result format " + rformat.toString() + " is not supported for output");
					break;
				default:
					assert false : "rformat unfound";
			}
	}

	@Override
	public void write(OptMapResultNode result) throws IOException {
		if (result != null && (result.isUsed() || writeunmap))
			switch (rformat) {
				case OMA:
					writeOMA(result);
					break;
				case OMD:
					writeOMD(result);
					break;
				case XMAP:
					writeXMAP(result);
					break;
				case VAL:
				case SOMA:
				case PSL:
					System.err.println("The selected result format is not supported for output");
					break;
				default:
					assert false : "rformat unfound";
			}
	}

	public void write(List<OptMapResultNode> resultlist) throws IOException {
		if (multiple)
			for (OptMapResultNode result : resultlist)
				write(result);
		else
			write(resultlist.get(0));
	}

	private void writeOMA(OptMapResultNode result) throws IOException {
		DataNode f = result.parentFrag;
		String finfo = "";
		if (writeinfo)
			finfo = f.getReflString();
		// OMA1.1
		if (result.isUsed())
			bw.write(String.format("%s\t%d\t%s\t%s\t%s\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", f.name, f.getTotalSegment(), finfo, result.mappedRegion.ref, result.mappedstrand == 1 ? "+"
					: result.mappedstrand == -1 ? "-" : "", result.mappedscore, result.confidence, result.subrefstart, result.subrefstop, result.subfragstart, result.subfragstop,
					result.mappedRegion.start, result.mappedRegion.stop, result.cigar));
		else
			bw.write(String.format("%s\t%d\t%s\t%s\n", f.name, f.getTotalSegment(), finfo, "Unmapped"));
		// OMA1.0
		// if (result.isUsed())
		// bw.write(String.format("%s\t%d\t%s\t%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", f.id, f.fragment.length, finfo, result.mappedRegion.ref, result.mappedstrand==1?"forward":result.mappedstrand==-1?"reverse":"", result.mappedscore, result.subrefstart, result.subrefstop, result.subfragstart, result.subfragstop, result.mappedRegion.start, result.mappedRegion.stop, result.cigar));
		// else
		// bw.write(String.format("%s\t%d\t%s\t%s\n", f.id, f.fragment.length, finfo, "Unmapped"));
	}

	private void writeOMD(OptMapResultNode result) throws IOException {

		DataNode f = result.parentFrag;
		String finfo = "";
		if (writeinfo)
			finfo = f.getReflString();
		if (result.isUsed())
			if (f.hasSimulationInfo())
				bw.write(String.format("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%f\t%d\t%d\t%f\t%f\t%f\t%b\n", f.name, f.simuInfo.simuRegion.ref,
						f.simuInfo.simuStrand == 1 ? "+" : f.simuInfo.simuStrand == -1 ? "-" : "", f.simuInfo.simuRegion.start, f.simuInfo.simuRegion.stop, result.length(), f.getTotalSegment(),
						finfo, result.mappedRegion.ref, result.mappedRegion.start, result.mappedRegion.stop, result.mappedstrand == 1 ? "+" : result.mappedstrand == -1 ? "-" : "", result.subrefstart,
						result.subrefstop, result.subfragstart, result.subfragstop, result.getSubFragRatio(), result.mappedscore, result.cigar, result.confidence, result.getFP(), result.getFN(),
						result.getMapScale(), result.getFPRate(), result.getFNRate(), result.correctlyMapped()));
			else
				bw.write(String.format("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%f\t%d\t%d\t%f\t%f\t%f\t%b\n", f.name, "", "", -1, -1, result.length(),
						f.getTotalSegment(), finfo, result.mappedRegion.ref, result.mappedRegion.start, result.mappedRegion.stop,
						result.mappedstrand == 1 ? "+" : result.mappedstrand == -1 ? "-" : "", result.subrefstart, result.subrefstop, result.subfragstart, result.subfragstop,
						result.getSubFragRatio(), result.mappedscore, result.cigar, result.confidence, result.getFP(), result.getFN(), result.getMapScale(), result.getFPRate(), result.getFNRate(),
						null));
		else if (f.hasSimulationInfo())
			bw.write(String.format("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n", f.name, f.simuInfo.simuRegion.ref, f.simuInfo.simuStrand == 1 ? "+" : f.simuInfo.simuStrand == -1 ? "-" : "",
					f.simuInfo.simuRegion.start, f.simuInfo.simuRegion.stop, result.length(), f.getTotalSegment(), finfo, "Unmapped"));
		else
			bw.write(String.format("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n", f.name, "", "", -1, -1, result.length(), f.getTotalSegment(), finfo, "Unmapped"));
	}

	private void writeXMAP(OptMapResultNode result) throws IOException {
		xmapDummyID++;
		DataNode f = result.parentFrag;
		bw.write(String.format("%s\t%s\t%s\t%d.0\t%d.0\t%d.0\t%d.0\t%s\t%.2f\t%s\n", Integer.toString(xmapDummyID), f.name, result.mappedRegion.ref, result.parentFrag.refp[result.subfragstart
				- (result.mappedstrand == 1 ? 1 : 0)], result.parentFrag.refp[result.subfragstop - (result.mappedstrand == 1 ? 0 : 1)], result.mappedRegion.start, result.mappedRegion.stop,
				result.mappedstrand == 1 ? "+" : result.mappedstrand == -1 ? "-" : "", result.mappedscore, result.cigar));
	}

	public static void assignOptions(ExtendOptionParser parser) {
		parser.addHeader("Result Writer Options", 1);
		parser.accepts("optresout", "Result file name").withRequiredArg().ofType(String.class);
		parser.accepts("optresoutformat", "Result file format " + ResultFormat.getFormatHelp()).withOptionalArg().ofType(Integer.class).defaultsTo(-1);
		parser.accepts("writeunmap", "Write discarded or unmapped molecules.").withOptionalArg().ofType(Boolean.class).defaultsTo(true);
		parser.accepts("multiple", "Write multiple maps for a molecule.").withOptionalArg().ofType(Boolean.class).defaultsTo(true);
		parser.accepts("writeinfo", "Write information of a molecule.").withOptionalArg().ofType(Boolean.class).defaultsTo(true);
	}
}
