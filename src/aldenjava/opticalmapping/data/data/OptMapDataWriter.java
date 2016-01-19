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

import java.io.IOException;
import java.util.LinkedHashMap;

import joptsimple.OptionSet;

import org.apache.commons.lang.StringUtils;

import aldenjava.opticalmapping.data.DataFormat;
import aldenjava.opticalmapping.data.OMWriter;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;

/**
 * Data Writer
 * 
 * @author Alden
 *
 */
public class OptMapDataWriter extends OMWriter<DataNode> {

	private DataFormat dformat;
	private boolean writtenSingleRefFile;

	public OptMapDataWriter(String filename, int format) throws IOException {
		this(filename, DataFormat.lookup(filename, format));
	}

	public OptMapDataWriter(String filename) throws IOException {
		this(filename, -1);
	}

	public OptMapDataWriter(String filename, DataFormat dformat) throws IOException {
		super(filename);
		if (dformat == null)
			throw new IOException("Unknown format.");
		this.dformat = dformat;
		// Initialize header only after dFormat is initialized
		initializeHeader();
		writtenSingleRefFile = false;

	}

	public OptMapDataWriter(OptionSet options) throws IOException {
		this((String) options.valueOf("optmapout"), (int) options.valueOf("optmapoutformat"));
	}

	@Override
	public void initializeHeader() throws IOException {
		if (dformat != null)
			switch (dformat) {
				case REF:
				case SILICO:
					break;
				case FA01:
					break;
				case SPOTS:
					break;
				case OPT:
					break;
				case SDATA:
					bw.write("#Fragment ID\tReference\tStrand\tStart\tEnd\tSize\tTotalSegments\tSegmentDetail\n");
					break;
				case DATA:
					bw.write("#Fragment ID\tSize\tTotalSegments\tSegmentDetail...\n");
					break;
				case BNX:
					bw.write("# BNX File Version:	0.1\n");
					bw.write("# Label Channels:	1\n");
					bw.write("# Nickase Recognition Site 1:\n");
					bw.write("# Quality Score QX01:	SNR\n");
					bw.write("# Quality Score QX02:	Ave Intensity\n");
					bw.write("# All above comments are fake\n");
					bw.write("#0h	LabelChannel	MapID	Length\n");
					bw.write("#0f	int	int	float\n");
					bw.write("#1h	LabelChannel	LabelPositions[N]\n");
					bw.write("#1f	int	float\n");
					bw.write("#2h	LabelChannel	LabelPositions[N]\n");
					bw.write("#2h	int	float\n");
					bw.write("#Qh	QualityScoreID	QualityScores[N]\n");
					bw.write("#Qf	str	float\n");
					break;
				case CMAP:
					bw.write("# \n");
					bw.write("# CompileDir= \n");
					bw.write("# CMAP File Version:	0.1\n");
					bw.write("# Label Channels:	1\n");
					bw.write("# Nickase Recognition Site 1:	unknown\n");
					bw.write("# Number of Consensus Nanomaps: 0\n");
					bw.write("#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence\n");
					bw.write("#f int	float	int	int	int	float	float	int	int\n");
					break;
				default:
					assert false : "dformat unfound";
			}
	}

	@Override
	public void write(DataNode data) throws IOException {
		if (data != null)
			switch (dformat) {
				case REF:
				case SILICO:
					bw.write(data.name + "\t" + Long.toString(data.size) + "\t" + Integer.toString(data.refp.length) + "\n");
					for (int i = 0; i < data.refp.length - 1; i++)
						bw.write(Long.toString(data.refp[i]) + "\t");
					if (data.refp.length > 0)
						bw.write(Long.toString(data.refp[data.refp.length - 1]) + "\n");
					else
						bw.write("\n");
					break;
				case FA01:
					bw.write(">");
					bw.write(data.name);
					bw.write("\n");
					for (int i = 0; i < data.refp.length + 1; i++) {
						bw.write(StringUtils.repeat("0", (int) data.getRefl(i)));
						if (i != data.refp.length + 1 - 1)
							bw.write("1");
					}

					bw.write("\n");
					break;
				case SPOTS:
					if (writtenSingleRefFile) {
						System.out.println("No more than one reference can be written on the same spot file.");
						return;
					}
					writtenSingleRefFile = true;
					bw.write("#\n");
					bw.write("#	Reference Size (Bp):	" + Long.toString(data.size) + "\n");
					bw.write("#	N Colors:	1\n");
					bw.write("NickID Color Location\n");

					long refp[] = data.refp;
					for (int i = 0; i < refp.length; i++) {
						bw.write(Integer.toString(i));
						bw.write(" 1 ");
						bw.write(Long.toString(refp[i]));
						bw.write("\n");
					}
					break;
				case OPT:
					if (writtenSingleRefFile) {
						System.out.println("No more than one reference can be written on the same opt file.");
						return;
					}
					writtenSingleRefFile = true;

					for (int i = 0; i < data.refp.length + 1; i++) {
						long r = data.getRefl(i);
						bw.write(String.format("%f\t%f\n", r / (double) 1000 + 0.001, (r / (double) 1000) * 0.05 + 0.001));
					}
					break;
				case DATA:
					StringBuilder s2 = new StringBuilder();
					for (int i = 0; i < data.getTotalSegment(); i++) {
						s2.append(Long.toString(data.getRefl(i)));
						if (i != data.getTotalSegment() - 1)
							s2.append(";");
					}

					bw.write(String.format("%s\t%d\t%d\t%s\n", data.name, data.length(), data.getTotalSegment(), s2.toString()));
					break;
				case SDATA:
					String r = "";
					long start = -1;
					long stop = -1;
					int simuStrand = 0;

					if (data.hasSimulationInfo()) {
						r = data.simuInfo.simuRegion.ref;
						start = data.simuInfo.simuRegion.start;
						stop = data.simuInfo.simuRegion.stop;
						simuStrand = data.simuInfo.simuStrand;
					}
					StringBuilder s = new StringBuilder();
					for (int i = 0; i < data.getTotalSegment(); i++) {
						s.append(Long.toString(data.getRefl(i)));
						if (i != data.getTotalSegment() - 1)
							s.append(";");
					}
					bw.write(String.format("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n", data.name, r, simuStrand == 1 ? "+" : simuStrand == -1 ? "-" : "", start, stop, data.length(), data.getTotalSegment(),
							s.toString()));
					break;
				case BNX:
					bw.write("0\t");
					bw.write(data.name);
					bw.write("\t");
					bw.write(Long.toString(data.length()));
					bw.write(".0");
					bw.write("\n");
					bw.write("1");
					for (int i = 0; i < data.refp.length; i++) {
						bw.write("\t");
						bw.write(Long.toString(data.refp[i]));
						bw.write(".0");
					}
					bw.write("\n");
					bw.write("QX01"); // Default as 10.0
					if (data instanceof BnxDataNode)
						for (int i = 0; i < data.refp.length; i++)
							bw.write("\t" + Double.toString(((BnxDataNode) data).snr[i]));
					else
						for (int i = 0; i < data.refp.length; i++)
							bw.write("\t10.0");
					bw.write("\n");
					bw.write("QX02"); // Default as 0.05
					if (data instanceof BnxDataNode)
						for (int i = 0; i < data.refp.length; i++)
							bw.write("\t" + Double.toString(((BnxDataNode) data).intensity[i]));
					else
						for (int i = 0; i < data.refp.length; i++)
							bw.write("\t0.05");
					bw.write("\n");
					break;
				case CMAP:

					String id = data.name;
					long size = data.size;
					int totalsites = data.getTotalSegment() - 1;
					// l[3] //site x
					int labelchannel = 1;
					double stddev = 0.0;
					int coverage = 1;
					int occurence = 1;

					for (int i = 0; i < data.refp.length + 1; i++) {
						long pos;
						if (i == data.refp.length) {
							pos = size;
							labelchannel = 0;
						} else
							pos = data.refp[i];
						bw.write(String.format("%s\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%d\n", id, size, totalsites, i + 1, labelchannel, pos, stddev, coverage, occurence));
					}
					break;
				default:
					assert false : "dformat unfound";
			}
	}

	public static void writeAll(String filename, int fileformat, LinkedHashMap<String, DataNode> fragmentmap) throws IOException {
		OptMapDataWriter omdw = new OptMapDataWriter(filename, fileformat);
		omdw.writeAll(fragmentmap);
		omdw.close();
	}

	public static void writeAll(OptionSet options, LinkedHashMap<String, DataNode> fragmentmap) throws IOException {
		OptMapDataWriter omdw = new OptMapDataWriter(options);
		omdw.writeAll(fragmentmap);
		omdw.close();
	}

	public static void assignOptions(ExtendOptionParser parser) {
		parser.addHeader("Data Writer Options", 1);
		parser.accepts("optmapout", "output file").withOptionalArg().ofType(String.class).defaultsTo("");
		parser.accepts("optmapoutformat", DataFormat.getFormatHelp()).withOptionalArg().ofType(Integer.class).defaultsTo(-1);
	}

}
