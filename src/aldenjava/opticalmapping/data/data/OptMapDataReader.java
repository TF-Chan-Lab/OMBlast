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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import javax.xml.stream.FactoryConfigurationError;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import joptsimple.OptionSet;

import org.apache.commons.lang.ArrayUtils;

import aldenjava.opticalmapping.GenomicPosNode;
import aldenjava.opticalmapping.data.DataFormat;
import aldenjava.opticalmapping.data.OMReader;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;

/**
 * The optical mapping data reader. It supports various input format
 * 
 * @author Alden
 *
 * @see DataFormat
 */
public class OptMapDataReader extends OMReader<DataNode> {

	/**
	 * The data format of the input file
	 */
	private DataFormat dformat;
	
	private boolean singleRecord;
	private String filename;
	private double bnxSNR = 3.0;
	private XMLStreamReader xmlReader;

	public OptMapDataReader(OptionSet options) throws IOException {
		this((String) options.valueOf("optmapin"), (int) options.valueOf("optmapinformat"), (double) options.valueOf("bnxsnr"));
	}

	public OptMapDataReader(String filename) throws IOException {
		this(filename, -1);
	}

	public OptMapDataReader(String filename, int format) throws IOException {
		this(filename, DataFormat.lookup(filename, format));
	}

	public OptMapDataReader(String filename, int format, double bnxSNR) throws IOException {
		this(filename, DataFormat.lookup(filename, format));
		this.bnxSNR = bnxSNR;
	}

	public OptMapDataReader(String filename, DataFormat dformat) throws IOException {
		super(filename);
		if (dformat == null)
			throw new NullPointerException("dformat");
		this.dformat = dformat;
		this.filename = filename;
		if (dformat == DataFormat.XML)
			createXMLReader();
	}

	public OptMapDataReader(InputStream stream, DataFormat dformat) throws IOException {
		super(stream);
		if (dformat == null)
			throw new NullPointerException("dformat");
		this.dformat = dformat;
		this.filename = "";
		if (dformat == DataFormat.XML)
			createXMLReader();

	}

	@Override
	public DataNode read() throws IOException {
		if (nextline == null)
			return null;
		else {
			switch (dformat) {
				case REF:
				case SILICO:
					return parseREF();
				case FA01:
					return parseFA01();
				case SPOTS:
					return parseSPOT();
				case OPT:
					return parseOPT();
				case DATA:
					return parseDATA();
				case SDATA:
					return parseSDATA();
				case BNX:
					return parseBNX();
				case CMAP:
					return parseCMAP();
				case XML:
					return parseXML();
				default:
					return null;

			}
		}
	}

	private DataNode parseREF() throws IOException {
		if (nextline == null)
			return null;
		String[] l = nextline.split("\t");
		String name = l[0];
		long size = Long.parseLong(l[1]);
		long refplen = Integer.parseInt(l[2]);
		proceedNextLine();
		if (refplen == 0)
			return new DataNode(name, size);
		if (nextline != null) {
			List<Long> refp = new ArrayList<Long>();

			String[] line = nextline.trim().split("\t");
			for (int i = 0; i < line.length; i++)
				refp.add(Long.parseLong(line[i]));
			if (refplen != refp.size())
				System.err.println("Warning: " + name + "\nNumber of labels stated (" + refplen + ") does not match true number of labels (" + refp.size() + ").");
			proceedNextLine();

			return new DataNode(name, size, ArrayUtils.toPrimitive(refp.toArray(new Long[refp.size()])));
		} else {
			System.err.println("Warning: Incomplete record found: " + name);
			return null;
		}
	}

	private DataNode parseFA01() throws IOException {
		if (nextline == null)
			return null;
		String name = nextline.substring(1);
		StringBuilder zeroOneString = new StringBuilder();
		String s;
		while ((s = br.readLine()) != null) {
			if (s.startsWith(">")) {
				nextline = s;
				break;
			} else
				zeroOneString = zeroOneString.append(s);
		}

		if (zeroOneString.length() == 0)
			return null;
		else {
			long ref_size = zeroOneString.length();
			long zero = 0;
			List<Long> refl = new ArrayList<Long>();
			List<Long> refp = new ArrayList<Long>();
			for (long i = 0; i < ref_size; i++) {
				if (zeroOneString.charAt((int) i) == '0')
					zero++;
				else if (zeroOneString.charAt((int) i) == '1') {
					refl.add(zero);
					refp.add(i);
					zero = 0;
				}
			}
			refl.add(zero);
			return new DataNode(name, zeroOneString.length(), ArrayUtils.toPrimitive(refp.toArray(new Long[refp.size()])));
		}
	}

	private DataNode parseSPOT() throws IOException {
		if (nextline == null)
			return null;
		if (singleRecord) {
			System.err.println("There is only one reference in spot file.");
			return null;
		}
		singleRecord = true;
		String s;
		long refsize = 0;
		br.reset();
		while ((s = br.readLine()).startsWith("#")) {
			if (s.contains("Reference Size")) {
				String[] l = s.split("\t");
				refsize = Long.parseLong(l[l.length - 1]);
			}
		}
		while (s.startsWith("NickID"))
			s = br.readLine();
		List<Long> refl = new ArrayList<Long>();
		List<Long> refp = new ArrayList<Long>();

		while (s != null) {
			String[] line = s.split("\\s+");
			refp.add(Long.parseLong(line[line.length - 1]));
			s = br.readLine();
		}
		br.close();
		refl.add(refp.get(0));
		for (int i = 0; i < refp.size() - 1; i++)
			refl.add(refp.get(i + 1) - refp.get(i));
		refl.add(refsize - refp.get(refp.size() - 1));

		return new DataNode(new File(filename).getName(), refsize, ArrayUtils.toPrimitive(refp.toArray(new Long[refp.size()])));
	}

	private DataNode parseOPT() throws IOException {
		if (nextline == null)
			return null;
		if (singleRecord) {
			System.err.println("There is only one reference in opt file.");
			return null;
		}
		singleRecord = true;

		List<Long> refl = new ArrayList<Long>();
		String s = br.readLine();
		while (s != null) {
			String[] line = s.split("\\s+");
			refl.add((long) (Double.parseDouble(line[0]) * 1000));
			s = br.readLine();
		}
		return new DataNode(new File(filename).getName(), ArrayUtils.toPrimitive(refl.toArray(new Long[refl.size()])));
	}

	private DataNode parseDATA() throws IOException {
		String s = this.nextline;
		s = s.trim();
		String[] l = s.split("\t");
		String name = l[0];
		long size = Long.parseLong(l[1]);
		int totalSegment = Integer.parseInt(l[2]);
		long[] refl = DataNode.parseReflInString(l[3], ";");

		// Validation
		if (refl.length != totalSegment) {
			System.err.println("Warning: Inconsistent total segments and fragment details");
		}
		long mysize = 0;
		for (int i = 0; i < refl.length - 1; i++) {
			mysize += refl[i] + 1; // need +1 here?
		}
		mysize += refl[refl.length - 1];
		if (size != mysize) {
			System.err.println("Warning: Inconsistent size and fragment details");
		}
		proceedNextLine();
		return new DataNode(name, refl);
	}

	private DataNode parseSDATA() throws IOException {
		String s = this.nextline;
		s = s.trim();
		String[] l = s.split("\t");
		String name = l[0];
		String ref = l[1];
		int simuStrand = (l[2].equalsIgnoreCase("forward") || l[2].equalsIgnoreCase("+")) ? 1 : (l[2].equalsIgnoreCase("reverse") || l[2].equalsIgnoreCase("-")) ? -1 : 0;
		long genomestart = Long.parseLong(l[3]);
		long genomestop = Long.parseLong(l[4]);
		long size = Long.parseLong(l[5]);
		int totalSegment = Integer.parseInt(l[6]);
		long[] refl = DataNode.parseReflInString(l[7], ";");

		// Validation
		if (refl.length != totalSegment) {
			System.err.println("Warning: Inconsistent total segments and fragment details");
		}
		long mysize = 0;
		for (int i = 0; i < refl.length - 1; i++) {
			mysize += refl[i] + 1; // need +1 here?
		}
		mysize += refl[refl.length - 1];
		if (size != mysize) {
			System.err.println("Warning: Inconsistent size and fragment details");
		}
		proceedNextLine();
		DataNode data = new DataNode(name, refl);
		if (SimulationInfo.checkInfoValid(ref, genomestart, genomestop, simuStrand))
			data.importSimulationInfo(new GenomicPosNode(ref, genomestart, genomestop), simuStrand);
		return data;

	}

	private DataNode parseBNX() throws IOException {

		// Currently unused information, Only support one (single-color) channel

		// int labelChannel;
		// double avgIntensity;
		// double snr;
		// int noOfLables;
		// String originalMoleculeID;
		// int scanNumber;
		// int scanDirection;
		// String chipID;
		// int flowCell;
		// int runID;
		// int globalScanNumber;

		// BNX File Version 1.0
		// #0h LabelChannel MoleculeId Length AvgIntensity SNR NumberofLabels OriginalMoleculeId ScanNumber ScanDirection ChipId Flowcell
		// #0f int int float float float int int int int string int
		// #1h LabelChannel LabelPositions[N]
		// #1f int float
		// #2h LabelChannel LabelPositions[N]
		// #2h int float
		// #Qh QualityScoreID QualityScores[N]
		// #Qf str float
		// BNX File Version 1.2
		// #0h LabelChannel MoleculeId Length AvgIntensity SNR NumberofLabels OriginalMoleculeId ScanNumber ScanDirection ChipId Flowcell RunId GlobalScanNumber
		// #0f int int float float float int int int int string int int int
		// #1h LabelChannel LabelPositions[N]
		// #1f int float
		// #2h LabelChannel LabelPositions[N]
		// #2h int float
		// #Qh QualityScoreID QualityScores[N]
		// #Qf str float
		String name = "";
		long size = -1;

		double[] snr = null;
		double[] intensity = null;

		long[] refp = null;
		boolean gotNameSizeInfo = false;
		boolean gotDetailInfo = false;
		boolean gotSNRInfo = false;
		boolean gotIntensityInfo = false;
		do {

			String s = this.nextline.trim();

			String[] l = s.split("\t");
			switch (l[0]) {
				case "0":
					name = l[1];
					size = (long) Double.parseDouble(l[2]);
					gotNameSizeInfo = true;
					break;
				case "1":
					refp = new long[l.length - 2]; // last element should be the size of molecule
					for (int i = 1; i < l.length - 1; i++)
						refp[i - 1] = (long) Double.parseDouble(l[i]);
					gotDetailInfo = true;
					break;
				case "QX01":
				case "QX11":
					snr = new double[l.length - 1];
					for (int i = 1; i < l.length; i++)
						snr[i - 1] = Double.parseDouble(l[i]);
					gotSNRInfo = true;
					break;
				case "QX02":
				case "QX12":
					intensity = new double[l.length - 1];
					for (int i = 1; i < l.length; i++)
						intensity[i - 1] = Double.parseDouble(l[i]);
					gotIntensityInfo = true;
					break;
				default:

			}
			proceedNextLine();
		} while ((this.nextline != null)
				&& (this.nextline.startsWith("1") || this.nextline.startsWith("QX01") || this.nextline.startsWith("QX02") || this.nextline.startsWith("QX11") || this.nextline.startsWith("QX12")));
		if (gotNameSizeInfo && gotDetailInfo && gotSNRInfo && gotIntensityInfo) {
			List<Long> refpList = new ArrayList<Long>();
			for (int i = 0; i < refp.length; i++)
				if (snr[i] >= bnxSNR)
					refpList.add(refp[i]);
			refp = ArrayUtils.toPrimitive(refpList.toArray(new Long[refpList.size()]));
			BnxDataNode data = new BnxDataNode(name, size, refp, intensity, snr);
			return data;
		} else {
			System.err.println("Warning: incomplete record is found.");
			return null;
		}
	}

	private DataNode parseCMAP() throws IOException {
		String s = this.nextline;
		s = s.trim();
		if (s.isEmpty())
			return null;
		String[] l = s.split("\\s+");
		String name = l[0];
		long size = (long) (Double.parseDouble(l[1]));
		int totalSignal = Integer.parseInt(l[2]);
		// l[3] //site x

		long[] refp = new long[totalSignal];
		for (int i = 0; i < totalSignal + 1; i++) {
			s = this.nextline;
			s = s.trim();
			l = s.split("\\s+");
			int labelchannel = Integer.parseInt(l[4]); // not used till supporting double color
			if (labelchannel != 0)
				refp[i] = (long) Double.parseDouble(l[5]);
			proceedNextLine();
		}

		return new DataNode(name, size, refp);
	}

	private void createXMLReader() throws IOException {
		try {
			xmlReader = XMLInputFactory.newInstance().createXMLStreamReader(super.br);
		} catch (XMLStreamException | FactoryConfigurationError e) {
			throw new IOException("XML parse exception");
		}
	}

	private DataNode parseXML() throws IOException {

		List<Long> reflList = null;
		long[] refl = null;
		String id = null;
		try {
			while (xmlReader.hasNext()) {

				int event = xmlReader.next();
				switch (event) {
					case XMLEvent.START_ELEMENT:
						switch (xmlReader.getLocalName()) {
							case "RESTRICTION_MAP":
								;
								id = xmlReader.getAttributeValue("", "ID");
								break;
							case "FRAGMENTS":
								reflList = new ArrayList<Long>();
								break;
							case "F":
								reflList.add(Long.parseLong(xmlReader.getAttributeValue("", "S")));
								break;
						}
						break;
					case XMLEvent.END_ELEMENT:
						switch (xmlReader.getLocalName()) {
							case "RESTRICTION_MAP":
								;
								return new DataNode(id, refl);
							case "FRAGMENTS":
								refl = ArrayUtils.toPrimitive(reflList.toArray(new Long[reflList.size()]));
								break;
							case "F":
								// Nothing to do
								break;
						}

						break;

					case XMLEvent.START_DOCUMENT:
						break;
					case XMLEvent.END_DOCUMENT:
						return null;

				}

			}

		} catch (XMLStreamException e) {
			throw new IOException("XML parse exception");
		}
		return null;
	}

	public LinkedHashMap<String, DataNode> readAllData() throws IOException {
		LinkedHashMap<String, DataNode> fragmentmap = new LinkedHashMap<String, DataNode>();
		DataNode fragment;
		do {
			fragment = read();
			if (fragment == null)
				break;
			else
				fragmentmap.put(fragment.name, fragment);
		} while (fragment != null);
		return fragmentmap;
	}

	public static LinkedHashMap<String, DataNode> readAllData(String filename) throws IOException {
		return readAllData(filename, -1);
	}

	public static LinkedHashMap<String, DataNode> readAllData(String filename, int format) throws IOException {
		OptMapDataReader omdr = new OptMapDataReader(filename, format);
		LinkedHashMap<String, DataNode> fragmentmap = omdr.readAllData();
		omdr.close();
		return fragmentmap;
	}

	public static LinkedHashMap<String, DataNode> readAllData(OptionSet options) throws IOException {
		OptMapDataReader omdr = new OptMapDataReader(options);
		LinkedHashMap<String, DataNode> fragmentmap = omdr.readAllData();
		omdr.close();
		return fragmentmap;
	}

	public static LinkedHashMap<Integer, LinkedHashMap<String, DataNode>> getLabelMap(OptionSet options) throws IOException {
		OptMapDataReader omdr = new OptMapDataReader(options);
		return getLabelMap(omdr);
	}

	public static LinkedHashMap<Integer, LinkedHashMap<String, DataNode>> getLabelMap(String filename, int format) throws IOException {
		OptMapDataReader omdr = new OptMapDataReader(filename, format);
		return getLabelMap(omdr);
	}

	private static LinkedHashMap<Integer, LinkedHashMap<String, DataNode>> getLabelMap(OptMapDataReader omdr) throws IOException {
		LinkedHashMap<Integer, LinkedHashMap<String, DataNode>> fragmentMapLabelMap = new LinkedHashMap<Integer, LinkedHashMap<String, DataNode>>();
		DataNode fragment;
		while ((fragment = omdr.read()) != null) {
			LinkedHashMap<String, DataNode> fragmentMap;
			if ((fragmentMap = fragmentMapLabelMap.get(fragment.getTotalSegment())) == null) {
				fragmentMap = new LinkedHashMap<String, DataNode>();
				fragmentMapLabelMap.put(fragment.getTotalSegment(), fragmentMap);
			}
			fragmentMap.put(fragment.name, fragment);
		}
		omdr.close();
		return fragmentMapLabelMap;

	}

	public static void assignOptions(ExtendOptionParser parser) {
		parser.addHeader("Data Reader Options", 1);
		parser.accepts("optmapin", "Optical map file input.").withRequiredArg().ofType(String.class);
		parser.accepts("optmapinformat", DataFormat.getFormatHelp()).withOptionalArg().ofType(Integer.class).defaultsTo(-1);
		parser.accepts("bnxsnr", "BNX SNR filter value").withOptionalArg().ofType(Double.class).defaultsTo(3.0);
	}

	public static int countData(String filename) throws IOException {
		OptMapDataReader omdr = new OptMapDataReader(filename);
		int totalCount = 0;
		while (omdr.read() != null)
			totalCount++;
		omdr.close();
		return totalCount;
	}
}