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

import java.util.HashMap;
import java.util.Map;

import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.commons.io.FilenameUtils;

import aldenjava.opticalmapping.miscellaneous.InvalidFileFormatException;

/**
 * Optical map alignment result format that can be parsed
 * 
 * @author Alden
 *
 */
public enum ResultFormat {
	OMA (0, "OM Alignment Format (OMA)", "oma"),
	OMD (1, "OM Detailed Alignment Format (OMD)", "omd"),
	XMAP (2, "XMAP format (XMAP)", "xmap"),
	VAL (3, "Valouev et al. format", "val"),
	SOMA (4, "SOMA v2 Unique Match Format", "somav2"),
	PSL (5, "Twin PSL Format", "psl");

	private final int format;
	private final String description;
	private final String extension;

	ResultFormat(int format, String description, String extension) {
		this.format = format;
		this.description = description;
		this.extension = extension;
	}

	public int getFormat() {
		return format;
	}

	public String getExtension() {
		return extension;
	}

	public static String[] getExtensions() {
		String[] extensions = new String[ResultFormat.values().length];
		int index = 0;
		for (ResultFormat rformat : ResultFormat.values())
			extensions[index++] = rformat.getExtension();
		return extensions;
	}

	public static final Map<Integer, ResultFormat> lookupmap = new HashMap<Integer, ResultFormat>();
	static {
		for (ResultFormat rformat : ResultFormat.values())
			lookupmap.put(rformat.getFormat(), rformat);
	}
	public static final Map<String, ResultFormat> lookupfileextmap = new HashMap<String, ResultFormat>();
	static {
		for (ResultFormat rformat : ResultFormat.values())
			lookupfileextmap.put(rformat.getExtension(), rformat);
	}

	public static final boolean isValidFormat(String extension) {
		return lookupfileextmap.containsKey(extension);
	}

	public static final ResultFormat lookup(String path, int format) {
		if (format == -1)
			return lookupfileext(FilenameUtils.getExtension(path));
		if (!lookupmap.containsKey(format))
			throw new InvalidFileFormatException();
		return lookupmap.get(format);
	}

	public static final ResultFormat lookup(int format) {
		if (!lookupmap.containsKey(format))
			throw new InvalidFileFormatException();
		return lookupmap.get(format);
	}

	public static final ResultFormat lookupfileext(String extension) {
		if (!lookupfileextmap.containsKey(extension))
			throw new InvalidFileFormatException();
		return lookupfileextmap.get(extension);
	}

	public static String getFormatHelp() {
		StringBuilder resultFormatHelp = new StringBuilder();
		resultFormatHelp.append(String.format("%d: %s; ", -1, "Auto-detected from file extension"));
		for (ResultFormat rformat : ResultFormat.values())
			resultFormatHelp.append(String.format("%d: %s; ", rformat.format, rformat.description));
		return resultFormatHelp.toString();
	}

	public static final FileNameExtensionFilter[] getFileNameExtensionFilter(boolean allowAllSupport) {
		FileNameExtensionFilter[] filters = new FileNameExtensionFilter[ResultFormat.values().length + (allowAllSupport ? 1 : 0)];
		int index = 0;
		if (allowAllSupport)
			filters[index++] = new FileNameExtensionFilter("All supported format", ResultFormat.getExtensions());
		for (ResultFormat rformat : ResultFormat.values())
			filters[index++] = new FileNameExtensionFilter(rformat.description, rformat.extension);
		return filters;
	}

}
