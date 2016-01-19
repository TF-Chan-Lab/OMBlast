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
import java.io.InputStream;
import java.util.LinkedHashMap;

import joptsimple.OptionSet;
import aldenjava.opticalmapping.data.DataFormat;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;
/**
 * A wrapper class redirecting to <code>OptMapDataReader</code> with options changed to refmapin
 * @author Alden
 *
 */
public class ReferenceReader extends OptMapDataReader {

	public ReferenceReader(OptionSet options) throws IOException {
		this((String) options.valueOf("refmapin"), (int) options.valueOf("refmapinformat"));
	}

	public ReferenceReader(String filename, int format) throws IOException {
		super(filename, format);
	}

	public ReferenceReader(String filename, DataFormat dformat) throws IOException {
		super(filename, dformat);
	}

	public ReferenceReader(InputStream stream, DataFormat dformat) throws IOException {
		super(stream, dformat);
	}

	public static LinkedHashMap<String, DataNode> readAllData(String filename) throws IOException {
		return readAllData(filename, -1);
	}

	public static LinkedHashMap<String, DataNode> readAllData(String filename, int format) throws IOException {
		ReferenceReader rr = new ReferenceReader(filename, format);
		LinkedHashMap<String, DataNode> optrefmap = rr.readAllData();
		rr.close();
		return optrefmap;
	}

	public static LinkedHashMap<String, DataNode> readAllData(OptionSet options) throws IOException {
		ReferenceReader rr = new ReferenceReader(options);
		LinkedHashMap<String, DataNode> optrefmap = rr.readAllData();
		rr.close();
		return optrefmap;
	}

	public static void assignOptions(ExtendOptionParser parser) {
		parser.addHeader("Reference Reader Options", 1);
		parser.accepts("refmapin", "reference map file").withRequiredArg().ofType(String.class);
		parser.accepts("refmapinformat", DataFormat.getFormatHelp()).withOptionalArg().ofType(Integer.class).defaultsTo(-1);
	}

}
