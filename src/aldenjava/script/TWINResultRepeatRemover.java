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


package aldenjava.script;

import java.io.IOException;
import java.util.List;

import joptsimple.OptionSet;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultNode;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultReader;
import aldenjava.opticalmapping.data.mappingresult.OptMapResultWriter;
import aldenjava.opticalmapping.miscellaneous.ExtendOptionParser;

public class TWINResultRepeatRemover {

	public static void main(String[] args) throws IOException {
		ExtendOptionParser parser = new ExtendOptionParser(TWINResultRepeatRemover.class.getSimpleName());
		OptMapResultReader.assignOptions(parser);
		OptMapResultWriter.assignOptions(parser);
		if (args.length == 0) {
			parser.printHelpOn(System.out);
			return;
		}
		OptionSet options = parser.parse(args);
		List<OptMapResultNode> resultlist;
		
		OptMapResultReader omrr = new OptMapResultReader(options);
		OptMapResultWriter omrw = new OptMapResultWriter(options);
		while ((resultlist = omrr.readNextList()) != null) {
			
			if (!resultlist.get(0).isUsed()) {
				omrw.write(resultlist.get(0));
				continue;
			}
			boolean overlap = true;
			LOOP:
			for (int i = 0; i < resultlist.size(); i++)
				for (int j = i + 1; j < resultlist.size(); j++) {
					if (!resultlist.get(i).mappedRegion.isClose(resultlist.get(j).mappedRegion, 0)) {
						overlap = false;
						break LOOP;
					}
				}
			if (overlap) {
				OptMapResultNode bestResult = resultlist.get(0);
				for (OptMapResultNode r : resultlist)
					if (r.mappedscore < bestResult.mappedscore)
						bestResult = r;
				omrw.write(bestResult);
			}
		}
		omrr.close();
		omrw.close();
	}

}
