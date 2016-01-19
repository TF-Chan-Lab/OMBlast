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


package aldenjava.opticalmapping.data;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * @author Alden
 * Base class for all Optical Mapping project writer.
 * Wrapping a buffered writer
 * @param <T>
 */
public abstract class OMWriter<T> implements Closeable{
	// Format handling is bad: initialize header fails in constructor
	private static String commands = null;
	protected BufferedWriter bw;
	public OMWriter(String filename) throws IOException {
		this(filename, true);
	}	

	public OMWriter(String filename, boolean autoInitializeHeader) throws IOException {
		bw = new BufferedWriter(new FileWriter(filename));
		if (autoInitializeHeader)
			initializeHeader();
	}
	
	protected void initializeHeader() throws IOException
	{
		// Unimplemented in any classes
		if (OMWriter.commands != null)
			bw.write("#" + OMWriter.commands + "\n");
		// Plan: Move this into a new method to be called before initialize header
		// Plan A: Use implement interface to state if a subclass support writing command (e.g. fasta writer does not support)
		// Plan B: Use an abstract method returning boolean to check if a subclass support
	}
	public abstract void write(T t) throws IOException;
			
	public void writeAll(List<T> tList) throws IOException
	{
		if (tList == null)
			throw new NullPointerException("tList");
		for (T t : tList)
			write(t);
	}
	public void writeAll(LinkedHashMap<?, T> tMap) throws IOException
	{
		for (T t : tMap.values())
			write(t);
	}
	@Override
	public final void close() throws IOException
	{
		bw.close();
	}
	public static void setCommands(String s)
	{
		OMWriter.commands = s;
	}
}
