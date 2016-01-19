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

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
/**
 * Abstract class for reading data of any type from a file or input stream using a buffered reader. The only method that a subclass must implement is read(). Subclasses can override commentReader() and proceedNextLine() based on the file format 
 * 
 * @param <T>	the type of elements to be read
 * 
 * @see	OMWriter<T>
 * 
 * @author Alden
 *
 * 
 */
public abstract class OMReader<T> implements Closeable {

	protected final BufferedReader br;
	protected String nextline;
	
	/**
	 * Creates a new reader on a file
	 * @param stream
	 * @throws IOException
	 */
	public OMReader(String filename) throws IOException {
		br = new BufferedReader(new FileReader(filename));
		commentReader();
	}
	/**
	 * Creates a new reader on an <code>InputStream</code>
	 * @param stream
	 * @throws IOException
	 */
	public OMReader(InputStream stream) throws IOException {
		br = new BufferedReader(new InputStreamReader(stream));
		commentReader();
	}
	/**
	 * Attempts to skip the headers in the data file. Subclasses can override this method depending on file format  
	 * @throws IOException
	 */
	protected void commentReader() throws IOException {
		do {
			nextline = br.readLine();
		}
		while (nextline != null && (nextline.startsWith("#") || nextline.isEmpty()));
	}
	/**
	 * Attempts to proceed to next line to read data. It is recommended to override and use this method to proceed to next line 
	 * @throws IOException	if an I/O error occurs
	 */
	protected void proceedNextLine() throws IOException	{
		do {
			nextline = br.readLine();
		}
		while (nextline != null && (nextline.startsWith("#") || nextline.isEmpty()));
	}
	/**
	 * Attempts to read for next data entry.  
	 * @return	the data of type <code>T</code> or <code>null</code> if the end of file is reached
	 * @throws IOException	if an I/O error occurs
	 */
	public abstract T read() throws IOException;
	/**
	 * Attempts to read for all data and output as list.   
	 * @return	the list of data of type <code>T</code> or an empty list if the there is no entry
	 * @throws IOException	if an I/O error occurs
	 */
	public List<T> readAll() throws IOException {
		List<T> tList = new ArrayList<T>();
		T t;
		while ((t = read()) != null)
			tList.add(t);
		return tList;
	}
//	public HashSet<String> readAllIdentifier() throws IOException
//	{
//		HashSet<String> identifiers = new HashSet<String>();
//		T t;
//		while ((t = read()) != null)
//		{
//			if (!(t instanceof Identifiable))
//				throw new UnsupportedOperationException();
//			identifiers.add(((Identifiable) t).getIdentifier());
//		}
//		return identifiers;
//		
//	}	
	@Override
	public void close() throws IOException {
		br.close();
	}

}


//interface Identifiable
//{
//	public String getIdentifier();
//}