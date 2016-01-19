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


package aldenjava.opticalmapping.miscellaneous;

import joptsimple.OptionSet;

/**
 * A class implements the <code>SelectableMode</code> interface to indicate it has various modes with different functions.
 * 
 * @author Alden
 *
 */
public interface SelectableMode {

	/**
	 * Select the mode according to <code>options</code>. 
	 * 
	 * @param options
	 * @throws IllegalArgumentException
	 *             if the selected mode is invalid
	 */
	public void setMode(OptionSet options) throws IllegalArgumentException;

	/**
	 * Select the mode according to <code>mode</code>
	 * 
	 * @param mode
	 *            the selected mode
	 * @throws IllegalArgumentException
	 *             if the selected mode is invalid
	 */
	public void setMode(int mode) throws IllegalArgumentException;

	/**
	 * Returns the selected mode
	 * 
	 * @return the selected mode
	 * @throws IllegalStateException
	 *             if the mode is not set or the selected mode is invalid. If a class that implements this interface has default mode, it should not throw this exception.
	 */
	public int getMode() throws IllegalStateException;

}
