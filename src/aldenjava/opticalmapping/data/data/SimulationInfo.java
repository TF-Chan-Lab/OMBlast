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

import java.util.List;

import aldenjava.opticalmapping.GenomicPosNode;

/**
 * Stores the simulation information of a simulated optical map molecule
 * @author Alden
 *
 */
public class SimulationInfo {
	public final GenomicPosNode simuRegion;
	public final int simuStrand;
	public SimulationInfo(GenomicPosNode simuRegion, int simuStrand) {
		this.simuRegion = simuRegion;
		this.simuStrand = simuStrand;
	}

	public SimulationInfo(SimulationInfo simuInfo) {
		this.simuRegion = simuInfo.simuRegion;
		this.simuStrand = simuInfo.simuStrand;
	}

	public static boolean checkInfoValid(String ref, long start, long stop, int strand) {
		return (ref != null && !ref.isEmpty()) && (start >= 1) && (stop >= 1) && ((strand == 1) || (strand == -1));
	}
}




// For future development
class ReferenceSignal {
	String ref;
	int refpPos; 
}

// For future development
class VirtualSignal {
	// We can recover which signal is not retained (FN)
	// We can recover which signal is extra (FP)
	// We can recover measuring error

	// We can recover scaling error
	
	
	List<VirtualSignal> sources; 
	ReferenceSignal refSig = null;
	
	public boolean isFPSignal() {
		for (VirtualSignal vs : sources)
			if (!vs.isFPSignal())
				return false;
		return true;
	}
	
}
