SeaSurfaceReconstruction
=========================

INPUT:
---------------
configs.txt

=================================================================================================

PARAMETERS in configs.txt:

- mode(1): <string>				mode of processing
							(1) "evaluate"	(default)
							(2) "reconstruct"

- candidates: <string>				type of triangle candidates
							- "empty" for Delaunay triangles
							- "kODT" for k-order Delaunay triangles (default)

- ev: <list epochs>/all				training epochs
			
- on: <list epochs>/all				reconstruction epochs

- projection: <string>				projection of the coordinates of TGs and gridded points
							- "", no projection (default)
							- "LCC" for Lambert Conformal Conic projection

- store: <string>				new folder in experiments-folder

- tc: <int>					number of k for k-order Delaunay triangles, only if parameter "candidates" is "kODT"

- rt: <int>					number of months indicating the period of SSA reconstruction based on the training periodonly,
					        only for reconstruction mode,
						for example "-rt: 12" defines the climatological reconstruction

- stationTraining: 				path to the file that contains the tide gauge data for training, only for reconstruction mode

- stationPath: 					path to the file that contains the tide gauge data

- gridDirectory: 				path to the file or folder that contains the altimetry data

- boundingPath: 				optional, path to a file containing coordinates that bound the area of interest


NOTE (1): 
	The program attempts to find some identifier in every point file/for every month in a d35-file. 
	More specifically in non-d35-files it will search for the first "_" in a file name and everything after 
	that except the file extension will be used as a identifier. (e.g. a file named points_january_2018.csv will 
	have the identifier january_2018 in the outputs)
	In the case of a d35-file the months will be used as ids. 
	The identifiers are used to match points to their rasters. For a d35-file that means there has to be a 
	raster file for every month one wants to triangulate. 
	For non-d35-files that means that if there are more than one raster files, the point-files have to have 
	a name that connects them to their raster. Basically the id of a point file has to appear in one of the 
	raster files names. Raster-files and point-files do not have to have the exactly same id. The raster-file name
	just has to contain the id of a point-file. 
	If there is only one raster file every set of points will be processed on that raster.

NOTE (2): 	
	See examples in "configs.txt" for evaluation and "configsRec.txt" for reconstruction.

NOTE (3): 	
	Do not change the order of the last three parameters (stationPath, gridDirectory, boundingPath).
 					
=================================================================================================
							
OUTPUT: 

For each of the different output-file-types the program creates a seperate subfolder in results.

interpolations:
	trainingEpoch_reconstructionEpoch_o.nc: 	contains the SSA reconstruction archieved with min-error k-OD triangulation
	
	trainingEpoch_reconstructionEpoch_d.nc:		contains the SSA reconstruction archieved with Delaunay triangulation

reconstructions: (only in reconstruction mode)
	trainingEpoch_reconstructionEpoch_o.nc: 	contains the SSA reconstruction archieved with min-error k-OD triangulation
	
	trainingEpoch_reconstructionEpoch_d.nc:		contains the SSA reconstruction archieved with Delaunay triangulation

simpleAnomalies:
	trainingEpoch_reconstructionEpoch_o.nc: 	contains the misfit of computed SSA by ME k-ODT and reference SSA
	
	trainingEpoch_reconstructionEpoch_d.nc:		contains the misfit of computed SSA by Delaunay T and reference SSA
											
squaredAnomalies: 
	trainingEpoch_reconstructionEpoch_o.nc:		contains the squared misfit of computed SSA by ME k-ODT and reference SSA
											
	trainingEpoch_reconstructionEpoch_d.nc: 	contains the squared misfit of computed SSA by ME k-ODT and reference SSA 
										
statistics: 
	statisticsLog.txt: 				contains some statistics for every training and reconstruction epoch
	
triangulationCSVs:
	trainingEpoch_reconstructionEpoch_o.csv: 	contains triangles of with min-error k-OD triangulation
	
	trainingEpoch_reconstructionEpoch_d.csv: 	contains triangles of Delaunay triangulation
				
validations: 
	validationLog.txt: 				contains the validations of the computed/evaluated triangulations,
							if there is no such folder then all triangulations are validate

times:
	timeLog.txt: 					contains the runtime of each triangulation seperated into 
							formulation time and solve time
	
=================================================================================================

RUN:

	java -jar ssr.jar configs.txt

=================================================================================================

LICENSE:

This code is released under the GNU Public License (GPL) 3.0.

To read the GPL 3.0, read the file COPYING in this directory.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

=================================================================================================

DATA:

The raw data can be downloaded at [1] (altimetry data) and [2] (gauge data)
For the processed data that is input for the program and the result data contact A. Förster (foerster@igg.uni-bonn.de).


[1]ESA, 2020. Sea Level CCI ECV dataset: Time series of gridded sea level anomalies (sla). 
URL: https://catalogue.ceda.ac.uk/uuid/142052b9dc754f6da47a631e35ec4609,
doi:10.5270/esa-sea\_level\_cci-msla-1993\_2015-v\_2.0-201612.
European Space Agency (ESA).

[2]PSMSL, 2020. Tide gauge data. Retrieved 18 Jan 2020
from http://www.psmsl.org/data/obtaining/, Tech. Rep. ,
Permanent Service for Mean Sea Level (PSMSL).

