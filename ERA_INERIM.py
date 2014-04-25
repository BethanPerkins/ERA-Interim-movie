#
#
# This code will create a movie of ERA Interim data through python.
#
#
# Built by bethan.perkins@assimila.eu May 2014
#
#

# Import libraries
import datetime
import os.path

def plotEIMovie(variable="temperature",start="2013-01-01", end="2013-01-31", resolution="1", latitudes=(-90, 90), longitudes=(-180, 180), level=1000):
	"""
	 The variables available are:
		-temperature
		-precipitation
		-geopotential
		-uWind
		-vWind
		-specificHumidity
		-cloudCover
		-relativeVorticity
		-soilTemperature
		-snowDepth
		-menaSeaLevelPressure
		-evaporation
	
	 Dates need to be defined using ISO 8601 Standard i.e. YYYY-MM-DD
	 (see http://en.wikipedia.org/wiki/ISO_8601 for more infomration on the standard)
	 In this syntax, the 30th of April, 2012 will be written 2012-04-30
	
	 Resolutions are given as fractions of a degree, available values for ERA Interim are:
		0.125, 0.25, 0.5, 0.75, 1.0, 1.125, 1.5, 2.0, 2.5, 3.0
	
	 latitudes and longitudes are defined as a tuple - (x,y) where x is the minimum latitude/longitude and y is the maximum.
	"""

	# hard-coded definition: the directory that the data is stored in. Global variable.
	filepath = "C:/Users/Bethan/Data/ERA_Interim"

	# do some basic checking of the vlaues.
	#
	# 
	# Is latitudes[0] smaller than latitudes[1]?
	# Is longitudes[0] smaller than longitudes[1]?
	# Is level one of the available pressure levels
	#

	#collect user-defined values and defaults together into a Query object
	thisQuery = Query(variable, start, end, resolution, latitudes, longitudes, level)

	#set the filename for this particular data query
	EIFilename = buildEIFilename(thisQuery, filepath)

	#if this file does not exist then download it
	if not os.path.isfile(EIFilename):
		out = downloadERAInterim(thisQuery, EIFilename)

	# check again that the file is now there
	#if not os.path.isfile(EIFilename):
		#some error.

	return out



	
def buildEIFilename(thisQuery,filepath):
	"""
	This function takes the Query object and the pre-defined filepath and creats a filename for the data relating to that query. 

	This is the file that subsequent searches will look for when checking to see if the data has been downloaded already.
	"""
	#start by setting the directory
	filename = filepath

	# add a slash if required
	if filename[-1] is not "/" or filename[-1] is not "\\":
		filename+="/"

	#add the common identifier "EI" to show that all files are ERA Interim data
 	filename+="EI_"

 	#add in the variable name
 	filename+=thisQuery.variable+"_"

 	#add in the level
 	filename+=str(thisQuery.level)+"hPa_"

 	#add in the lat-long
 	filename+="lats"+str(thisQuery.latitude_min)+str(thisQuery.latitude_max)
 	filename+="_"
 	filename+="lons"+str(thisQuery.longitude_min)+str(thisQuery.longitude_max)
	filename+="_"

 	#add in the time
 	#filename+=datetime.date.strftime(thisQuery.start,"%Y%m%d")
 	filename+=thisQuery.start
 	filename+="_to_"
 	#filename+=datetime.date.strftime(thisQuery.end,"%Y%m%d")
 	filename+=thisQuery.end

 	#add the netcdf extension
 	filename+=".nc"

	return filename

def downloadERAInterim(thisQuery, EIFilename):

	#
	# Create eraInterimDictionary
	# The eraInterimDictionary holds supplementary information on each of the supported variables. This provides the information that the ECMWF api needs in order to access the archive. For further information on these and many more variables go to http://www.ecmwf.int/publications/manuals/d/gribapi/param
	#
	eraInterimDictionary = {
		"temperature": (130, "pl", "an"),
		"precipitation": (228, "sfc", "fc"),
		"geopotential": (129, "pl", "an"),
		"uWind": (131, "pl", "an"),
		"vWind": (132, "pl", "an"),
		"specificHumidity": (133, "pl", "an"),
		"cloudCover": (164, "sfc", "an"),
		"relativeVorticity": (138, "pl", "an"),
		"soilTemperature": (139, "sfc", "an"),
		"snowDepth": (141, "sfc", "an"),
		"menaSeaLevelPressure": (151, "sfc", "an"),
		"evaporation": (182, "sfc", "fc"),
	 } 

	out= thisQuery.getDataECMWF(eraInterimDictionary["temperature"], EIFilename)

	return out

# --> called only if the output to previous is false
# download it if it doesn't exist
# inputs - user-block object, filename
# pre-defined: a database of the variable numbers and other preripheral information required in order to download the data
# output - none

# load a dictionary containing references to each of the available datasets
# inputs -> search terms
# outputs -> dictionary entry for this variable defining: 
# stream: oper 
# levtype: sfc (surface), pl (pressure levels)
# type: fc (forecast), an (analysis) 
#def searchDictionary(search ):




# load data extract data and captions etc from the NetCDF file
# inputs - filename
# process: open file, define new object, populate new object with information from the file.
# outputs - object with the required data (lat, long, time), dates (time), lat extent (lat), lon extent (lon), Name (long_name), units (units)


#build a simple movie from it 
# inputs data object
# process: define the size and shape of the movie frame etc. Unpack the data into that framework
# output: movie file (play movie also)



#define the user's required data block
#inputs - variable (default: temperature), startdate (default:today-1 year-30days), enddate (default:today- 1year), lat-lon resolution (default: 1degree), lat extents (default:-90, 90), -lon extents (-180, 180), level (default:0)
# This class holds the user definition of a geospatial data block.
class Query:
	"holds the user definition of the required paramter and required geotemporal data block size."
	def __init__(self, variable, start, end, resolution, latitudes, longitudes, level):
		self.variable = variable
		self.start = start #datetime.datetime.strptime(start,"%Y-%m-%d")
		self.end = end #datetime.datetime.strptime(end,"%Y-%m-%d")
		self.resolution = resolution
		self.latitude_min = latitudes[0]
		self.latitude_max = latitudes[1]
		self.longitude_min = longitudes[0]
		self.longitude_max = longitudes[1]
		self.level = level

	def getDataECMWF(self, dictionaryData, filename):
		
		#from ecmwfapi import ECMWFDataServer
		#server = ECMWFDataServer()
  
		#build a dictironary that provdes the MARS like information tothe ECMWF servers
		# start with fixed fields that are the same for all queries
		MARSCode = {
			'stream'  : "oper",
		    'dataset' : "interim",
		    'param'   : str(dictionaryData[0]),
		    'date'    : self.start+"/to/"+self.end,
		    'area'    : str(self.latitude_max)+"/"+str(self.longitude_min)+"/"+str(self.latitude_min)+"/"+str(self.longitude_max),
		    'time'    : "00/12",
		    'grid'    : str(self.resolution)+"/"+str(self.resolution),
		    'target'  : filename
		    'levtype' : str(dictionaryData[1])
		    'type'    : str(dictionaryData[0])
		    'number'  : "all"
		}

		if dictionarydata[2] is "fc":
			MARSCode['step'] = "12"

		elif dictionaryData[2] is "an":
			MARSCode['step'] = '0'

			if dictionaryData[1] is "pl":
				MARSCode['level'] = str(self.level)

    	# 

		##retrieve the data
		#server.retrieve(MARScode)

		return MARSCode


##############################################################################
out = plotEIMovie()

print out

