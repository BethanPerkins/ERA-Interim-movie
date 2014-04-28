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
import netCDF4

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

    # hard-coded definition: the directory that the data is stored in and the directory that the movie's are written to.
    dataFilepath = "/media/Data_/ERA_Interim"
    movieFilepath = "/home/assimila/Videos"

    # do some basic checking of the vlaues. Outsource this to function?
    #
    # 
    # Is latitudes[0] smaller than latitudes[1]?
    # Is longitudes[0] smaller than longitudes[1]?
    # Is level one of the available pressure levels
    # are the times relistic. does the first time come before the last?
    #

    #collect user-defined values and defaults together into a Query object
    thisQuery = Query(variable, start, end, resolution, latitudes, longitudes, level)

    #set the filename for this particular data query
    EIFilename = buildEIFilename(thisQuery, dataFilepath, ".nc")

    #if this file does not exist then download it
    if not os.path.isfile(EIFilename):
        downloadERAInterim(thisQuery, EIFilename)

    #open the file and pull out the required information. Bundle this all up into an object?
    # info that we need:
    # variable data to sit in the movie
    # units
    # axes values (lat, long and time).
    movieData = loadMovieData(thisQuery, EIFilename)


    #make a movie of that information. This should call the existing MovieData object, so that the user can re-define thing in here on the fly.
    #makeMovie(movieData)    
    
    #return out



	
def buildEIFilename(thisQuery,filepath, extension):
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
 	filename+=extension

	return filename

def downloadERAInterim(thisQuery, EIFilename):
    """
         Create eraInterimDictionary
         The eraInterimDictionary holds supplementary information on each of the supported variables. 
         This provides the information that the ECMWF api needs in order to access the archive. 
         For further information on these and many more variables go to http://www.ecmwf.int/publications/manuals/d/gribapi/param
    """
    eraInterimDictionary = {
	    "temperature": (130, "pl", "an", "t"),
	    "precipitation": (228, "sfc", "fc", "tp"),
	    "geopotential": (129, "pl", "an", "z"),
	    "uWind": (131, "pl", "an", "u"),
	    "vWind": (132, "pl", "an", "v"),
	    "specificHumidity": (133, "pl", "an", "q"),
	    "cloudCover": (164, "sfc", "an", "tcc"),
	    "relativeVorticity": (138, "pl", "an", "vo"),
	    "soilTemperature": (139, "sfc", "an", "stl1"),
	    "snowDepth": (141, "sfc", "an", "sd"),
	    "menaSeaLevelPressure": (151, "sfc", "an", "msl"),
	    "evaporation": (182, "sfc", "fc", "e"),
     } 

    thisQuery.addDictionaryData(eraInterimDictionary["temperature"])

    thisQuery.getDataECMWF(EIFilename)




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


def loadMovieData(query, dataFile):
    """
        this function sets up the MovieData object and puts information from the netcdf file into this object.
    """
    
    movieData = MovieData

    ncobject = netCDF4.Dataset(dataFile)

    movieData.variable = ncobject

    
    print "ping"    

    return "this will be movie data"

# load data extract data and captions etc from the NetCDF file
# inputs - filename
# process: open file, define new object, populate new object with information from the file.
# outputs - object with the required data (lat, long, time), dates (time), lat extent (lat), lon extent (lon), Name (long_name), units (units)


#build a simple movie from it 
# inputs data object
# process: define the size and shape of the movie frame etc. Unpack the data into that framework
# output: movie file (play movie also)

def MakeMovie(movieData):
    """
        Build a simple movie using the movieDataObject
    """
    # things to think about
    # wrapping longitudes..? Nah! we're only producing a flat map
    #

#define the user's required data block
#inputs - variable (default: temperature), startdate (default:today-1 year-30days), enddate (default:today- 1year), lat-lon resolution (default: 1degree), lat extents (default:-90, 90), -lon extents (-180, 180), level (default:0)
# This class holds the user definition of a geospatial data block.
class Query:

    """
        Holds the user definition of the required paramter and required geotemporal data block size. 
        Has a method to download data from ECMWF servers using this information
    """

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

    def addDictionaryData(self, dictionaryData):
        self.param = dictionaryData[0]
        self.levtype = dictionaryData[1]
        self.type = dictionaryData[2]
        self.shortName = dictionaryData[3]

    def getDataECMWF(self, filename):
	
	    from ecmwfapi import ECMWFDataServer
	    server = ECMWFDataServer()

	    #build a dictironary that provdes the MARS like information to the ECMWF servers
	    # start with fixed fields that are the same for all queries
	    self.MARSCode = {
		    'stream'  : "oper",
	        'dataset' : "interim",
	        'param'   : str(self.param),
	        'date'    : self.start+"/to/"+self.end,
	        'area'    : str(self.latitude_max)+"/"+str(self.longitude_min)+"/"+str(self.latitude_min)+"/"+str(self.longitude_max),
	        'time'    : "00/12",
	        'grid'    : str(self.resolution)+"/"+str(self.resolution),
	        'target'  : filename,
	        'levtype' : str(self.levtype),
	        'type'    : str(self.type),
            'format'  : "netcdf"
	    }

	    if dictionaryData[2] is "fc":
		    self.MARSCode['step'] = "12"

	    elif dictionaryData[2] is "an":
		    self.MARSCode['step'] = '0'

		    if dictionaryData[1] == "pl":
			    self.MARSCode['levelist'] = str(self.level)

	    ##retrieve the data
	    server.retrieve(self.MARSCode)



class MovieData:
    """
        Holds the data required to construct a movie using latitude-longitude data and animated in time
    """
    def __init__(self):

        # import the required libraries
        import numpy as np
        import matplotlip.pyplot as plt

        # define empty containers 
        self.variable = np.array([])
        self.variableUnits = ""
        self.time = np.array([])
        self.latitude = np.array([])
        self.longitude = np.array([])
        self.colourRange = np.array([])

        #define defaults
        self.colourmap = plt.cm.gist_ncar



    def checkMovieData(self):
        """
            Method which checks that all the essential Movie Data
        """
        print "checkMovieData still to write"
        #write this



    def defineFilename(self, query, path, extension):
        """
            Build the filename of the movie file
        """
        self.filename = buildEIFilename(query, path, ".mov")



    def setColourmap(self, newMapName):
        """
            Defines a new colormap, newMapName must be the name of a matplotlib colormap. 
            All colormaps can be seen here: http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
        """
        self.colourmap = plt.get_cmap(newMapName)



    def setColourRange(self, colourbarMinimum, colourbarMaximum):
        """
            Manually set the colourbar maximum and minimum values
        """
        self.colourRange = np.array([colorbarMinimum, colourbarMaximum])


    


##############################################################################
out = plotEIMovie()

#print out

