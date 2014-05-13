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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
from ecmwfapi import ECMWFDataServer


#########################################################################################
# OTHER DEPENDENCIES ##
#
# - movie-building function (makeMovie) requires FFMpeg to be working on the machine.
#
#
#########################################################################################

def plotEIMovie(variable="temperature",start="2013-01-01", end="2013-01-31", resolution="1", latitudes=(-90, 90), longitudes=(-180, 180), level=1000):
    """
        Function to download data from the ECMWF ERA Interim data archive and render as a movie.

        plotEIMovie(variable=[variable name], start=[movie start date], 
        
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
        -meanSeaLevelPressure
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
    print "data directory is "+dataFilepath
    movieFilepath = "/home/assimila/Videos"
    print "movie will be created in "+movieFilepath

    # check that the variable is equal to one of the allowed variable
    allowedVariables = ("temperature", "precipitation","geopotential","uWind","vWind","specificHumidity","cloudCover","relativeVorticity","soilTemperature","snowDepth","meanSeaLevelPressure","evaporation")
    if variable not in allowedVariables:
        raise ValueError("variable must be one of the available values: "+str(allowedValues))

    # Two values of longitude
    if len(longitudes) is not 2:
        raise ValueError("please provide two values of longitude: (minimum, maximum)")

    # Negative westerly longitudes
    if max(longitudes) > 180:
        raise ValueError("westerly longitudes should be input as negative values. For example, use '-90' instead of '270'")

    # Westerly longitude comes first
    if latitudes[0] > latitudes[1]:
        raise ValueError("first value of longitude must be the most westerly value (we don't wrap around the dateline)")

    # Two values of latitude
    if len(latitudes) is not 2:
        raise ValueError("please provide two values of latitude: (minimum, maximum)")

    # Sensible values of latitude
    if latitudes[0] > latitudes[1]:
        raise ValueError("first value of latitude must be the most the most southerly value")

    # Just one value for level
    try:
        level = int(level)
    except TypeError:
        raise ValueError("please provide only one value for level")

    # Level one of the available levels
    availableLevels=(1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1)
    if level not in availableLevels:
        raise ValueError("level must be one of: "+ str(availableLevels[:]))

    # Just one value of resolution
    try:
        resolution = float(resolution)
    except TypeError:
        raise ValueError("please provide only one value for resolution") 

    # Resolution one of the available resolutions
    availableResolutions=(3.0, 2.5, 2.0, 1.5, 1.125, 1.0, 0.75, 0.5, 0.25, 0.125)
    if resolution not in availableResolutions:
        raise ValueError("resolution must be one of the available values: "+ str(availableResolutions[:]))    

    # check that the start date is a real date
    try:
        startDate = datetime.datetime.strptime(start, "%Y-%m-%d")
    except ValueError:
        raise ValueError("Start date is invalid. Dates need to be defined using ISO 8601 Standard i.e. YYYY-MM-DD")
        
    # check that the end date is a real date
    try:
        endDate = datetime.datetime.strptime(end, "%Y-%m-%d")
    except ValueError:
        raise ValueError("End date is invalid. Dates need to be defined using ISO 8601 Standard i.e. YYYY-MM-DD")

    #check that star comes beofre end
    if startDate-endDate > datetime.timedelta(0):
        raise ValueError("start date must come before end date")

    #collect user-defined values and defaults together into a Query object
    thisQuery = Query(variable, start, end, resolution, latitudes, longitudes, level)

    #set the filename for this particular data query
    EIFilename = buildEIFilename(thisQuery, dataFilepath, ".nc")

    #if this file does not exist then download it
    if not os.path.isfile(EIFilename):
        #download the data
        thisQuery.getDataECMWF(EIFilename)        

    #define the movie filename
    movieFilename = buildEIFilename(thisQuery, movieFilepath, ".mp4")

    #open the file and pull out the required information.
    movieData = loadMovieData(thisQuery, EIFilename, movieFilename)

    #make a movie 
    makeMovie(movieData)    




	
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
 	filename+=thisQuery.start
 	filename+="_to_"
 	filename+=thisQuery.end

 	#add the netcdf extension
 	filename+=extension

	return filename


def loadMovieData(query, dataFile, filename):
    """
        this function sets up the MovieData object and puts information from the netcdf file into this object.
        Open the file and pull out the required information. Bundle this all up into an object?
        info that we need:
            variable data to sit in the movie
            units
            axes values (lat, long and time).

    """
    # Create a new MovieData object
    movieData = MovieData()
    
    # Open the netCDF file
    ncobject = netCDF4.Dataset(dataFile)

    # extract variable/parameter date (values and units)
    movieData.variable = ncobject.variables[query.shortName][:]
    movieData.variableName = ncobject.variables[query.shortName].long_name
    movieData.variableUnits = ncobject.variables[query.shortName].units

    # extract time data (values and units). Convert values into datetime objects.
    #movieData.time = ncobject.variables['time'][:]
    relativeTimeUnits = ncobject.variables['time'].units
    relativeTimes = ncobject.variables['time'][:]
    movieData.time = netCDF4.num2date(relativeTimes, units=relativeTimeUnits, calendar="standard")
    

    # extract latitude and longitude data
    movieData.latitude = ncobject.variables['latitude'][:]
    movieData.longitude = ncobject.variables['longitude'][:]

    # set colour range to be equal to whole variable data range.
    movieData.setColourRange(np.min(movieData.variable), np.max(movieData.variable))

    # set output filename
    movieData.filename = filename

    #close the netCDF file
    ncobject.close()

    return movieData

def makeMovie(movieData):
    """
        Build movie using the movieData object
    """

    # Build the figure environment   
    fig = plt.figure(figsize=(12.0, 6.75))
    ax = fig.add_axes([0.1,0.1,0.85,0.75])

    # Build the foundation map using Basemap  
    m = Basemap(llcrnrlat=movieData.latitude[-1], llcrnrlon=movieData.longitude[0], urcrnrlat=movieData.latitude[0],urcrnrlon=movieData.longitude[-1], projection="cyl", lat_0 = 0, lon_0 = 0, resolution = 'l')    

    # Turn lat-long data into a grid and translate this for the map projection.  
    LONG, LAT = np.meshgrid(movieData.longitude, movieData.latitude)
    longitude,latitude = m(LONG, LAT) 

    #create a function which takes a value which takes a value of time and updates the figure with the next frame. This is used by matplotlib.animation (imported as "animation")
    def makeMovieFrames(time):
        global plot
        ax.cla()
        plot = m.contourf(longitude, latitude, movieData.variable[time,:,:], movieData.contours, cmap=movieData.colourmap, extend="both")
        m.drawcoastlines(linewidth=0.5)

        # draw parallels and meridians
        pll = m.drawparallels(np.arange(-90,90,20),labels=[1,1,0,1], latmax=90, color="white")
        mrd = m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1], color="white")
        
        # add in the colorbar and label it with the units of the variable
        cbar = m.colorbar(plot,location='bottom',pad="20%")
        cbar.set_label(movieData.variableUnits, color="white")

        #set the parallels, meridians and colorbar labels to be white
        setLabelsColor(pll, "white")
        setLabelsColor(mrd, "white")
        for t in cbar.ax.get_xticklabels():
            t.set_color("white") 

        # annotate the plot with the time and variale name
        plt.annotate(datetime.datetime.strftime(movieData.time[time], "%Y %m %d  %HZ"), xy=(-0.2, 1.1), xycoords='axes fraction', color="white", fontsize=16)
        plt.annotate(movieData.variableName, xy=(0.9, 1.1), xycoords='axes fraction', color="white", fontsize=16)

    #draw the first frame
    makeMovieFrames(0)

    #plt.savefig(movieData.filename+".png", facecolor="black")

    #use this function to build the movie
    movie = animation.FuncAnimation(fig, makeMovieFrames, frames=movieData.time.shape[0], interval=20, blit=False)

    #define the writer for the movie
    writer = animation.FFMpegWriter(fps=6, bitrate=10000)

    #save the movie
    movie.save(movieData.filename, codec="mp4", writer=writer, dpi=500, extra_args=['-vcodec', 'libx264'], savefig_kwargs={'facecolor':'black'}) 


def setLabelsColor(objectOutputDict, color):
    """
        A function to set the color of plot labels created by "drawmeridians" and "drawparallels"
        modified from: https://github.com/matplotlib/basemap/issues/145
    """
    for key in objectOutputDict:
        for label in objectOutputDict[key][1]:
            label.set_color(color)

def eraInterimDictionary():
    """
         The eraInterimDictionary holds supplementary information on each of the supported variables. 
         This provides the information that the ECMWF api needs in order to access the archive. 
         For further information on these and many more variables go to http://www.ecmwf.int/publications/manuals/d/gribapi/param
    """
    dictionary = {"temperature": (130, "pl", "an", "t"),
    "precipitation": (228, "sfc", "fc", "tp"),
    "geopotential": (129, "pl", "an", "z"),
    "uWind": (131, "pl", "an", "u"),
    "vWind": (132, "pl", "an", "v"),
    "specificHumidity": (133, "pl", "an", "q"),
    "cloudCover": (164, "sfc", "an", "tcc"),
    "relativeVorticity": (138, "pl", "an", "vo"),
    "soilTemperature": (139, "sfc", "an", "stl1"),
    "snowDepth": (141, "sfc", "an", "sd"),
    "meanSeaLevelPressure": (151, "sfc", "an", "msl"),
    "evaporation": (182, "sfc", "fc", "e"),
    }

    return dictionary


class Query:

    """
        Holds the user definition of the required paramter and required geotemporal data block size. 
        Has a method to download data from ECMWF servers using this information
    """

    def __init__(self, variable, start, end, resolution, latitudes, longitudes, level):
        self.variable = variable
        self.start = start
        self.end = end
        self.resolution = resolution
        self.latitude_min = latitudes[0]
        self.latitude_max = latitudes[1]
        self.longitude_min = longitudes[0]
        self.longitude_max = longitudes[1]
        self.level = level

        dictionary = eraInterimDictionary()

        self.param = dictionary[variable][0]
        self.levtype = dictionary[variable][1]
        self.type = dictionary[variable][2]
        self.shortName = dictionary[variable][3]

    def getDataECMWF(self, filename):

        # create a new ECMWFServer object
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

        # if the filed is a forecast, then we want the 12-hour cumulated values.
	    if self.type is "fc":
		    self.MARSCode['step'] = "12"

        # if the field is analysis then we want the instantaneous values
	    elif self.type is "an":
		    self.MARSCode['step'] = '0'
            # if the analysed field is also on pressure levels, then define the level that we want.
		    if self.levtype is "pl":
			    self.MARSCode['levelist'] = str(self.level)

	    ##retrieve the data
	    server.retrieve(self.MARSCode)



class MovieData:
    """
        Holds the data required to construct a movie using latitude-longitude data and animated in time
    """
    def __init__(self):

        # define empty containers 
        self.variable = np.array([])
        self.variableUnits = ""
        self.time = np.array([])
        self.timeUnits = np.array([])
        self.latitude = np.array([])      
        self.longitude = np.array([])
        self.level = np.array([])
        self.colourRange = np.array([])
        self.contours = np.array([])

        #define defaults
        self.colourmap = plt.cm.gist_ncar

    def defineFilename(self, query, path, extension):
        """
            Build the filename of the movie file
        """
        self.filename = buildEIFilename(query, path, ".mov")



    def setContours(self, newContours):
        """
            Defines a new colormap, newMapName must be the name of a matplotlib colormap. 
            All colormaps can be seen here: http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
        """
        self.contours = newContours


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
        self.colourRange = np.array([colourbarMinimum, colourbarMaximum])

        #re-define contours as well
        self.setContours(np.linspace(self.colourRange[0], self.colourRange[1], 20))
    

##############################################################################
out = plotEIMovie()


