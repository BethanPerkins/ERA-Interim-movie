#
#
# This code will create a movie of ERA Interim data through python.
#
#
# Built by bethan.perkins@assimila.eu May 2014
#
#
# Help must provide a list of available variables

# hard-coded definition: the directory tha the data is stored in. Global variable


#define the user's required data block
#inputs - variable, startdate, enddate, resolution (we're sticking to the same for lat and lon), lat-lon extents, level
#output - an object which defines all these things 

#build a filename from the user-block
# inputs - user-block
# outputs - filename

#check to see if this already exists
# inputs - filename
# outputs - true/false, 

# --> called only if the output to previous is false
# download it if it doesn't exist
# inputs - user-block object, filename
# pre-defined: a database of the variable numbers and other preripheral information required in order to download the data
# output - none


# load data extract data and captions etc from the NetCDF file
# inputs - filename
# process: open file, define new object, populate new object with information from the file.
# outputs - object with the required data (lat, long, time), dates (time), lat extent (lat), lon extent (lon), Name (long_name), units (units)


#build a simple movie from it 



#save that movie somewhere useful

