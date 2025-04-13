####################################################################################################################################################################
# Importing Libraries
import numpy as np
####################################################################################################################################################################

# Function 	: getSystemConfig
# Details	: Fetches the configs that are necessary for the system.
def getSystemConfig():

	## Satellite Details
	# Orbital Velocity of the Satellite
	satelliteVelocity 			= 1000 #m/s
	# Satellite to Centre of Earth Distance
	distSatelliteEarthCentre	= 7000 #m
	# Satellite's 1st Location
	satelliteLoc1 				= []
	# Satellite's 2nd Location
	satelliteLoc2 				= []

	## UE Details
	# UE Location
	ueLocation 					= []
	# UE Azimuth wrt. North
	ueAzimuthNorth