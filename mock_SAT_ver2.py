## Version 2: Contains:
# 5 SOSO satellites
# Positional calculations
# Time interval calculation
# Plotting for 2d mapping of satellites
# Conversion of positional data to json for future use in cesium
# Is the satellite in eclipse or in sunlight?
# (Additional_Ver1.1) Power restriction based on eclipse of sunlight
# (Additional_Ver2) FOV ground track to satellite calculation (References: Earth Coverage Paper) [In Progress]
# (Additional_Ver2) Image Order validation based on FOV ground track to satellite [To Do]
# (Additional_Ver2) Power degradation documentation (References: ) [To Do]

## Imports
from skyfield.api import load, EarthSatellite, Topos
from datetime import timedelta, datetime
import matplotlib.pyplot as plt
import json
import numpy as np
from math import atan, degrees

## Step 1: Load the TLE files
ts = load.timescale() # Create timescale object for TLE computation
satellites = [] # Empty list of satellite objects (SOSO-1, SOSO-2, etc)

# Create empty lists for the latitudes and longitudes of each satellite
latitudes = [[] for _ in range(5)]
longitudes = [[] for _ in range(5)]

for i in range(1, 6): # Iterate over numbers 1 to 5
    try: # Used to catch and handle exceptions in the code
        with open(f'SOSO-{i}_TLE.json') as f: # For-loop and f-string used to open the TLE files for SOSO-1, SOSO-2, etc.
            data = json.load(f) # Load the JSON data from the file
            name = data['name']
            line1 = data['line1']
            line2 = data['line2']
        satellite = EarthSatellite(line1, line2, name, ts) # Create new satellite object where line 1 = tle[1], line 2 = tle[2], title = tle[0], and ts for timescale
        satellites.append(satellite) # Add each satellite to the empty list
    except IndexError: # Handles TLE files without title line or missing lines
        print(f"Error: TLE file for satellite {i} is not formatted correctly.") # Output of error

## Step 4: (Maintenance) Is the satellite in eclipse or in sunlight?
eph = load('de421.bsp')  # Load the JPL ephemeris DE421
sun = eph['sun']  # Get the 'sun' object from the ephemeris
earth = eph['earth']  # Get the 'earth' object from the ephemeris

# Power Management Example
P_sunlit = 1000 # in Watts during Sunlight
P_eclipse = P_sunlit * 0.4 # in Watts during Eclipse (assuming 40% of power is used)

## Step 2: Select time interval for satellite and ground station accesses.

# Ask the user for the start and end times
start_time_str = input("Enter the start time (YYYY-MM-DD HH:MM:SS): ")
end_time_str = input("Enter the end time (YYYY-MM-DD HH:MM:SS): ")

# Convert the input strings to datetime objects
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S")

# Convert the datetime objects to skyfield Time objects
start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

## Step 3: Get live data of the satellites' position (i.e. x, y, z coordinates, latitude, longitude)
for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
    current_time_skyfield = start_time_skyfield
    while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
        position_current = satellite.at(current_time_skyfield).position.km  # Plain (x, y, z) coordinates at the current time
        subpoint_current = satellite.at(current_time_skyfield).subpoint()
        latitude_current = subpoint_current.latitude.degrees  # Latitude at the current time
        longitude_current = subpoint_current.longitude.degrees  # Longitude at the current time
        
        # Get the positions of the Earth, Sun, and satellite
        earth_pos = earth.at(current_time_skyfield).position.km
        sun_pos = sun.at(current_time_skyfield).position.km
        satellite_pos = (earth + satellite).at(current_time_skyfield).position.km

        # Calculate the vectors from the satellite to the Earth and Sun
        satellite_to_earth = earth_pos - satellite_pos
        satellite_to_sun = sun_pos - satellite_pos

        # Calculate the angle between these vectors
        angle = satellite_to_earth.angle(satellite_to_sun).degrees

        # Append the latitude and longitude to their respective lists
        latitudes[i].append(latitude_current)
        longitudes[i].append(longitude_current)

        # Calculate altitude from position data
        altitude_current = np.linalg.norm(position_current) - 6371  # Earth's radius is approximately 6371 km

        # Calculate FOV
        fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))

        print(f"SOSO-{i + 1} at Current Time:")
        print(f"  Position: {position_current}")
        print(f"  Latitude: {latitude_current}")
        print(f"  Longitude: {longitude_current}")
        print(f"  Altitude: {altitude_current} km")
        print(f"  Field of View: {fov_current} degrees")
        print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

        current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1))

for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
    current_time_skyfield = start_time_skyfield
    while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
        sat_pos_current = (earth + satellite).at(current_time_skyfield).position.km  # Get the position of the satellite relative to Earth at the current time

        is_sunlit_current = satellite.at(current_time_skyfield).is_sunlit(eph) # Check if satellite is sunlit at current time

        if is_sunlit_current:
            print(f"SOSO-{i + 1} at Current Time: {is_sunlit_current}")
            print(f"  The satellite is in sunlight. Power is unrestricted = {P_sunlit} W.")
            print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")
        else:
            print(f"SOSO-{i + 1} at Current Time: {is_sunlit_current}")
            print(f"  The satellite is in eclipse. All activities must use less power than a predefined limit = {P_eclipse} W.")
            print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

        current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1))

###################################################################################################################################################################################################
## FRONT-END
###################################################################################################################################################################################################

# Convert lists to numpy arrays after collecting all data
latitudes_soso1 = np.array(latitudes[0])
longitudes_soso1 = np.array(longitudes[0])

latitudes_soso2 = np.array(latitudes[1])
longitudes_soso2 = np.array(longitudes[1])

latitudes_soso3 = np.array(latitudes[2])
longitudes_soso3 = np.array(longitudes[2])

latitudes_soso4 = np.array(latitudes[3])
longitudes_soso4 = np.array(longitudes[3])

latitudes_soso5 = np.array(latitudes[4])
longitudes_soso5 = np.array(longitudes[4])

## [FRONT-END] Cesium Front-end (Satellite View) -> Parsing data into json file for js (cesium_SAT.js)
# Create a dictionary to hold your data
data = {
    'soso1': {
        'latitudes': latitudes_soso1.tolist(),
        'longitudes': longitudes_soso1.tolist()
    },
    'soso2': {
        'latitudes': latitudes_soso2.tolist(),
        'longitudes': longitudes_soso2.tolist()
    },
    'soso3': {
        'latitudes': latitudes_soso3.tolist(),
        'longitudes': longitudes_soso3.tolist()
    },
    'soso4': {
        'latitudes': latitudes_soso4.tolist(),
        'longitudes': longitudes_soso4.tolist()
    },
    'soso5': {
        'latitudes': latitudes_soso5.tolist(),
        'longitudes': longitudes_soso5.tolist()
    },
    # ... continue this for all your satellites ...
}

# Write the data to a JSON file
with open('satellite_data.json', 'w') as f:
    json.dump(data, f)
    
# Create a new figure
fig = plt.figure()

# Add a subplot with a projection of 'mollweide'
ax = fig.add_subplot(111, projection='mollweide')

# Convert degrees to radians as required by matplotlib for mollweide projection
latitudes_soso1_rad = np.deg2rad(latitudes_soso1)
longitudes_soso1_rad = np.deg2rad(longitudes_soso1)

latitudes_soso2_rad = np.deg2rad(latitudes_soso2)
longitudes_soso2_rad = np.deg2rad(longitudes_soso2)

latitudes_soso3_rad = np.deg2rad(latitudes_soso3)
longitudes_soso3_rad = np.deg2rad(longitudes_soso3)

latitudes_soso4_rad = np.deg2rad(latitudes_soso4)
longitudes_soso4_rad = np.deg2rad(longitudes_soso4)

latitudes_soso5_rad = np.deg2rad(latitudes_soso5)
longitudes_soso5_rad = np.deg2rad(longitudes_soso5)

# Plot the trajectories of the satellites
ax.plot(longitudes_soso1_rad, latitudes_soso1_rad, color='red', label='SOSO-1')
ax.plot(longitudes_soso2_rad, latitudes_soso2_rad, color='yellow', label='SOSO-2')
ax.plot(longitudes_soso3_rad, latitudes_soso3_rad, color='pink', label='SOSO-3')
ax.plot(longitudes_soso4_rad, latitudes_soso4_rad, color='blue', label='SOSO-4')
ax.plot(longitudes_soso5_rad, latitudes_soso5_rad, color='green', label='SOSO-5')

# Add a legend
ax.legend()

# Show the plot
plt.show()
