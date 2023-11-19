## Version 2: Contains:
# 5 SOSO satellites
# Positional calculations
# Time interval calculation
# Plotting for 2d mapping of satellites
# Conversion of positional data to json for future use in cesium
# Is the satellite in eclipse or in sunlight?
# (Additional_Ver1.1) Power restriction based on eclipse of sunlight
# (Additional_Ver2) FOV ground track to satellite calculation (References: Earth Coverage Paper) [In Progress]
# (Additional_Ver2) Power degradation documentation (References: ) [To Do]

## Imports
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from math import atan, degrees
from skyfield.api import load, EarthSatellite, Topos
from datetime import timedelta, datetime
from skyfield.sgp4lib import EarthSatellite

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

## Step 2: (Maintenance) Is the satellite in eclipse or in sunlight?
eph = load('de421.bsp')  # Load the JPL ephemeris DE421
sun = eph['sun']  # Get the 'sun' object from the ephemeris
earth = eph['earth']  # Get the 'earth' object from the ephemeris

# Power Management Example
P_sunlit = 500 # in Watts during Sunlight
# 200-800 Watts for research sat.
# 1000-1500 Watts for commercial sat.
P_eclipse = P_sunlit * 0.4 # in Watts during Eclipse (assuming 40% of power is used)

## Step 3: Select time interval for satellite and ground station accesses.

# Ask the user for the start and end times
start_time_str = input("Enter the start time (YYYY-MM-DD HH:MM:SS): ")
end_time_str = input("Enter the end time (YYYY-MM-DD HH:MM:SS): ")

# Convert the input strings to datetime objects
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S")

# Convert the datetime objects to skyfield Time objects
start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

## Step 4: Get live data of the satellites' position (i.e. x, y, z coordinates, latitude, longitude)
for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
    current_time_skyfield = start_time_skyfield
    while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
        position_current = satellite.at(current_time_skyfield).position.km  # Plain (x, y, z) coordinates at the current time (Center of Earth)
        subpoint_current = satellite.at(current_time_skyfield).subpoint()
        latitude_current = subpoint_current.latitude.degrees  # Latitude at the current time
        longitude_current = subpoint_current.longitude.degrees  # Longitude at the current time
        
        # Get the positions of the Earth, Sun, and satellite
        earth_pos = earth.at(current_time_skyfield).position.km
        sun_pos = sun.at(current_time_skyfield).position.km

        # Append the latitude and longitude to their respective lists
        latitudes[i].append(latitude_current)
        longitudes[i].append(longitude_current)

        # Calculate altitude from position data
        semi_major_axis_km = satellite.model.a * 6378.137  # Get the semi-major axis in kilometers
        altitude_current = semi_major_axis_km - 6378.137  # The altitude is the semi-major axis minus the Earth's radius

        # Calculate FOV
        # fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))
        # No need to calculate FOV since it is given to us.

        print(f"SOSO-{i + 1} at Current Time:")
        print(f"  Position: {position_current} km")
        print(f"  Latitude: {latitude_current} degrees")
        print(f"  Longitude: {longitude_current} degrees")
        print(f"  Altitude: {altitude_current} km")
        # print(f"  Field of View: {fov_current} degrees")
        print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

        current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1)) # Print all variables every minute from start and end times.

for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
    current_time_skyfield = start_time_skyfield
    while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
        sat_pos_current = satellite.at(current_time_skyfield).position.km  # Get the position of the satellite relative to Earth at the current time

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

## Step 5: Load the Image Orders & Validation Methodology

# Spotlight image = 10 km height x 10 km width square area, 120 second transfer time, 512 MB.
# Medium image = 40 km height x 20 km width square area, 45 second transfer time, 256 MB.
# Low image = 40 km height x 20 km width square area, 20 second transfer time, 128 MB.

# Given: Satellite FOV -> Viewing Angle = 30 degrees, Full View = 60 degrees.
# If FOV square area of 577 km height x 577 km width propagated on ground track doesn't cover square area of image in image order list, then unacceptable.
# Variables to consider in image order: Latitude, Longitude, Image Start Time, Image End Time, Revisit Time (True or False).
# Image order square area should be within the satellite FOV square area throughout the transfer time during between the start and end time of the image.

# False Image Orders: Order 37 for, 

# # Define the image types and their dimensions
# image_types = {
#     "Spotlight": [10, 10],
#     "Medium": [40, 20],
#     "Low": [40, 20]
# }

# # Iterate over the image orders
# for i in range(1, 51):
#     with open(f'SampleOrders/Order_{i}.json') as f:
#         order = json.load(f)
#     image_type = order["ImageType"]
#     lat, lon = order["Latitude"], order["Longitude"]
#     start_time = datetime.strptime(order["ImageStartTime"], "%Y-%m-%dT%H:%M:%S")
#     end_time = datetime.strptime(order["ImageEndTime"], "%Y-%m-%dT%H:%M:%S")

#     # Calculate the square area for the image type
#     image_area = image_types[image_type]

#     # Iterate over the satellites
#     for satellite in satellites:
#         # Calculate the satellite's FOV square area
#         satellite_area = [577, 577]  # Replace with actual calculation

#         # Note 2: Dependent on altitude. 
        
#         # Check if the image area falls within the satellite's FOV
#         if image_area[0] <= satellite_area[0] and image_area[1] <= satellite_area[1]:
#             # Calculate the transfer time
#             transfer_time = end_time - start_time
#             # Note 1: Change for 120, 45, 20 sec for transfer time.

#             # Check if the satellite stays within the FOV during the transfer time
#             t = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
#             geocentric = satellite.at(t)
#             subpoint = geocentric.subpoint()
#             start_lat, start_lon = subpoint.latitude.degrees, subpoint.longitude.degrees

#             t = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)
#             geocentric = satellite.at(t)
#             subpoint = geocentric.subpoint()
#             end_lat, end_lon = subpoint.latitude.degrees, subpoint.longitude.degrees

#             # If the satellite stays within the FOV, the image is acceptable
#             if start_lat <= lat <= end_lat and start_lon <= lon <= end_lon:
#                 print(f"Image order {i} is unacceptable for satellite {satellite.name}")
#             else:
#                 print(f"Image order {i} is acceptable for satellite {satellite.name}")
#         else:
#             print(f"Image order {i} is acceptable for satellite {satellite.name}")

#         # Plot the FOV and image area
#         fig, ax = plt.subplots()
#         ax.add_patch(plt.Rectangle((lon - image_area[1] / 2, lat - image_area[0] / 2), image_area[1], image_area[0], fill=None, edgecolor='r'))
#         ax.add_patch(plt.Rectangle((start_lon - satellite_area[1] / 2, start_lat - satellite_area[0] / 2), satellite_area[1], satellite_area[0], fill=None, edgecolor='b'))
#         ax.set_xlim([min(lon - image_area[1] / 2, start_lon - satellite_area[1] / 2) - 10, max(lon + image_area[1] / 2, start_lon + satellite_area[1] / 2) + 10])
#         ax.set_ylim([min(lat - image_area[0] / 2, start_lat - satellite_area[0] / 2) - 10, max(lat + image_area[0] / 2, start_lat + satellite_area[0] / 2) + 10])
#         plt.show()
        
#         # Note 3: add the lat and lon for the time for when image is acceptable.

# Define image types and their properties
image_types = {
    "Spotlight": {"height": 10, "width": 10, "transfer_time": 120},
    "Medium": {"height": 40, "width": 20, "transfer_time": 45},
    "Low": {"height": 40, "width": 20, "transfer_time": 20}
}

# Define the viewing angle
viewing_angle = 30  # degrees

# Iterate over the image orders
for i in range(1, 50):
    with open(f'SampleOrders/Order_{i}.json') as f:
        order = json.load(f)

    # Get the image type
    image_type = image_types[order["ImageType"]]

    # Calculate the FOV square area
    for satellite in satellites:
        semi_major_axis_km = satellite.model.a * 6378.137  # Get the semi-major axis in kilometers
        altitude_current = semi_major_axis_km - 6378.137  # The altitude is the semi-major axis minus the Earth's radius
        fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))  # Calculate FOV

        # Check if the image order falls within the FOV
        if image_type["height"] <= fov_current and image_type["width"] <= fov_current:
            # Calculate the start and end times
            start_time = datetime.strptime(order["ImageStartTime"], "%Y-%m-%dT%H:%M:%S")
            end_time = datetime.strptime(order["ImageEndTime"], "%Y-%m-%dT%H:%M:%S")
            
            # Convert the datetime objects to skyfield Time objects
            start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
            end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

            # Get the satellite's position at the start and end times
            start_position = satellite.at(start_time_skyfield).position.km
            end_position = satellite.at(end_time_skyfield).position.km

            # Check if the satellite stays within the FOV during the transfer time
            transfer_end_time = start_time + timedelta(seconds=image_type["transfer_time"])
            if start_time <= transfer_end_time <= end_time:
                if start_position[0] <= end_position[0] and start_position[1] <= end_position[1]:
                    print(f"Order {i} is acceptable for {satellite.name}.")
                    print(f"Transfer completed at time {transfer_end_time}, lat {order['Latitude']}, lon {order['Longitude']}")
                    
                    # Plot the satellite's position
                    plt.figure(figsize=(10, 10))
                    plt.plot([order['Longitude']], [order['Latitude']], 'ro')
                    plt.xlim(order['Longitude'] - 1, order['Longitude'] + 1)
                    plt.ylim(order['Latitude'] - 1, order['Latitude'] + 1)
                    plt.grid(True)
                    plt.title(f"Satellite {satellite.name} Position for Order {i}")
                    plt.xlabel("Longitude")
                    plt.ylabel("Latitude")
                    plt.show()
                else:
                    print(f"Order {i} is unacceptable for {satellite.name}.")
            else:
                print(f"Order {i} is unacceptable for {satellite.name}.")
        else:
            print(f"Order {i} is unacceptable for {satellite.name}.")

###################################################################################################################################################################################################
################################ FUNCTIONS
###################################################################################################################################################################################################

# Function 1: Read TLE files for SOSO constellation
def load_satellites():
    ts = load.timescale()
    satellites = []

    for i in range(1, 6):
        try:
            with open(f'SOSO-{i}_TLE.json') as f:
                data = json.load(f)
                name = data['name']
                line1 = data['line1']
                line2 = data['line2']
            satellite = EarthSatellite(line1, line2, name, ts)
            satellites.append(satellite)
        except IndexError:
            print(f"Error: TLE file for satellite {i} is not formatted correctly.")

    # return ts, satellites
    pass

# Function 2: Ephemeris Data, Power Draw, and User-defined Start and End Times for Schedule
def ephemeris_power_schedule_start_end():
    ## Step 2: (Maintenance) Is the satellite in eclipse or in sunlight?
    eph = load('de421.bsp')  # Load the JPL ephemeris DE421
    sun = eph['sun']  # Get the 'sun' object from the ephemeris
    earth = eph['earth']  # Get the 'earth' object from the ephemeris

    # Power Management Example
    P_sunlit = 500 # in Watts during Sunlight
    # 200-800 Watts for research sat.
    # 1000-1500 Watts for commercial sat.
    P_eclipse = P_sunlit * 0.4 # in Watts during Eclipse (assuming 40% of power is used)

    ## Step 3: Select time interval for satellite and ground station accesses.

    # Ask the user for the start and end times
    start_time_str = input("Enter the start time (YYYY-MM-DD HH:MM:SS): ")
    end_time_str = input("Enter the end time (YYYY-MM-DD HH:MM:SS): ")

    # Convert the input strings to datetime objects
    start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S")

    # Convert the datetime objects to skyfield Time objects
    start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
    end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)
    pass

# Function 3: Location of the Satellites
def location_satellite():
    ## Step 4: Get live data of the satellites' position (i.e. x, y, z coordinates, latitude, longitude)
    for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
        current_time_skyfield = start_time_skyfield
        while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
            position_current = satellite.at(current_time_skyfield).position.km  # Plain (x, y, z) coordinates at the current time (Center of Earth)
            subpoint_current = satellite.at(current_time_skyfield).subpoint()
            latitude_current = subpoint_current.latitude.degrees  # Latitude at the current time
            longitude_current = subpoint_current.longitude.degrees  # Longitude at the current time
            
            # Get the positions of the Earth, Sun, and satellite
            earth_pos = earth.at(current_time_skyfield).position.km
            sun_pos = sun.at(current_time_skyfield).position.km

            # Append the latitude and longitude to their respective lists
            latitudes[i].append(latitude_current)
            longitudes[i].append(longitude_current)

            # Calculate altitude from position data
            semi_major_axis_km = satellite.model.a * 6378.137  # Get the semi-major axis in kilometers
            altitude_current = semi_major_axis_km - 6378.137  # The altitude is the semi-major axis minus the Earth's radius

            # Calculate FOV
            # fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))
            # No need to calculate FOV since it is given to us.

            print(f"SOSO-{i + 1} at Current Time:")
            print(f"  Position: {position_current} km")
            print(f"  Latitude: {latitude_current} degrees")
            print(f"  Longitude: {longitude_current} degrees")
            print(f"  Altitude: {altitude_current} km")
            # print(f"  Field of View: {fov_current} degrees")
            print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

            current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1)) # Print all variables every minute from start and end times.
    pass

# Function 4: Power Management
def eclipse_sunlight():
    for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
        current_time_skyfield = start_time_skyfield
        while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
            sat_pos_current = satellite.at(current_time_skyfield).position.km  # Get the position of the satellite relative to Earth at the current time

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
    pass

# Function 5: Image Validation
def image_validation():
    # Define image types and their properties
    image_types = {
        "Spotlight": {"height": 10, "width": 10, "transfer_time": 120},
        "Medium": {"height": 40, "width": 20, "transfer_time": 45},
        "Low": {"height": 40, "width": 20, "transfer_time": 20}
    }

    # Define the viewing angle
    viewing_angle = 30  # degrees

    # Iterate over the image orders
    for i in range(1, 50):
        with open(f'SampleOrders/Order_{i}.json') as f:
            order = json.load(f)

        # Get the image type
        image_type = image_types[order["ImageType"]]

        # Calculate the FOV square area
        for satellite in satellites:
            semi_major_axis_km = satellite.model.a * 6378.137  # Get the semi-major axis in kilometers
            altitude_current = semi_major_axis_km - 6378.137  # The altitude is the semi-major axis minus the Earth's radius
            fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))  # Calculate FOV

            # Check if the image order falls within the FOV
            if image_type["height"] <= fov_current and image_type["width"] <= fov_current:
                # Calculate the start and end times
                start_time = datetime.strptime(order["ImageStartTime"], "%Y-%m-%dT%H:%M:%S")
                end_time = datetime.strptime(order["ImageEndTime"], "%Y-%m-%dT%H:%M:%S")
                
                # Convert the datetime objects to skyfield Time objects
                start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
                end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

                # Get the satellite's position at the start and end times
                start_position = satellite.at(start_time_skyfield).position.km
                end_position = satellite.at(end_time_skyfield).position.km

                # Check if the satellite stays within the FOV during the transfer time
                transfer_end_time = start_time + timedelta(seconds=image_type["transfer_time"])
                if start_time <= transfer_end_time <= end_time:
                    if start_position[0] <= end_position[0] and start_position[1] <= end_position[1]:
                        print(f"Order {i} is acceptable for {satellite.name}.")
                        print(f"Transfer completed at time {transfer_end_time}, lat {order['Latitude']}, lon {order['Longitude']}")
                        
                        # Plot the satellite's position
                        plt.figure(figsize=(10, 10))
                        plt.plot([order['Longitude']], [order['Latitude']], 'ro')
                        plt.xlim(order['Longitude'] - 1, order['Longitude'] + 1)
                        plt.ylim(order['Latitude'] - 1, order['Latitude'] + 1)
                        plt.grid(True)
                        plt.title(f"Satellite {satellite.name} Position for Order {i}")
                        plt.xlabel("Longitude")
                        plt.ylabel("Latitude")
                        plt.show()
                    else:
                        print(f"Order {i} is unacceptable for {satellite.name}.")
                else:
                    print(f"Order {i} is unacceptable for {satellite.name}.")
            else:
                print(f"Order {i} is unacceptable for {satellite.name}.")
    pass

###################################################################################################################################################################################################
################################ FRONT-END
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
fig = plt.figure(1)

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

plt.figure(2)
plt.title('Top Down View of Satellite FOV and Image Order')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()