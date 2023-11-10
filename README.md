# SOSO-Satellite-Mock
Repo for mock satellite.
mock_SAT_ver2.py - Main file.

## Current Functionality:
1. Read TLE files in json.
2. Apply ephemeris data for the Sun and Earth.
3. Power Management during eclipse or sunlight condition (1000 W - sunlight, 400 W - eclipse).
4. Time Interval.
5. Positional calculations for the satellite (latitude, longitude, x, y, z).
6. Image order validation using propagated FOV.
7. 2D plot visualization.

## Next Steps:
1. Utilize FOV calculation to determine if image orders are acceptable or not.
2. Implement downlink rate from satellite to ground station.
