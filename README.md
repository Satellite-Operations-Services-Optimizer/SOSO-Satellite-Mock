# SOSO-Satellite-Mock
Repo for mock satellites.
mock_SAT_ver2.py - Main file.

## Prerequisites - Libraries


## Summary:
This repository contains the code used for mapping of the satellites used in the database, scheduler microservice, and more. It is meant to be used as a template for a satellite and can be expanded on by including more theories of satellite engineering. The satellites are defined through TLE files in a .json format (see below).

### Example of .json format for TLE file:
code()

## Current Functionality:
1. Read TLE files in json.

2. Apply ephemeris data for the Sun and Earth.
3. Power Management during eclipse or sunlight condition (1000 W - sunlight, 400 W - eclipse).
4. Time Interval.
5. Positional calculations for the satellite (latitude, longitude, x, y, z).
6. Image order validation using propagated FOV.
7. 2D plot visualization.

## References:

- Power Management:
- FOV Ground Track:

