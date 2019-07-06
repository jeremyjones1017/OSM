# OSM - Oblate Star Model
Models oblate stars to determine their ages (among other parameters)

Generates a model star with the following parameters (referred to as 'free parameters' below):
-Equatorial Radius
-Equartorial Rotation Velocity
-Inclination
-Polar Temperature
-Position Angle

Synthetic interferometric and photometric data are created for this model star given the following fixed parameters:
-Distance
-Stellar Mass
-Gravity darkening law

We use a markov chain monte carlo (mcmc) simulation to generate the breadth of free parameters that fit the observed data for the given star.

We use evolutionary models to determine the age, mass, and initial rotation velocity for each point in the distribution.


Original project (rotmodel) began ~2014

This updated version began ~2016
