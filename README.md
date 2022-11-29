# AstroPhysics-Projects

A series of projects from an AstroPhysics lab in 2020, where multiple sets of telescope data were analyzed. 

Project list: 
* Project 1: Assessing the accuracy of the modern detector's ability to count photons, which helps get a gauge for the error we will have when taking actual measurements with the telescope. Compared the error distributions to theoretical Gaussian and Poisson distributions, and computed the read noise and gain. 

* Project 2: Took data from a spectrometer and found the centroids (peaks) of our observed spectra using a centroid finding function that we wrote, compared them to theoretical wavelengths, estimating the error, and used this wavelength calibration to produce a final spectra. This could be used in further experiments to identify the main elemental components of different objects in space (such as galaxies). Primarily used image recognition functions that we wrote and linear least squares fit.

* Project 3: Analyzed two images taken by Nickel 1-m telescope at Lick Observatory to characterize the motion of the asteroid Parthenope (large, bright, main belt asteroid). Removing bias from the data was critical so we could successfully use image recognition and processing functions to identify the position of the asteroid. We then queried the USNO database to match up our asteroid centroids with the existing data. Data correction also had to be done to account for the fact that we could not assume an ideal camera, and we also needed to correct for aspects such as shear and distortion. 

Primary packages used: numpy, pandas, astropy.io.fits, matplotlib. Most functions (even primitive ones) were written in house for practice. 

Note: The code offers little explanation, looking at the report pdfs will allows you to interpret what is going on much better. 
