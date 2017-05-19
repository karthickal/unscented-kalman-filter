# Unscented Kalman Filter

This project is part of Udacity's Self-Driving Car Engineer Nanodegree Program. For more information [read here](https://in.udacity.com/course/self-driving-car-engineer-nanodegree--nd013/)

## Overview

This repository contains a c++ implementation of an unscented Kalman Filter that predicts state from RADAR and LIDAR data. The prediction step calculates the state and covariance by predicting sigma points for the non-linear CRTV model. The update step calculates the posterior for RADAR and LIDAR separately using the cross correlation matrix.

For the sample data the calculated RMSE is - 
* 0.0720901
* 0.0894626
* 0.381186
* 0.261678
## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`

