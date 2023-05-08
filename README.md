# SessileDropContactAngle

This repository contains a Python script that calculates the contact angle of a sessile droplet from molecular dynamics simulations.
Description

The script uses the MDAnalysis library to read in the simulation trajectory and topology files. It then selects the water atoms and calculates the z-density and central density of the water molecules in the sessile droplet. The script calculates contact angles with linear extrapolation and spline fitting. It outputs the mean linear and spline contact angles and saves these contact angles into .xvg files.

## Requirements
 - Python (3.7 or above)
 - MDAnalysis
 - numpy
 - scipy

## Installation

Clone this repository to your local machine using:
```
git clone https://github.com/TobiasMaterzok/SessileDropContactAngle.git
```

## Input

The script takes command-line arguments which include the names of the trajectory and topology files, the output file name, the prefix for the output file, the size of the z-bin, and the start and end steps for the analysis.

## Output

The script outputs two .xvg files: one for the contact angles calculated with linear curve fitting and one for the contact angles calculated with spline fitting. Each file contains the timestep in picoseconds and the corresponding contact angle.

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.
