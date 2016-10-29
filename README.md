# Fast Weather Simulation for Inverse Procedural Design of 3D Urban Models
This is the source code for our paper **Fast Weather Simulation for Inverse Procedural Design of 3D Urban Models**, currently *accepted with minor revision* in ACM Transactions on Graphics TOG (10/2016).

## Abstract 
We present the first realistic, physically-based, fully coupled, real-time
weather design tool for use in urban procedural modeling. We merge de-
signing of a 3D urban model with a controlled long-lasting spatiotemporal
interactive simulation of weather. Starting from the fundamental dynami-
cal equations similar to those used in state-of-the-art weather models, we
present a novel simplified urban weather model for interactive graphics.
Control of physically-based weather phenomena is accomplished via an in-
verse modeling methodology. In our results, we present several scenarios of
forward design, inverse design with high-level and detailed-level weather
control and optimization, as well as comparisons of our method against
well-known weather simulation results and systems.

## Platforms

Linux (tested in Ubuntu) and Windows (tested on Windows 10) are currently supported but it should work with any platform that compiles C++.

## Install

The only requirement is GCC/G++. In Ubuntu, the default installed version (including live usb) should work.

## Run (Ubuntu)

To run, just compule all cpp files with C++11 flags (this is not a hard constraint), and run the code.

    $ g++ *.cpp *.h -std=c++11 -o weather
    $ chmod +x weather
    $ ./weather
    
## Run (Windows)

We have tested using Visual Studio 2010 and 2013. Just create a project that include all files and run.

## Result

The default run demo:

1. Load the initial weather state from a sounding. Note: We include code to load two types of sounding files.
2. Load the ground variebles. Here we randomly assign to each bottommost gridcell a land use, and then, we populate the ground variables from the provided files (Supplemental Material B). Note: In our implementation we load this from images but here it is random to avoid to have additional dependencies.
3. All the variables are initialized.
4. The system runs 1000 iterations.
