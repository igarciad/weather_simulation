# Fast Weather Simulation for Inverse Procedural Design of 3D Urban Models
This is the source code for our paper **Fast Weather Simulation for Inverse Procedural Design of 3D Urban Models** by Ignacio Garcia-Dorado, Daniel Aliaga, Saiprasanth Bhalachandran, Paul Schmid, and Dev Niyogi, currently *accepted with minor revision* in ACM Transactions on Graphics TOG (10/2016).

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

## Instalation

The only requirement is GCC/G++. In Ubuntu, the default installed version (including live usb) should work. In Windows, we have used Visual Studio.

## Compilation (Ubuntu)

To run, just compule all cpp files with C++11 flags (this is not a hard constraint), and run the code.

    $ g++ *.cpp *.h -std=c++11 -o weather
    $ chmod +x weather
    $ ./weather
    
## Compilation (Windows)

We have tested using Visual Studio 2010 and 2013. Just create a project that include all files and run.

## Usage

The default run demo:

1. Load the initial weather state from a sounding. Note: We include code to load two types of sounding files.
2. Load the ground variebles. Here we randomly assign to each bottommost gridcell a land use, and then, we populate the ground variables from the provided files (Supplemental Material B). Note: In our implementation we load this from images but here it is random to avoid to have additional dependencies.
3. All the variables are initialized.
4. The system runs 1000 iterations.

## License

Code of the paper: Fast  Weather  Simulation Coupled with Inverse Procedural Design of 3D Urban Models
This project is licensed under the [MIT license](https://opensource.org/licenses/MIT)

Copyright (c) 2011-2016 Ignacio Garcia-Dorado, Daniel Aliaga.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

