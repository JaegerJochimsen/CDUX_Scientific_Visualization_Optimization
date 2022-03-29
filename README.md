# Optimizations for Post-Hoc Interpolation in Flow Visualization Using Delaunay Triangulations
## Background:
Consider a set of basis flows S, where a flow is denoted by a start and an end point,
S = {(s<sub>i</sub>, e<sub>i</sub>) | 0 <= i <= |S|}.

Our <strong>goal</strong> is to interpolate end positions for a set of query points Q, where 
Q = {(q<sub>i</sub>, ?) | 0 <= i <= |Q|}.

All points containted within S and Q, are scattered over a cube-shaped domain, and that domain
is split into 64 subdomains (each of uniform size and shape).

## Dependencies
1. CGAL
2. VTK

## File Manifest
*Software Files*
1. reconstruct.cxx      - source code for the project
2. genPts.py            - a script used for generating .vtk test data of variable length 

*Data File(s)*
1. Groundtruth_0_0.vtk  - a set of test start and end points for basis flows

*Config File(s)*
1. CMakeLists.txt       - a cmake file used to generate a makefile for the project
## Known Bugs and Fixes 
1. 
## Credits
1. Hank Childs, PhD
2. Sudhanshu Sane, PhD
## What's Next
1.  
