# RUS
Resonant Ultrasound Spectroscopy

This project contains tools used in performing resonant ultrasound spectroscopy on rocks and other materials.

## C Implementation
### Setup
Detailed instructions can be found in the [Installation and Help Guide](Zadler%20Computer%20Code%20Manual.pdf).

### System Requirements
This implementation is designed for Linux and requires
installation of several other components. Please follow
the installation guide carefully or consider using the
Python implementation.

## Python Implementation
### Motivation
The C implementation depends on a large number of
libraries and installation can be difficult.
It was for this reason that a Python implementation
was desired.

### System Requirements
The Python implementation is designed for cross-platform
support and should run on Windows, MacOS, and Linux. It
should support both Python 2 and Python 3. Any platform
related issues with the Python implementation should be
reported as bugs.

This version depends on the SciPy library for its
linear algebra module. However, this module is free
and readily available from www.scipy.org.

### Performance
Although this Python implementation has been optimized
where possible, there is a noteable decrease in performance
from the original C version. This was anticipated.

If the Python version is unable to meet your speed
requirements, you should be able to an order of magnitude
increase in speed by using the C version.

## Credits
### Original C implementation by
Jerome H.L. Le Rousseau (jerome@dix.mines.edu)  
Center for Wave Phenomena / Physical Acoustic Laboratory  
Colorado School of Mines  
Golden, CO 80401 USA

### Updated in 2014 by
Leighton Watson (lwat054@stanford.edu)  
Physical Acoustic Laboratory  
University of Auckland, Auckland, 1010, New Zealand

### Translated to Python in 2015 by
Paul Freeman (pfre484@aucklanduni.ac.nz)  
Computer Science Department  
University of Auckland, Auckland, 1010, New Zealand

