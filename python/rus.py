# Author: Jerome H.L. Le Rousseau (jerome@dix.mines.edu)
# Center for Wave Phenomena / Physical Acoustic Laboratory
# Colorado School of Mines
# Golden, CO 80401 USA

# Updated 2014
# By: Leighton Watson (lwat054@aucklanduni.ac.nz)
# Physical Acoustic Laboratory
# University of Auckland, Auckland, 1010, New Zealand

# Translated to Python 2015
# By: Paul Freeman (pfre484@aucklanduni.ac.nz)
# Computer Science Department
# University of Auckland, Auckland, 1010, New Zealand

import rus_parser as parser
import rus_forward as forward
import rus_inverse as inverse

if __name__ == "__main__":

    args = parser.start()

    if args.subcommand == 'forward':
        forward.start(args)

    if args.subcommand == 'inverse':
        inverse.start(args)

