# cython: profile=True

import rus_parser as parser
import rus_forward as forward
import rus_inverse as inverse
import profile
import cProfile

print ("accessed outside")
#if __name__ == "__main__":
args = parser.start()
print ("Accessed INSIDE")
if args.subcommand == 'forward':
    forward.start(args)
    print ("Accessed 1")

if args.subcommand == 'inverse':
    inverse.start(args)
    print ("Accessed 2")