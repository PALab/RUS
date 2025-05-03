# cython: profile=True

import rus_parser as parser
import rus_forward as forward
import rus_inverse as inverse
import profile
import cProfile

#if __name__ == "__main__":
args = parser.start()
if args.subcommand == 'forward':
    forward.start(args)

if args.subcommand == 'inverse':
    inverse.start(args)