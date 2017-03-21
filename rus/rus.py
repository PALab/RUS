'''
Start a RUS experiment
'''
import rus.rus_parser as parser
import rus.rus_forward as forward
import rus.rus_inverse as inverse

def main():
    '''Start the forward or reverse code, depending on the subcommand.
    '''

    args = parser.start()

    if args.subcommand == 'forward':
        forward.start(args)

    if args.subcommand == 'inverse':
        inverse.start(args)

if __name__ == "__main__":

    main()

