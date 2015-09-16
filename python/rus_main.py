import sys
import rus
import rus_parser
import numpy

def start(argv):

    # parser provided arguments
    args = rus_parser.inverse_parser(argv)

    # half sample dimensions are used in calculations
    dimension1 = args.d1 / 2.0
    dimension2 = args.d2 / 2.0
    dimension3 = args.d3 / 2.0
    dimensions = [dimension1, dimension2, dimension3]

    # dimension of the problem
    a = args.order + 1
    b = args.order + 2
    c = args.order + 3
    problem_size = 3 * a * b * c / 6

    # relationship between ir and l,m,n - filling tables
    tabs, irk = rus.index_relationship(args.order, problem_size)

    # read the input file
    # remove frequencies above or below the max or min
    freq_list = read_input(args.input_file)
    y = [(i,j) for (i,j) in freq_list if i > args.freqmin and i < args.freqmax]

    # calculate number of frequencies used for inversion
    num_freqs = len(y)

    print_args(args, num_freqs)

    ia  = [1 for i in range(num_freqs)]

        # sig is not a true 'error' but the weighting from freq_data
        # is multipled by the difference in the misfit function/chisq so modes with a zero weighting are ignored


    covar  = numpy.identity(args.ns)
    alpha  = numpy.identity(args.ns)
    alamda = -1.0
    orig_cxx_values = args.a.copy()
    for i in range(args.iterations):
        rus.mrqmin(args.order,problem_size,tabs,irk,dimensions,args.rho,args.shape,args.freqmin,y,num_freqs,args.a,ia,args.ns,covar,alpha,args.hextype,alamda)
        print('iter #{}'.format(i))
        for k in sorted(args.a.keys()):
            print('{0:.6f}'.format(100 * args.a[k])) # print estimated cxx values

    print('\nThis calculation can be executed again with the following command:')
    print('python {} --order {} --shape {} --ns {} --hextype {} --d1 {} --d2 {} --d3 {} --rho {} --freqmin {} --freqmax {} --iterations {} {}'.format(sys.argv[0], args.order, args.shape, args.ns, args.hextype, args.d1, args.d2, args.d3, args.rho, args.freqmin, args.freqmax, args.iterations, ' '.join(('--' + k + ' ' + str(v*100)) for k,v in orig_cxx_values.iteritems())))

def read_input(infile):

    if infile == 'sample/default_frequencies':
        print('!! no frequency file specified - sample file will be used !!')

    # get measured frequencies from file
    try:
        file_handle = open(infile, "rU")
    except IOError:
        print('Could not open frequency file.')
        sys.exit(-1)
    else:
        number_of_freqs = int(file_handle.readline())
        print('nfreq={}'.format(number_of_freqs))
        freq_list   = []
        for i in range(number_of_freqs):
            line = file_handle.readline()
            data = line.split(None, 1)
            if len(data) != 2:
                print('Could not parse some lines in the frequency file.')
                file_handle.close()
                sys.exit(-1)
            freq_list.append((float(data[0]), float(data[1])))
            print(' freq={0:.6f}'.format(float(data[0])))
        file_handle.close()

    if len(freq_list) != number_of_freqs:
        print('Unexpected number of frequencies.')
        print('Expected {} but read {}'.format(number_of_freqs, len(freq_list)))
        sys.exit(-1)

    return freq_list



def print_args(args, ndata):

    print('d={}'.format(args.order))
    print('shape={}'.format(args.shape))
    print('ns={}'.format(args.ns))
    print('hextype={}'.format(args.hextype))
    print('d1={0:.6f}'.format(args.d1))
    print('d2={0:.6f}'.format(args.d2))
    print('d3={0:.6f}'.format(args.d3))
    print('rho={0:.6f}'.format(args.rho))
    print('freqmin={0:.6f}'.format(args.freqmin))
    print('freqmax={0:.6f}'.format(args.freqmax))

    for key, value in args.a.iteritems():
        print('{0:.6f}'.format(value))
        args.a[key] = value / 100

    print('ndata={}'.format(str(ndata)))
