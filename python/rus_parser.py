import argparse

def inverse_parser(params):
    # Parser help text
    d_help  = 'order of polynomials used to estimate the eigenvectors (Default: 8)'
    s_help  = '0=sphere, 1=cylinder, 2=parallelepiped (Default: 1)'
    n_help  = 'number of cijs (Default: 2)'
    t_help  = 'hextype - 1=VTI, 2=HTI. Type of hexagonal symetry (Only matters for ns=5) (Default: 1)'
    x_help  = 'dimension 1 in cm (diameter for cyl. or sphere) (Default: 4.420)'
    y_help  = 'dimension 2 in cm (diameter for cyl. or sphere) (Default: 4.420)'
    z_help  = 'dimension 3 in cm (height for cyl. diameter for sphere) (Default: 6.414)'
    r_help  = 'density in grams/cm^3 (Default: 2.713)'
    l_help  = 'lower frequency bound for inversion in MHz (set >1 KHz lower than your lowest measured value) (Default: 0.020)'
    u_help  = 'upper frequency bound for inversion in MHz (set >5 or 10KHz higher than your highest value used for THIS particular fit as defined by Line 1 of freq_data) (Default: 0.110)'

    parser = argparse.ArgumentParser(description='Inverse Algorithm')

    parser.add_argument(
        '-d', '--order',
        type = int,
        default = 8,
        help = d_help)

    parser.add_argument(
        '-s', '--shape',
        type = int,
        default = 1,
        choices = [0,1,2],
        help = s_help)

    parser.add_argument(
        '-n', '--ns',
        type = int,
        default = 2,
        choices = [2,3,5,9],
        help = n_help)

    parser.add_argument(
        '-t', '--hextype',
        type = int,
        default = 1,
        choices = [1,2],
        help = t_help)

    parser.add_argument(
        '-x', '--d1',
        type = float,
        default = 4.420,
        help = x_help)

    parser.add_argument(
        '-y', '--d2',
        type = float,
        default = 4.420,
        help = y_help)

    parser.add_argument(
        '-z', '--d3',
        type = float,
        default = 6.414,
        help = z_help)

    parser.add_argument(
        '-r', '--rho',
        type = float,
        default = 2.713,
        help = r_help)

    parser.add_argument(
        '-l', '--freqmin',
        type = float,
        default = 0.020,
        help = l_help)

    parser.add_argument(
        '-u', '--freqmax',
        type = float,
        default = 0.110,
        help = u_help)

    parser.add_argument('--c11', type = float)
    parser.add_argument('--c12', type = float)
    parser.add_argument('--c13', type = float)
    parser.add_argument('--c22', type = float)
    parser.add_argument('--c23', type = float)
    parser.add_argument('--c33', type = float)
    parser.add_argument('--c44', type = float)
    parser.add_argument('--c55', type = float)
    parser.add_argument('--c66', type = float)

    args = parser.parse_args(params[1:])

    args.a = dict()
    if args.c11: args.a['c11'] = args.c11
    if args.c12: args.a['c12'] = args.c12
    if args.c13: args.a['c13'] = args.c13
    if args.c22: args.a['c22'] = args.c22
    if args.c23: args.a['c23'] = args.c23
    if args.c33: args.a['c33'] = args.c33
    if args.c44: args.a['c44'] = args.c44
    if args.c55: args.a['c55'] = args.c55
    if args.c66: args.a['c66'] = args.c66

    # set default a
    if len(args.a) == 0:
        args.a['c11'] = 110.0
        args.a['c44'] = 26.0
    
    return args
