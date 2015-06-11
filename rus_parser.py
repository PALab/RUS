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
    c_help  = 'The order of cxxs depends on the symmetry type. Isotropic: c11, c44. Cubic: c11, c12, c44. Hexagonal VTI: c33, c23, c12, c44, c66. Hexagonal HTI: c11, c33, c12, c44, c66. Orthorhombic: c11, c22, c33, c23, c13, c12, c44, c55, c66. (Default: [110.0, 26.0])'

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

    parser.add_argument(
        '-c', '--cxxs',
        type = float,
        nargs = '+',
        default = [110.0, 26.0],
        help = c_help)

    args = parser.parse_args(params[1:])

    if args.ns == 2 and len(args.cxxs) != 2:
        parser.error('Please provide cxx values for c11 and c44')
    if args.ns == 3 and len(args.cxxs) != 3:
        parser.error('Please provide cxx values for c11, c12, and c44')
    if args.ns == 5 and len(args.cxxs) != 5 and args.hextype == 1:
        parser.error('Please provide cxx values for c33, c23, c12, c44, and c66')
    if args.ns == 5 and len(args.cxxs) != 5 and args.hextype == 2:
        parser.error('Please provide cxx values for c11, c33, c12, c44, and c66')
    if args.ns == 9 and len(args.cxxs) != 9:
        parser.error('Please provide cxx values for c11, c22, c33, c23, c13, c12, c44, c55, and c66')

    return args
