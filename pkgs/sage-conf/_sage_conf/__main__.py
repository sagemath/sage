# Entry point 'sage-config'.  It does not depend on any packages.

def _main():
    from argparse import ArgumentParser
    from sys import exit, stdout

    import sage_conf

    parser = ArgumentParser(prog='sage-config')
    try:
        from sage_conf import VERSION
    except ImportError:
        VERSION = 'none'

    parser.add_argument('--version', help="show version", action="version",
                       version='%(prog)s ' + VERSION)
    parser.add_argument("VARIABLE", nargs='?', help="output the value of VARIABLE")
    args = parser.parse_args()
    if args.VARIABLE:
        stdout.write('{}\n'.format(getattr(sage_conf, args.VARIABLE)))
    else:
        for k in dir(sage_conf):
            if not k.startswith('_'):
                stdout.write('{}={}\n'.format(k, getattr(sage_conf, k)))
