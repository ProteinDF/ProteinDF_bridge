import sys
import argparse
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge

def main():
    # parse args
    parser = argparse.ArgumentParser(description='crystallize molecules')
    parser.add_argument('INPUT_BRD_PATH',
                        nargs=1,
                        help='input bridge file path')
    parser.add_argument('OUTPUT_BRD_PATH',
                        nargs=1,
                        help='output bridge file path')

    parser.add_argument("-x",
                        nargs=1,
                        help="number of molecules for x")
    parser.add_argument("-y",
                        nargs=1,
                        help="number of molecules for y")
    parser.add_argument("-z",
                        nargs=1,
                        help="number of molecules for z")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    args = parser.parse_args()

    print(args)


if __name__ == '__main__':
    main()
