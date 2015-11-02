#!/usr/bin/env python

import os
import sys
import argparse

def validate(args):
    """

    process & validates the output FASTQC output folders.
    prints the lines with the status

    :param args: args for the validator
    """

    print("running validator with args: [%s]" % str(args))

    root = args.fastqc_dir
    output_dirs = [os.path.join(root, d) for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]
    for dir in output_dirs:
        name = os.path.basename(dir)
        input_file = os.path.join(dir, args.summary_file)

        print("\nSUMMARY: Name: [%s] with status: [%s]" % (name, args.status))
        with open(input_file) as f:
            lines = f.readlines()

            for line in lines:
                line = line.strip()

                if args.status:
                    if args.status in line:
                        print line
                else:
                    print line

def main(argv):
    """
    parse input arguments & run the tool.

    :return:
    """

    parser = argparse.ArgumentParser(prog='fastqc_validator', description='fastqc validator')
    parser.add_argument('-d', '--fastqc-dir', required=True, help='FASTQC output directory')
    parser.add_argument('-f', '--summary-file', default='summary.txt', help='summary file name')
    parser.add_argument('-s', '--status', default=None, help='status of summary lines')

    try:
        # parse arguments
        args = parser.parse_args(argv[1:])

        # validate
        validate(args)
    except Exception as e:
        print("error encountered while running validator: [%s]" % e.message)
        return 1

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))