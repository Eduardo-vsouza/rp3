import os
import sys


def check_hits(utps, blast, out):
    entries = []
    with open(utps, 'r') as handler:
        lines = handler.readlines()
        for line in lines:
            entries.append(line.rstrip())

    blast_utps = []
    with open(blast, 'r') as blast_handler:
        lines = blast_handler.readlines()
        for line in lines:
            if any(entry in line for entry in entries):
                blast_utps.append(line)
    with open(out, 'w') as outfile:
        outfile.writelines(blast_utps)


if __name__ == '__main__':
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("usage: .py <utps> <blast> <out>")
    else:
        check_hits(utps=sys.argv[1], blast=sys.argv[2], out=sys.argv[3])