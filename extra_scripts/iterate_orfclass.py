import os
import sys


def iterate(results):
    for result in results:
        # cmd = f'python3 ~/PycharmProjects/rp3/rp3.py anno --outdir {result} --rescored --genomeAssembly hg38 --orfClass --threads 128'
        cmd = f'python3 ~/PycharmProjects/rp3/rp3.py homology --outdir {result} --threads 128 --genomeAssembly hg38 --tBlastN'
        os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: .py <results1> <results2> ... <results_n>")
    else:
        iterate(results=sys.argv[1:])