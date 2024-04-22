import os
import sys


def move_ind(folder, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    subdirs = os.listdir(folder)
    for subdir in subdirs:
        os.system(f'mv {folder}/{subdir}/*ind {outdir}/.')


if __name__ == '__main__':
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("usage: .py ")
    else:
        move_ind(folder=sys.argv[1], outdir=sys.argv[2])