import os


def check_multiple_dirs(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)