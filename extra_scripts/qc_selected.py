import os
import sys

import pandas as pd


class Redo:
    def __init__(self, target_files, fastq_folder, outdir):
        self.targetFiles = pd.read_csv(target_files, sep='\t')
        self.targetFileNames = self.__get_target_files()

        self.fastqFolder = fastq_folder
        self.outdir = outdir
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    def __get_target_files(self):
        return self.targetFiles["file"].tolist()

    def run_fastqc(self):
        files = os.listdir(self.fastqFolder)
        input_files = ''
        i = 0
        for file in files:
            if any(f in file for f in self.targetFileNames):
                i += 1
                input_files += f' {self.fastqFolder}/{file}'
        cmd = f'/home/microway/miniconda3/bin/fastqc --threads {i} -o {self.outdir}{input_files}'
        print(cmd)
        os.system(cmd)


    def move_files_to_realign(self):
        files = os.listdir(self.fastqFolder)
        for file in files:
            if any(f in file for f in self.targetFileNames):
                cmd = f'mv {self.fastqFolder}/{file} {self.outdir}/{file}'
                os.system(cmd)


if __name__ == '__main__':
    # data = Redo(target_files='brendan_files_to_realign',
    #             fastq_folder='/ceph/pbla/brendan/human_brain_riboseq/fastq',
    #             outdir='/ceph/pbla/brendan/human_brain_riboseq/fastq/files_to_redo')
    # data.run_fastqc()

    ### move files to re-align
    data = Redo(target_files='brendan_files_to_realign',
                fastq_folder='/ceph/pbla/brendan/human_brain_riboseq/fastq',
                outdir='/ceph/pbla/brendan/human_brain_riboseq/fastq/25-01-28_files_to_realign')
    data.move_files_to_realign()
