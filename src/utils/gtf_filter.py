import os
import sys

from Bio import SeqIO
from ..pipeline_config import PipelineStructure


class GTFFilter(PipelineStructure):
    def __init__(self, args, gtf, fasta):
        super().__init__(args=args)
        self.gtf = gtf
        self.fasta = fasta

        self.entries = []
        self.tmpDir = f'{self.outdir}/tmp_gtf'

    def get_entries(self):
        records = SeqIO.parse(self.fasta, 'fasta')
        for record in records:
            self.entries.append(str(record.description))

    def filter_gtf(self, output):
        new_lines = []
        with open(self.gtf, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if any(entry in line for entry in self.entries):
                    new_lines.append(line)
        with open(output, 'w') as outfile:
            outfile.writelines(new_lines)

    def filter_gtf_tmp(self):
        self.check_dirs([self.tmpDir])
        for i, entry in enumerate(self.entries):
            cmd = f'grep "{entry}" {self.gtf} > {self.tmpDir}/entry_{i}_tmp.txt'
            os.system(cmd)
        print("merging tmp files")
        cmd_cat = f'cat {self.tmpDir}/* > {self.rescoredMicroproteinsGTF}'
        os.system(cmd_cat)
        cmd_rm = f'rm {self.tmpDir}/*'
        os.system(cmd_rm)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("Filters a GTF based on entries present in a fasta file. \n"
              "usage: .py <gtf> <fasta> <outfile>")

    else:
        data = GTFFilter(gtf=sys.argv[1], fasta=sys.argv[2])
        data.get_entries()
        data.filter_gtf(output=sys.argv[3])


