import os
import sys

import pandas as pd
from Bio import SeqIO


class RP3Intersection:
    def __init__(self, outfile, outdirs):
        self.outfile = outfile
        self.outdirs = outdirs

        self.pathoMicroproteins = {}
        self.hostMicroproteins = {}

        self.hostSharedMPs = {}


    def gather_microproteins(self, host_pattern='chr', patho_pattern='PCMN'):
        for outdir in self.outdirs:
            if os.path.exists(f'{outdir}/duo_proteogenomics/microproteins_report.txt'):
                duodir = f'{outdir}/duo_proteogenomics'
                host_records = SeqIO.parse(f'{duodir}/{host_pattern}_microproteins.fasta',
                                           'fasta')
                patho_records = SeqIO.parse(f'{duodir}/{patho_pattern}_microproteins.fasta', 'fasta')
                self.__iterate_fasta(records=host_records, outdir=outdir, dictio=self.hostMicroproteins)
                self.__iterate_fasta(records=patho_records, outdir=outdir, dictio=self.pathoMicroproteins)
            else:
                patho_records = SeqIO.parse(f'{outdir}/summarized_results/merged/microproteins_150.fasta',
                                            'fasta')
                self.__iterate_fasta(records=patho_records, outdir=outdir, dictio=self.pathoMicroproteins)

    def __iterate_fasta(self, records, outdir, dictio):
        for record in records:
            if str(record.seq) not in dictio:
                dictio[str(record.seq)] = []
            dictio[str(record.seq)].append(outdir)

    def generate_data_frames(self, outfile_pattern):

        def feed_table(dictio):
            data = {'microproteins': []}
            for outdir in self.outdirs:
                data[outdir] = []

            for mp in dictio:
                data['microproteins'].append(mp)
                for outdir in self.outdirs:
                    if outdir in dictio[mp]:
                        data[outdir].append('TRUE')
                    else:
                        data[outdir].append('FALSE')
            return data

        host_data = feed_table(self.hostMicroproteins)
        patho_data = feed_table(self.pathoMicroproteins)
        host_df = pd.DataFrame(data=host_data)
        patho_df = pd.DataFrame(data=patho_data)
        host_df.to_csv(f'{outfile_pattern}_host_microproteins.xls', sep='\t', index=False)
        patho_df.to_csv(f'{outfile_pattern}_viral_microproteins.xls', sep='\t', index=False)


if __name__ == '__main__':
    data = RP3Intersection(outfile=sys.argv[1], outdirs=sys.argv[2:])
    data.gather_microproteins()
    data.generate_data_frames(outfile_pattern=sys.argv[1])


