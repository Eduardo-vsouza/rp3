import os
import sys

import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure


class MHCDetector(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.mhcDir = f'{self.outdir}/mhcFlurry'
        self.check_dirs([self.mhcDir])
        self.mhcMicroproteins = f'{self.mhcDir}/mhc_microproteins.fasta'

    def run_mhc_flurry(self):
        if os.path.exists(self.rescoredMicroproteinsFasta):
            fasta = self.rescoredMicroproteinsFasta
        else:
            fasta = self.uniqueMicroproteins
        alleles = ('HLA-A*02:01,HLA-A*03:01,HLA-B*57:01,HLA-B*45:01,HLA-C*02:02,HLA-C*07:02'
                   ' HLA-A*01:01,HLA-A*02:06,HLA-B*44:02,HLA-B*07:02,HLA-C*01:02,HLA-C*03:01')
        # 7-12 is the length in MSBooster paper for HLA peptidomics
        cmd = (f'mhcflurry-predict-scan {fasta} --alleles {alleles} --peptide-lengths 7-12 --out'
               f' {self.mhcDir}/predictions.txt')
        os.system(cmd)

    def filter_results(self):
        df = pd.read_csv(f'{self.mhcDir}/predictions.txt', sep=',')
        df = df[df["affinity"] <= float(self.args.affinity)]
        df = df[df["affinity_percentile"] <= float(self.args.affinityPercentile)]
        smorfs = df["sequence_name"].tolist()

        fasta = self.select_fasta()
        records = SeqIO.parse(fasta, 'fasta')
        out_fasta = []
        seqs = []
        for record in records:
            entry = str(record.description)
            if entry in smorfs:
                seqs.append(str(record.seq))
                out_fasta.append(f'>{entry}\n{str(record.seq)}\n')
        with open(self.mhcMicroproteins, 'w') as outfile:
            outfile.writelines(out_fasta)
        print(f"A total of {len(set(seqs))} containing immunogenic peptides with affinity <= {self.args.affinity} "
              f"and an affinity percentile <= {self.args.affinityPercentile} were identified.")
        self.__save_report()

        if self.args.filterPipeResults:
            cmd = f'mv {fasta} {fasta}_unfiltered_MHC.fasta'
            os.system(cmd)
            cmd = f'cp {self.mhcMicroproteins} {fasta}'
            os.system(cmd)

    def __save_report(self):
        lines = []
        records = SeqIO.parse(self.mhcMicroproteins, 'fasta')
        seqs = []
        for record in records:
            if str(record.seq) not in seqs:
                seqs.append(str(record.seq))
        unfilt_seqs = []
        records = SeqIO.parse(self.select_fasta(), 'fasta')
        for record in records:
            seq = str(record.seq)
            if seq not in unfilt_seqs:
                unfilt_seqs.append(seq)
        lines = [f'Number of identified non-redundant proteins before MHC filtering: {len(set(unfilt_seqs))}\n',
                 f'Number of identified non-redundant proteins after MHC filtering: {len(set(seqs))}\n']
        with open(f'{self.mhcDir}/analysis_report.txt', 'w') as out:
            out.writelines(lines)