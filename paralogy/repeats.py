import os
import sys

import pandas as pd

from ..pipeline_config import PipelineStructure


class RepeatsDatabase(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.check_dirs(self.repeatsDBDir)
        self.check_dirs([self.repeatsProteogenomicsDBs])

        self.repeatsGTF = f'{self.repeatsDBDir}/repeats.gtf'

    def convert_repeats_to_gtf(self):
        df = pd.read_csv(self.args.repeatsFile, sep='\t')
        dict_df = df.to_dict(orient="list")
        gtf = []
        for i, chrom in enumerate(dict_df["query_seq"]):
            start = dict_df["query_start"][i]
            end = dict_df["query_end"][i]
            strand = dict_df["strand"][i]
            name = f'{dict_df["matching_repeat"][i]}_{dict_df["repeat_class_family"][i]}_{dict_df["id"][i]}'
            attrs = f'gene_id "{name}"; transcript_id "{name}";'
            line = f'{chrom}\tRp3_repeatMasker\ttranscript\t{start}\t{end}\t.\t{strand.replace("C", "-")}\t.\t{attrs}\n'
            gtf.append(line)
        with open(self.repeatsGTF, 'w') as handler:
            handler.writelines(gtf)

    # def