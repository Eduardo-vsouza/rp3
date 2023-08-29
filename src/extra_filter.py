import os
import sys

from .utils import group_folder_generator, BlastParser
from .pipeline_config import PipelineStructure


class ExtraFilter(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.refSeq = self.args.refseq
        self.blastXML = f'{self.resultsDir}/unique_microproteins_blasted_to_refseq.xml'

    def blast_filter(self):
        file = f'{self.summarizedResultsDir}/merged/microproteins_150.fasta'
        cmd = f'{self.toolPaths["blastp"]} -query {file} -subject {self.refSeq} -outfmt 5 -out' \
              f' {self.blastXML} -num_threads {self.args.threads}'
        os.system(cmd)

    def filter_results(self):
        data = BlastParser(xml=self.blastXML)
        data.parse(evalue=1, score=1)
        folders = os.listdir(self.summarizedResultsDir)
        for folder in folders:
            if os.path.isdir(f'{self.summarizedResultsDir}/{folder}'):
                db_dir = f'{self.summarizedResultsDir}/{folder}'
                files = os.listdir(db_dir)
                for file in files:
                    if 'blast' not in file:
                        data.filter_peptide_fasta(fasta=f'{db_dir}/{file}', output=f'{db_dir}/{file}_blast_filt.fasta')
        # data.filter_peptide_fasta(fasta=f'{self.uniqueMicroproteins}', output=f'{self.summarizedResultsDir}/microproteins_150_blast_filt.fasta')
        # data.filter_peptide_fasta(fasta=f'{self.summarizedResultsDir}/microproteins_utps_150.fasta',
        #                           output=f'{self.summarizedResultsDir}/microproteins_utps_150_blast_filt.fasta')
