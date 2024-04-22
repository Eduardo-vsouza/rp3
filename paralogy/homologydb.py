import os
import sys

import pandas as pd
from Bio import SeqIO

from ..pipeline_config import PipelineStructure
from ..utils import BlastParser


class HomologyDB(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        # directories
        self.blastDBDir = f'{self.homologyDBDir}/blastdb'
        self.proteogenomicsDatabases = f'{self.homologyDBDir}/proteogenomics_databases'
        self.check_dirs([self.homologyDBDir, self.blastDBDir, self.proteogenomicsDatabases])

        # files
        self.genomeBlastDB = f'{self.blastDBDir}/genome_db'
        self.smorfsBlastedToGenome = f'{self.homologyDBDir}/smorfs_blasted_to_genome.xml'
        self.homologyDF = f'{self.homologyDBDir}/smorfs_genome_homologs.txt'

    def create_blast_db(self):
        print(f"--Generating blastdb for the genome.")
        if not os.path.exists(self.genomeBlastDB):
            cmd = f'makeblastdb -in {self.args.genome} -dbtype nucl -out {self.genomeBlastDB}'
            os.system(cmd)
        else:
            print(f"Found genome blastdb. Skipping blastdb generation...")

    def blast(self):
        print(f"--Identifying highly homologous sequences for 3-frame translated smORFs.")
        assemblies = os.listdir(self.translationDir)
        for ass in assemblies:
            files = os.listdir(f'{self.translationDir}/{ass}')
            for file in files:
                if file.endswith("split_nuc"):
                    if not os.path.exists(self.smorfsBlastedToGenome):
                        smorfs = f'{self.translationDir}/{ass}/{file}'
                        print(f"---Blasting predicted smORFs against the genome.")
                        cmd = (f'blastn -query {smorfs} -db {self.genomeBlastDB} -evalue 0.001 -outfmt 5 -out '
                               f'{self.smorfsBlastedToGenome} -task blastn-short -soft_masking false -dust no '
                               f'-num_threads {self.args.threads}')
                        os.system(cmd)
                    else:
                        print(f"---Found blast XML for {file}. Skipping step...")

    def identify_homologs(self):
        print(f"--Parsing blast results")
        if not os.path.exists(self.homologyDF):
            blast = BlastParser(xml=self.smorfsBlastedToGenome)
            blast.parse(pc_id=0, qcov=0, paralogs=True)
            blast.save_paralogs(outfile=self.homologyDF)
        else:
            print(f"---Found homology data frame. Skipping step...")

    def filter_highly_homologous(self):
        print(f"--Selecting smORFs with high homology to the rest of the genome.")
        df = pd.read_csv(self.homologyDF, sep='\t')
        df = df[df["number_of_homologs"] >= self.args.minHomologs]
        smorfs = df["smorfs"].tolist()
        assemblies = os.listdir(f'{self.translationDir}')
        for ass in assemblies:
            files = os.listdir(f'{self.translationDir}/{ass}')
            for file in files:
                if file.endswith(".pep"):
                    out_ass_dir = f'{self.proteogenomicsDatabases}/{ass}'
                    self.check_dirs([out_ass_dir])
                    out_db = f'{self.proteogenomicsDatabases}/{ass}/{file[:-4]}_highHomologyDB.pep'
                    if not os.path.exists(out_db):
                        out_fasta = []
                        records = SeqIO.parse(f'{self.translationDir}/{ass}/{file}', 'fasta')
                        for record in records:
                            entry = str(record.description)
                            if entry in smorfs:
                                out_fasta.append(f'>{entry}\n{str(record.seq)}\n')
                        with open(out_db, 'w') as out:
                            out.writelines(out_fasta)
                    else:
                        print(f"---Found homology DB for {file}. Skipping step...")


