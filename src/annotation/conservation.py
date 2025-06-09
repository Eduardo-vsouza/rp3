import collections
import os
import sys

from Bio import SeqIO
import pandas as pd

from ..pipeline_config import PipelineStructure
from ..utils import BlastParser


class Conservation(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)


        self.check_dirs([self.phyloDir, self.blastDir])

        self.blastDatabase = self.args.blast_db

        self.blastedMicroproteins = f'{self.blastDir}/microproteins_blasted_to_db.xml'
        self.homologsPerSpecies = f'{self.phyloDir}/homologs_per_species.xls'
        self.nonRedundantFasta = f'{self.mergedResults}/microproteins_blast_filt_non_redundant.fasta'

    def generate_non_redundant_fasta(self):
        print("--Generating non-redundant fasta\n")
        fasta = []
        checker = []
        file = self.microproteinsBlast
        if self.args.rescored:
            file = self.rescoredMicroproteinsFasta
        if self.args.externalFasta is not None:
            file = self.args.externalFasta
            self.nonRedundantFasta = f'{self.outdir}/external_fasta_non_redundant.fasta'
        records = SeqIO.parse(f'{file}', 'fasta')
        for record in records:
            seq = str(record.seq)
            if '>' in seq:
                continue
                # print(str(record.description), seq)
            if seq not in checker:
                checker.append(seq)
                fasta.append(f'>{str(record.description)}\n{seq}\n')
        with open(self.nonRedundantFasta, 'w') as out:
            out.writelines(fasta)

    def blast_microproteins(self):
        print("--BLASTing microproteins to conservation database\n")
        if self.args.blastType == 'tblastn':
            cmd = f'{self.toolPaths["tblastn"]} -outfmt 5 -out {self.blastedMicroproteins} -evalue 0.001 ' \
                f' -num_threads {self.args.threads} -query {self.nonRedundantFasta} ' \
                f'-db {self.blastDatabase}'
        elif self.args.blastType == 'diamond':
            print(f"--Running diamondBlast")
            cmd = f'{self.toolPaths["diamond"]} blastp -q {self.nonRedundantFasta} -d {self.args.diamondDB} ' \
                f'--outfmt 5 --out {self.blastedMicroproteins} --threads {self.args.threads} --ultra-sensitive'
        print(f"command: {cmd}")
        os.system(cmd)

    def parse_blast_results(self):
        print("--Parsing blast results\n")
        blast = BlastParser(xml=self.blastedMicroproteins)
        blast.parse(evalue=0.001, score=50, pc_id=0, qcov=0, conservation=True, format=self.args.blastType)
        blast.save_conserved(output=self.homologsPerSpecies)
        blast.create_spreadsheet(output=f'{self.phyloDir}/smorfs_entries_per_species.xls')
        print("Finished parsing\n")

    def create_evolview_input(self):
        print(f"--Creating evolview input\n")
        evol_lines = ['!groups	group1\n',
                      '!colors	#81DBFD\n',
                      '!PlotWidth	500\n',
                      '!grid\n',
                      "!axis\n",
                      '!itemHeightPCT	80\n',
                      '!showdataValue	show=1,fontsize=10,fontitalic=0,textalign=end,fontcolor=black\n']
        with open(self.homologsPerSpecies, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if 'Species' not in line:
                    evol_lines.append(line)
        with open(f'{self.phyloDir}/homologs_per_species_evolview.txt', 'w') as outfile:
            outfile.writelines(evol_lines)

    def classify_conservation_by_mapping_groups(self):
        groups_df = pd.read_csv(f'{self.microproteinMappingGroupsForPlotsUnion}', sep='\t')
        classification = {}
        smorfs, groups = groups_df["smorf"].tolist(), groups_df["group"].tolist()
        mus_no_cov = []
        mus_rp3 = []

        for smorf, group in zip(smorfs, groups):
            if group == "No coverage":
                classification[smorf] = group
                mus_no_cov.append(smorf)
            else:
                mus_rp3.append(smorf)
                classification[smorf] = "RP3"
        homologs = {}

        homologs_df = pd.read_csv(f'{self.phyloDir}/smorfs_entries_per_species.xls', sep='\t')
        species, smorfs = homologs_df["species"].tolist(), homologs_df["smorf"].tolist()
        for sp, smorf in zip(species, smorfs):
            if sp not in homologs:
                homologs[sp] = {'No coverage': [], "RP3": []}
            homologs[sp][classification[smorf]].append(smorf)
        print(homologs)
        evol_no_cov = ['!groups	No coverage\n',
                       '!colors	#81DBFD\n',
                       '!PlotWidth	500\n',
                       '!grid\n',
                       "!axis\n",
                       '!itemHeightPCT	80\n',
                       '!showdataValue	show=1,fontsize=10,fontitalic=0,textalign=end,fontcolor=black\n']
        evol_rp3 = ['!groups	RP3\n',
                    '!colors	#81DBFD\n',
                    '!PlotWidth	500\n',
                    '!grid\n',
                    "!axis\n",
                    '!itemHeightPCT	80\n',
                    '!showdataValue	show=1,fontsize=10,fontitalic=0,textalign=end,fontcolor=black\n']
        for homolog in homologs:
            for group in homologs[homolog]:
                if homolog == 'Homo sapiens':
                    if group == 'No coverage':
                        evol_no_cov.append(f'{homolog}\t{len(mus_no_cov)}\n')
                    else:
                        evol_rp3.append(f'{homolog}\t{len(mus_rp3)}\n')
                else:
                    print(homolog, group)
                    print(len(homologs[homolog][group]))
                    if group == 'No coverage':
                        evol_no_cov.append(f'{homolog}\t{len(homologs[homolog][group])}\n')
                    else:
                        evol_rp3.append(f'{homolog}\t{len(homologs[homolog][group])}\n')
        with open(f'{self.phyloDir}/evolview_no_cov.txt', 'w') as out_no_cov, open(f'{self.phyloDir}/evolview_rp3.txt', 'w') as out_rp3:
            out_no_cov.writelines(evol_no_cov)
            out_rp3.writelines(evol_rp3)
