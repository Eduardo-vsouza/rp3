import os
import sys

import pandas
import pandas as pd
from Bio import SeqIO

from src.utils import group_folder_generator, check_multiple_dirs
from src.pipeline_config import PipelineStructure


class FormylSummarizer(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.modsDir = f'{self.outdir}/mods'
        check_multiple_dirs([self.modsDir])

        self.databases = {}
        self.__get_gene_names()

    def get_modified_proteins(self, mod='28.0104', qvalue=0.01):
        generator = group_folder_generator(self.postProcessDir)
        for content in generator:
            peptides = f'{content.dbDir}/peptides_fixed.txt'
            df = pd.read_csv(peptides, sep='\t')
            df = df[df["peptide"].str.contains(mod)]
            df = df[df["proteinIds"].str.contains("contaminant") == False]
            df = df[df["proteinIds"].str.startswith("rev_") == False]
            proteins = df["proteinIds"].tolist()
            genes = self.__get_protein_gene(db=content.db, proteins=proteins)
            df.insert(4, "gene", genes)
            group = content.group
            if not self.args.include_high_fdr:
                df = df[df["q-value"] <= qvalue]
            viral = df[df["proteinIds"].str.contains("_ANNO") == False]
            human = df[df["proteinIds"].str.contains("_ANNO") == True]
            db_dir = f'{self.modsDir}/{content.db}'
            if not os.path.exists(db_dir):
                os.mkdir(db_dir)
            viral.to_csv(f"{db_dir}/{group}_viral_formyl_proteins.xls", sep='\t', index=False)
            human.to_csv(f"{db_dir}/{group}_human_formyl_proteins.xls", sep='\t', index=False)

    def __get_protein_gene(self, db, proteins):
        genes = []
        for prot in proteins:
            if prot in self.databases[db]:
                gene = self.databases[db][prot]['gene']
                genes.append(gene)
            else:
                genes.append(prot)
        return genes

    def __get_gene_names(self):
        dbs = os.listdir(f'{self.databaseDir}')
        for db in dbs:
            if db.endswith(".fasta"):
                gene_names = {}
                # print(db)
                records = SeqIO.parse(f'{self.databaseDir}/{db}', 'fasta')
                for record in records:
                    entry = str(record.id)
                    full = str(record.description)
                    # print(full)
                    if 'GN=' in full:
                        if '_ANNO' in full:
                            try:
                                gene = full.split("_")[-4].split("=")[1]
                            except:
                                gene = ''
                        else:
                            if ' ' in full:
                                pattern = ' '
                            else:
                                pattern = '_'
                            splat = full.split(pattern)
                            if 'GN=' in splat:
                                gene = full.split(pattern)[-3].split("=")[1]
                            else:
                                gene = full.split(pattern)[-4].split("=")[1]

                            # gene = ''

                    else:
                        gene = 'not found'

                    gene_names[entry] = {'full': full, 'gene': gene}
                self.databases[db] = gene_names

    def summarize_data(self):
        dbs = os.listdir(self.modsDir)
        for db in dbs:
            db_dir = f'{self.modsDir}/{db}'
            if os.path.isdir(db_dir):
                files = os.listdir(db_dir)
                data = {'source_mass_spec': [], 'ids_human': [], 'ids_virus': [], 'ids_human_high_fdr': [],
                        'ids_virus_high_fdr': []}

                for file in files:
                    full = f'{self.modsDir}/{db}/{file}'
                    if not file.startswith("summarized") and file.endswith(".xls"):
                        print(file)
                        df = pd.read_csv(full, sep='\t')
                        source = file.split("_")[0]
                        organism = file.split("_")[1]
                        # print(df)
                        df = df.drop_duplicates(subset=["proteinIds"])
                        # print(df)
                        high_fdr_df = df[df["q-value"] > 0.01]
                        human_high, virus_high = self.__add_info(high_fdr_df)

                        df = df[df["q-value"] <= 0.01]
                        human, virus = self.__add_info(df)


                        data['source_mass_spec'].append(source)
                        data['ids_human'].append(human)
                        data['ids_virus'].append(virus)
                        data['ids_human_high_fdr'].append(human_high)
                        data['ids_virus_high_fdr'].append(virus_high)
                df = pd.DataFrame(data=data)
                df.to_csv(f'{db_dir}/summarized_data.xls', sep='\t', index=False)

    def __add_info(self, df):
        human = 0
        virus = 0
        proteins = df["proteinIds"].tolist()
        for prot in proteins:
            if '_ANNO' not in prot:
                virus += 1
            else:
                human += 1
        return human, virus






