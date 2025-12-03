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
        if self.args.diamondBlastDBFasta:
            blast.parse(evalue=0.001, score=50, pc_id=0, qcov=0, format='diamond', return_id=True, conservation=True)
            blast.add_sequences_from_db(diamonblastdb_fasta=self.args.diamondBlastDBFasta, outdir=self.phyloDir)
        else:
            blast.parse(evalue=0.001, score=50, pc_id=0, qcov=0, conservation=True, format=self.args.blastType)

        # blast.parse(evalue=0.001, score=50, pc_id=0, qcov=0, conservation=True, format=self.args.blastType)

        blast.save_conserved(output=self.homologsPerSpecies)
        blast.create_spreadsheet(output=f'{self.phyloDir}/smorfs_entries_per_species.xls')
        # blast.add_sequences_from_db(output=f'{self.phyloDir}/smorfs_entries_per_species_with_seqs.xls')
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

    def generate_data_frame(self, external_df=None):
        if self.args.externalDFtoMerge is not None and external_df is None:
            # Load external DF once; keep a copy for the DB-seq merge
            external_df = pd.read_csv(self.args.externalDFtoMerge, sep='\t', low_memory=False)
            ext_df_for_db = external_df.copy()
        else:
            ext_df_for_db = external_df.copy() if external_df is not None else None

        df = pd.read_csv(f'{self.phyloDir}/smorfs_entries_per_species.xls', sep='\t')
        # Create a pivot table with smorf as rows and species as columns
        df_pivot = pd.crosstab(df['smorf'], df['species'])
        # Convert counts to boolean values (True when a species is present, False otherwise)
        df_boolean = df_pivot.applymap(lambda x: True if x > 0 else False)
        sp = self.args.spOrigin.replace("_", " ").capitalize()
        if sp in df_boolean.columns:
            df_boolean.drop(columns=[sp], inplace=True)
        df_boolean.to_csv(f'{self.phyloDir}/smorfs_entries_per_species_boolean.csv', sep='\t')
        
        # --- Second dataframe with smorf → comma-separated species ---
        df_summary = (
            df.groupby('smorf')['species']
              .apply(lambda x: ', '.join(sorted(set(x))))
              .reset_index()
        )
        df_summary.columns = ['smorf', 'conserved_species']
        
        # --- If external_df is provided, merge on provided id col (species summary) ---
        if external_df is not None:
            id_col = getattr(self.args, 'externalDFidCol', 'ID')
            if id_col in external_df.columns:
                # normalize keys to strings
                external_df[id_col] = (
                    external_df[id_col].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
                )
                df_summary['smorf'] = df_summary['smorf'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)

                merged_species = external_df.merge(df_summary, left_on=id_col, right_on='smorf', how='left')
                merged_species.to_csv(
                    f'{self.phyloDir}/merged_with_species_summary.csv',
                    sep='\t',
                    index=False
                )
                print(f"--Merged conservation summary with external DataFrame → merged_with_species_summary.csv")
            else:
                print("Warning: external_df does not have the provided id column; skipping species-summary merge.")
        
        # Always write the standalone summary
        df_summary.to_csv(
            f'{self.phyloDir}/smorfs_species_summary.csv',
            sep='\t',
            index=False
        )

        # --- NEW: Merge the DB sequences table (smorf, hit_id, sequence) with the external DF by id col ---
        seq_path = f'{self.phyloDir}/smorfs_conserved_sequences_from_db.tsv'
        if ext_df_for_db is not None and os.path.exists(seq_path):
            id_col = getattr(self.args, 'externalDFidCol', 'ID')

            seq_df = pd.read_csv(seq_path, sep='\t', low_memory=False)
            # Ensure expected columns exist; if not, try to coerce first 3
            expected = {'smorf', 'hit_id', 'sequence'}
            if not expected.issubset(set(seq_df.columns)):
                if seq_df.shape[1] >= 3:
                    seq_df = seq_df.iloc[:, :3]
                    seq_df.columns = ['smorf', 'hit_id', 'sequence']
                else:
                    print(f"Warning: {seq_path} has unexpected columns; skipping DB-seq merge.")
                    return

            # Normalize types
            for c in ('smorf', 'hit_id', 'sequence'):
                seq_df[c] = seq_df[c].astype(str).str.strip()

            # Aggregate to one row per smorf to avoid row explosion
            agg_df = (
                seq_df.groupby('smorf')
                      .agg(
                          db_hit_ids=('hit_id', lambda s: '; '.join(sorted(set(s.dropna())))),
                          db_sequences=('sequence', lambda s: '; '.join(sorted(set(s.dropna())))),
                      )
                      .reset_index()
            )

            # Normalize keys for merge
            agg_df['smorf'] = agg_df['smorf'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
            if id_col in ext_df_for_db.columns:
                ext_df_for_db[id_col] = (
                    ext_df_for_db[id_col].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
                )
                merged_db = ext_df_for_db.merge(agg_df, left_on=id_col, right_on='smorf', how='left')
                # Drop redundant smorf if different key name
                if id_col != 'smorf' and 'smorf' in merged_db.columns:
                    merged_db.drop(columns=['smorf'], inplace=True)
                merged_db.to_csv(f'{self.phyloDir}/merged_with_db_sequences.tsv', sep='\t', index=False)
                print(f"--Merged external DataFrame with DB sequences → merged_with_db_sequences.tsv")
            else:
                print(f"Warning: external_df does not have id column '{id_col}'; skipping DB-seq merge.")
        elif ext_df_for_db is not None:
            print(f"Warning: sequences table not found at {seq_path}; skipping DB-seq merge.")

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
