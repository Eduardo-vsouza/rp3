import os
import sys
import inspect

import pandas as pd
from Bio import SeqIO

from .decoy import Decoy
from .paralogy import HomologyDB
from .pipeline_config import PipelineStructure


class Database(PipelineStructure):
    def __init__(self, translation_folder, reference_proteome, outdir, args, external_database=None):
        super().__init__(args=args)
        self.print_row(word="Database")
        self.translationFolder = translation_folder
        self.proteome = reference_proteome
        self.externalDatabase = external_database
        self.targetDatabases = {}
        self.args = args
        self.outdir = outdir
        self.databaseDir = f'{outdir}/databases'
        self.__check_dir()

        # params
        self.mode = 'database'
        self.params = []

        # uniprot annotation levels
        if self.args.uniprotAnnotation is not None:
            self.annotations = self.__get_annotation_levels()
        else:
            self.annotations = {}
        
        self.__save_ref_proteome()

    def __save_ref_proteome(self):
        self.print_state(message='Saving reference proteome to the database directory', level=1, color='yellow')
        ref_proteome = []
        records = SeqIO.parse(self.args.proteome, 'fasta')
        for record in records:
            seq = str(record.seq)
            entry = str(record.description).replace(" ", "_").replace(",", "_")
            ref_proteome.append(f'>{entry}\n{seq}\n')
        with open(self.refProteome, 'w') as out:
            out.writelines(ref_proteome)
        # print(f"--Reference proteome saved to {self.databaseDir}/proteome.fasta")
        self.print_state(message=f"Reference proteome saved to {self.refProteome}", level=2,
                         color='green')

        if self.args.cascade:
            # print(f"--Cascade mode detected. Generating decoy sequences for the reference proteome.")
            self.print_state(message='Cascade mode detected. Generating decoy sequences for the reference proteome.',
                             level=1, color='yellow')
            decoy = Decoy(db=self.refProteome)
            decoy.reverse_sequences().to_fasta(output=f'{self.refProteomeWithDecoy}',
                                                pattern='rev', merge=True)
            self.print_state(message=f"Reference proteome decoy sequences saved to {self.refProteomeWithDecoy} for using in Cascade mode",
                             color='green', level=2)

    def __check_dir(self):
        if not os.path.exists(self.databaseDir):
            os.mkdir(self.databaseDir)

    def prepare_external_database(self):

        if self.externalDatabase is not None:
            print(f"--Preparing external database.")
            if '/' in self.externalDatabase:
                file = self.externalDatabase.split("/")[-1]
            else:
                file = self.externalDatabase
            db_dir = f'{self.translationFolder}/{file.replace(".fasta", "")}'
            if not os.path.exists(db_dir):
                os.mkdir(db_dir)
            out_fasta = []
            records = SeqIO.parse(self.externalDatabase, 'fasta')
            for record in records:
                seq = str(record.seq)
                entry = str(record.description).replace(" ", "_").replace(",", "_")
                out_fasta.append(f'>{entry}\n{seq}\n')
            with open(f'{db_dir}/{file.replace(".fasta", "")}.pep', 'w') as out:
                out.writelines(out_fasta)
            # cmd = f'cp {self.externalDatabase} {db_dir}/{file.replace(".fasta", "")}.pep'
            # os.system(cmd)

    def unzip_assemblies(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        assemblies = os.listdir(self.translationFolder)
        for assembly in assemblies:
            cmd = f'gzip -d {self.translationFolder}/{assembly}/*.gz'
            self.params.append(cmd)
            os.system(cmd)

    def select_highly_homologous(self):
        if self.args.highHomologyDB:
            db = HomologyDB(args=self.args)
            db.create_blast_db()
            db.blast()
            db.identify_homologs()
            db.filter_highly_homologous()

    def append_reference(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        folder = self.translationFolder
        if self.args.highHomologyDB:
            folder = self.homologyDBDir
        added_seqs = []
        assemblies = os.listdir(folder)
        for assembly in assemblies:
            if assembly not in self.targetDatabases:
                self.targetDatabases[assembly] = []
            pep = f'{folder}/{assembly}/{assembly}.pep'
            predicted = SeqIO.parse(pep, 'fasta')
            if not self.args.cascade:
                reference = SeqIO.parse(self.proteome, 'fasta')
                for record in reference:
                    # check if including uncharacterized microproteins or not
                    # only MPs below the length threshold and annotation level will be considered
                    if self.args.uniprotAnnotation is not None:
                        seq = str(record.seq)
                        if seq in self.annotations:
                            if int(self.annotations[seq]) <= int(self.args.annotationLevel):
                                self.targetDatabases[assembly].append(
                                    f'>{str(record.description).replace(" ", "_").replace(",", "_")}_UNCH\n{str(record.seq)}\n')
                            else:
                                self.targetDatabases[assembly].append(
                                    f'>{str(record.description).replace(" ", "_").replace(",", "_")}_ANNO\n{str(record.seq)}\n')
                        else:
                            self.targetDatabases[assembly].append(
                                f'>{str(record.description).replace(" ", "_").replace(",", "_")}_ANNO\n{str(record.seq)}\n')
                            
                    else:
                        if str(record.seq) not in added_seqs:
                            added_seqs.append(str(record.seq))
                        self.targetDatabases[assembly].append(f'>{str(record.description).replace(" ", "_").replace(",", "_")}_ANNO\n{str(record.seq)}\n')

            for record in predicted:
                if str(record.seq) not in added_seqs:
                    self.targetDatabases[assembly].append(f'>{str(record.description)}\n{str(record.seq)}\n')


    def save_target_dbs(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        for assembly in self.targetDatabases:
            with open(f'{self.databaseDir}/{assembly}_target_database.fasta', 'w') as out:
                out.writelines(self.targetDatabases[assembly])

    def create_decoy_dbs(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        targets = os.listdir(self.databaseDir)
        for target in targets:
            if os.path.isdir(f'{self.databaseDir}/{target}'):
                continue
            if not target.endswith("_target_database.fasta"):
                continue
            decoy = Decoy(db=f'{self.databaseDir}/{target}')
            # if self.args.cat:
            decoy.reverse_sequences().to_fasta(output=f'{self.databaseDir}/{target}'.replace("_target_", "_target_decoy_"),
                                                pattern='rev', merge=True)
            # decoy.add_contaminants()
            # else:
            #     decoy.reverse_sequences().to_fasta(output=f'{self.databaseDir}/{target}'.replace("_target_", "_decoy_"),
                                                #    pattern='rev', merge=False)
        print(f"--Finished generating databases. You can safely ignore the numpy warning.")

    def __get_annotation_levels(self):
        annotations = {}
        df = pd.read_csv(self.args.uniprotAnnotation, sep='\t')
        seqs, annos = df["Sequence"].tolist(), df["Annotation"].tolist()
        for seq, anno in zip(seqs, annos):
            if len(seq) <= self.args.maxLength:
                if seq not in annotations:
                    annotations[seq] = int(anno)
                else:
                    if int(anno) > int(annotations[seq]):
                        annotations[seq] = int(anno)
        return annotations


    def split_databases(self):

        def split_db(db, outdir):
            out_fastas = {0: []}
            records = SeqIO.parse(db, 'fasta')
            i = 0
            split = 0
            for record in records:
                i += 1
                if i == self.args.splitDatabase + 1:
                    i = 0
                    split += 1
                    if split not in out_fastas:
                        out_fastas[split] = []
                out_fastas[split].append(f'>{str(record.description)}\n{str(record.seq)}\n')
            for key in out_fastas:
                with open(f'{outdir}/db_{key}.fasta', 'w') as out:
                    out.writelines(out_fastas[key])
                decoy = Decoy(db=f'{outdir}/db_{key}.fasta')
                decoy.reverse_sequences().to_fasta(output=f"{outdir}/db_{key}_decoy.fasta",
                                                   pattern='rev', merge=True, contaminants=False)
                cmd = f'rm {outdir}/db_{key}.fasta'
                os.system(cmd)

        if self.args.splitDatabase is not None:
            self.print_row(word="database splitting", character='-')
            self.check_dirs([self.splitDbDir, self.splitDbProteogenomicsDir, self.splitDbProteomeDir])
            # print(f"--Splitting databases into {self.args.splitDatabase} parts")
            print(f"--Splitting reference proteome")
            db = self.select_database(decoy=False, proteome=True)
            split_db(db=db, outdir=self.splitDbProteomeDir)

            db = self.select_database(decoy=False, proteome=False)
            print(f"--Splitting proteogenomics database")
            split_db(db=db, outdir=self.splitDbProteogenomicsDir)

        conta = Decoy()
        self.check_dirs([self.splitDbDir])
        conta.to_fasta(self.fullContaminantsDb)
        self.print_row(word="Finished")


        