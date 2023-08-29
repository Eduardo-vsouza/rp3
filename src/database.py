import os
import sys
import inspect

from Bio import SeqIO

from .decoy import Decoy


class Database:
    def __init__(self, translation_folder, reference_proteome, outdir, args, external_database=None):
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

    def __check_dir(self):
        if not os.path.exists(self.databaseDir):
            os.mkdir(self.databaseDir)

    def prepare_external_database(self):
        if self.externalDatabase is not None:
            if '/' in self.externalDatabase:
                file = self.externalDatabase.split("/")[-1]
            else:
                file = self.externalDatabase
            db_dir = f'{self.translationFolder}/{file.replace(".fasta", "")}'
            if not os.path.exists(db_dir):
                os.mkdir(db_dir)
            cmd = f'cp {self.externalDatabase} {db_dir}/{file.replace(".fasta", "")}.pep'
            os.system(cmd)

    def unzip_assemblies(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        assemblies = os.listdir(self.translationFolder)
        for assembly in assemblies:
            cmd = f'gzip -d {self.translationFolder}/{assembly}/*.gz'
            self.params.append(cmd)
            os.system(cmd)

    def append_reference(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        assemblies = os.listdir(self.translationFolder)
        for assembly in assemblies:
            if assembly not in self.targetDatabases:
                self.targetDatabases[assembly] = []
            pep = f'{self.translationFolder}/{assembly}/{assembly}.pep'
            predicted = SeqIO.parse(pep, 'fasta')
            for record in predicted:
                self.targetDatabases[assembly].append(f'>{str(record.description)}\n{str(record.seq)}\n')
            reference = SeqIO.parse(self.proteome, 'fasta')
            for record in reference:
                self.targetDatabases[assembly].append(f'>{str(record.description).replace(" ", "_")}_ANNO\n{str(record.seq)}\n')

    def save_target_dbs(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        for assembly in self.targetDatabases:
            with open(f'{self.databaseDir}/{assembly}_target_database.fasta', 'w') as out:
                out.writelines(self.targetDatabases[assembly])

    def create_decoy_dbs(self):
        self.params.append(f'## {self.__class__.__name__}.{inspect.currentframe().f_code.co_name}')
        targets = os.listdir(self.databaseDir)
        for target in targets:
            decoy = Decoy(db=f'{self.databaseDir}/{target}')
            if self.args.cat:
                decoy.reverse_sequences().to_fasta(output=f'{self.databaseDir}/{target}'.replace("_target_", "_target_decoy_"),
                                                   pattern='rev', merge=True)
            # decoy.add_contaminants()
            else:
                decoy.reverse_sequences().to_fasta(output=f'{self.databaseDir}/{target}'.replace("_target_", "_decoy_"),
                                                   pattern='rev', merge=False)