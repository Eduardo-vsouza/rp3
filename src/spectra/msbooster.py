import os
import sys

import pandas as pd

from ..pipeline_config import PipelineStructure


class Booster(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.msBoosterParams = f'{self.boosterDir}/msbooster_params.txt'
        self.boosterTmpDir = f'{self.boosterDir}/tmp'
        self.check_dirs([self.boosterDir, self.boosterPinDir, self.boosterTmpDir])

        self.pinFolder, self.boosterSearchDir = self.__define_pin_folder()

    def __define_pin_folder(self):
        if self.args.rescore:
            pin_folder = f'{self.rescoreSearchDir}'
            search_dir = self.rescoreSearchDir
        else:
            pin_folder = None
            search_dir = f'{self.searchDir}/group'
            db_dirs = os.listdir(search_dir)
            for subdir in db_dirs:
                if subdir.endswith("target_decoy_database.fasta"):
                    pin_folder = f'{search_dir}/{subdir}'
                    break
        return pin_folder, search_dir

    def prepare_pin_files(self):
        print("Preparing pin files.")
        subdirs = os.listdir(self.boosterSearchDir)
        for subdir in subdirs:
            if not self.args.rescore:
                if subdir.endswith("target_decoy_database.fasta"):
                    pin_files = os.listdir(f'{self.boosterSearchDir}/{subdir}')
                    for file in pin_files:
                        pin = f'{self.boosterSearchDir}/{subdir}/{file}'
                        cmd = f'mv {pin} {pin.replace("_target.pin", ".pin")}'
                        os.system(cmd)
            else:
                pin_files = os.listdir(f'{self.boosterSearchDir}/{subdir}')
                for file in pin_files:
                    pin = f'{self.rescoreSearchDir}/{subdir}/{file}'
                    cmd = f'mv {pin} {pin.replace("_target.pin", ".pin")}'
                    os.system(cmd)
        # pin_files = {}

        # with open(f'{self.rescoreDir}/all_pins.pin', 'r') as handler:
        #     lines = handler.readlines()
        #     header = ''
        #     for line in lines:
        #         if 'SpecId' in line:
        #             header = line
        #         else:
        #             cols = line.split('\t')
        #             specs = cols[0]
        #             file = specs.split(".")[0]
        #             if file not in pin_files:
        #                 pin_files[file] = [header]
        #             pin_files[file].append(line)
        #
        # for file in pin_files:
        #     with open(f'{self.boosterPinDir}/{file}.pin', 'w') as outpin:
        #         outpin.writelines(pin_files[file])

    def configure_parameters(self):
        print("Configuring parameters for MSBooster.")
        params = []
        print(self.toolPaths)

        folder = '/'.join(self.toolPaths["MSBooster"].split("/")[:-1])
        mzmls = os.listdir(self.args.mzml)
        mzml_folders = ''
        for f in mzmls:
            mzml_folders += f'{self.args.mzml}/{f} '
        with open(f'{folder}/msbooster_params.txt', 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if 'numThreads' in line:
                    line = f'numThreads = {self.args.threads}\n'
                if 'mzmlDirectory =' in line:
                    line = f'mzmlDirectory = {mzml_folders[:-1]}\n'
                if 'pinPepXMLDirectory =' in line:
                    line = 'pinPepXMLDirectory ='
                    search_files = os.listdir(self.boosterSearchDir)
                    for subdir in search_files:
                        line += f' {self.boosterSearchDir}/{subdir}/'
                    line += '\n'
                params.append(line)
        with open(self.msBoosterParams, 'w') as outfile:
            outfile.writelines(params)

    def run(self):
        print(f"--Running MSBooster on your pin files before rescoring.")
        cmd = f'java -jar {self.toolPaths["MSBooster"]} --paramsList {self.msBoosterParams}'
        os.system(cmd)

    def merge_pin_files(self):
        cats = ''
        subdirs = os.listdir(self.boosterSearchDir)
        for subdir in subdirs:
            pins = os.listdir(f'{self.boosterSearchDir}/{subdir}')
            for pin in pins:
                if pin.endswith('_edited.pin'):
                    cats += f' {self.boosterSearchDir}/{subdir}/{pin}'
        cmd = f'cat {cats} > {self.boosterTmpDir}/merged_edited.pin'
        os.system(cmd)  # merge pin files
        cmd = f"awk 'FNR<2' {self.boosterTmpDir}/merged_edited.pin > {self.boosterTmpDir}/header.txt"
        os.system(cmd)  # get pin header
        cmd = f"grep -v 'SpecId' {self.boosterTmpDir}/merged_edited.pin > {self.boosterTmpDir}/no_header.pin"
        os.system(cmd)  # remove duplicated headers
        cmd = (f"cat {self.boosterTmpDir}/header.txt {self.boosterTmpDir}/no_header.pin > "
               f"{self.boosterTmpDir}/merged_pin_spectralEntropy.pin")
        os.system(cmd)  # concatenate header and remainder of pin file
        new_lines = []
        with open(f"{self.boosterTmpDir}/merged_pin_spectralEntropy.pin", 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split('\t')
                del cols[19]
                new_lines.append('\t'.join(cols))
        with open(self.mergedBoosterPin, 'w') as outpin:
            outpin.writelines(new_lines)