import os
import sys
import json

from Bio import SeqIO

from ..pipeline_config import PipelineStructure


class RpHub(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.print_row(word='RpHub')
        print('\n')

        self.rpHubDir = f'{self.args.rpHubDir}'
        self.resultsHub = f'{self.rpHubDir}/results_hub'
        self.__setup_rphub()

        self.projectsFile = f'{self.rpHubDir}/projects.json'
        self.mpsFile = f'{self.rpHubDir}/mps.json'
        self.projects = {}
        self.__check_project()

    def __setup_rphub(self):
        self.check_dirs([self.rpHubDir, self.resultsHub])

    def __check_project(self):
        project = self.args.project
        projectdir = f'{self.resultsHub}/{project}'
        self.check_dirs([projectdir])
        if os.path.exists(self.projectsFile):
            # Opening JSON file
            with open(self.projectsFile, 'r') as handler:
                # Reading from json file
                self.projects = json.load(handler)
        else:
            self.projects[project] = {}


    def update_project(self, result):
        items = ['paths', 'runs']
        if self.args.project not in self.projects:
            self.projects[self.args.project] = {}
        for item in items:

            if item not in self.projects[self.args.project]:
                self.projects[self.args.project][item] = []

        if os.path.abspath(result) not in self.projects[self.args.project]['paths']:
            self.projects[self.args.project]['paths'].append(os.path.abspath(result))
        if result not in self.projects[self.args.project]['runs']:
            self.projects[self.args.project]['runs'].append(result)
        # print(self.projects)


    def integrate_results(self):
        self.print_row(word="Integrating results")
        for result in self.args.results:
            print(f"--Integrating results from {result} into RpHub.")
            outdir = f'{self.resultsHub}/{self.args.project}/{result}'
            self.outdir = result
            self.check_dirs([outdir])
            fasta = self.retrieve_file(result, filetype='fasta')
            # print(fasta)
            gtf = self.retrieve_file(result, filetype='gtf')
            peptides = self.retrieve_file(result, filetype='peptide')
            args = f'{result}/args.txt'
            self.update_project(result)

            files = [fasta, gtf, peptides, args]
            for file in files:
                if file != '':
                    if os.path.exists(file):
                        cmd = f'cp {file} {outdir}/.'
                        os.system(cmd)
            print(f"--Done integrating results into RpHub at {outdir}.\n")

    def retrieve_file(self, outdir, filetype='fasta'):
        file = ''
        if filetype == 'fasta':

            not_rescored = f'{outdir}/summarized_results/merged/microproteins_150.fasta'
            rescored = f'{outdir}/rescore/summarized_results/filtered_rescored_microproteins_150.fasta'
            if os.path.exists(rescored):
                file = rescored
            else:
                if os.path.exists(not_rescored):
                    file = not_rescored

        elif filetype == 'peptide':
            rescored = f'{outdir}/rescore/post_processing/group/peptides_fixed.txt'
            not_rescored = f'{outdir}/post_processing/group/db/peptides_fixed.txt'
            if os.path.exists(rescored):
                file = rescored
            else:
                if os.path.exists(not_rescored):
                    file = not_rescored

        elif filetype == 'gtf':
            rescored = f'{outdir}/rescore/summarized_results/filtered_rescored_microproteins_150.gtf'
            not_rescored = f'{outdir}/summarized_results/merged/microproteins_150.gtf'
            if os.path.exists(rescored):
                file = rescored
            else:
                if os.path.exists(not_rescored):
                    file = not_rescored

        return file


    def save(self):
        self.print_row()
        with open(self.projectsFile, "w") as outfile:
            json.dump(self.projects, outfile)
        print(f"--Saving Rp3 projects to {self.projectsFile}.")

    def generate_summary(self):
        project_num = len(self.projects)
        self.print_row(word="Projects", character="_")
        print('\n')

        print(f"--Number of projects: {project_num}")
        projects = ''
        runs = {}
        for project in self.projects:
            projects += f'{project}\n'
            runs[project] = len(self.projects[project]['runs'])

        # print(projects)
        for project in runs:
            print(f"{project}: {runs[project]} runs")

        self.count_microproteins()

    def count_microproteins(self):
        print('\n')

        self.print_row(word="Microproteins", character="_")
        print('\n')

        mps = {}
        for project in self.projects:
            if project not in mps:
                mps[project] = {}
                # mps[project] = {'runs': [], 'mps': []}
            for i, run in enumerate(self.projects[project]['paths']):
                name = self.projects[project]['runs'][i]
                if name not in mps[project]:
                    mps[project][name] = {}
                    mps[project][name]['mps'] = []
                    mps[project][name]['entries'] = {}

                fasta = self.retrieve_file(outdir=run, filetype='fasta')
                if os.path.exists(fasta):
                    if fasta != '':
                        records = SeqIO.parse(fasta, 'fasta')
                        for record in records:
                            seq = str(record.seq)
                            entry = str(record.description)
                            if entry not in mps[project][name]['entries']:
                                mps[project][name]['entries'][entry] = {'seq': seq, 'peptides': []}
                            if seq not in mps[project][name]['mps']:
                                mps[project][name]['mps'].append(seq)

        with open(self.mpsFile, "w") as outfile:
            json.dump(mps, outfile)
        # print(f"--Saving Rp3 projects to {}.")
        total_overall = []
        # for project in na
        for project in mps:
            total = []
            # print(f"{project}")
            self.print_row(word=project, character='-')
            # print(mps[project])
            for run in mps[project]:
                # print(mps[project][run])
                print(run, ': ',len(mps[project][run]['mps']))
                for mp in mps[project][run]['mps']:
                    if mp not in total_overall:
                        total_overall.append(mp)
                    if mp not in total:
                        total.append(mp)
                # total += len(mps[project][run]['mps'])
                # for mp in mps[project][run]['mps']:
                # print(project, run)
                # print(mps[project][run])
                # print(f"{run}: {len(mps[project][run]['mps'])}\n")

            print(f"--Total microproteins (nr): {len(set(total))}")
            print('\n')
        self.print_row(word="General summary", character='-')
        print(f"--Total microproteins across all projects (nr): {len(set(total_overall))}")


    def fetch_protein_seq(self):
        if self.args.proteinSeq is not None:
            protein = self.__fetch_protein(self.args.proteinSeq)
        

    def __fetch_protein(self, protein):
        found = False
        if os.path.exists(self.mpsFile):
            with open(self.mpsFile, 'r') as handler:
                # Reading from json file
                self.hubMps = json.load(handler)
            for project in self.hubMps:
                for run in self.hubMps[project]:
                    if protein in self.hubMps[project][run]['mps']:
                        self.print_row(word=run, character='*')
                        print(protein)
                        found = True
        if not found:
            print(f"{self.args.proteinSeq} not found.")
