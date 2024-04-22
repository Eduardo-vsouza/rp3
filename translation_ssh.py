import os
import sys
import inspect


class GTFtoFasta:
    def __init__(self, folder, genome, local_outdir):
        self.folder = folder
        self.genome = genome
        self.localOutdir = local_outdir

        self.__check_dirs()

        # self.GTFFiles = os.listdir(self.folder)
        self.GTFFiles = self.folder
        # params
        self.mode = 'ribocov'
        self.params = []

    def __check_dirs(self):
        if not os.path.exists(self.localOutdir):
            os.mkdir(self.localOutdir)

    def translate(self):
        genome_short_path = self.genome.split("/")[-1]
        # cmd = f'cp {self.genome} {self.localOutdir}/.'
        # self.params.append(cmd)
        # os.system(cmd)
        if os.path.isdir(self.GTFFiles):
            gtf_files = os.listdir(self.GTFFiles)
            for file in gtf_files:
                if file.endswith("gtf"):

                    local_dir = f'{self.localOutdir}/{file[:-4]}'
                    outfile = f'{local_dir}/{file[:-4]}_ORFs.gtf'
                    print(outfile)

                    if not os.path.exists(outfile):

                        if not os.path.exists(local_dir):
                            os.mkdir(local_dir)

                        cmd = f'cp {self.GTFFiles}/{file} {local_dir}/.'
                        self.params.append(cmd)
                        os.system(cmd)
                        # self.__check_dirs([])
                        cmd_gtf = f'{sys.path[0]}/dependencies/GTFtoFasta/GTFtoFasta {local_dir}/{file} {self.genome} met'
                        self.params.append(cmd_gtf)
                        os.system(cmd_gtf)

                        cmd_gzip = f'gzip -d {local_dir}/*.gz'
                        self.params.append(cmd_gzip)
                        os.system(cmd_gzip)
        else:
            file = self.GTFFiles.split("/")[-1]
            fullfile = self.GTFFiles
            local_dir = f'{self.localOutdir}/{file[:-4]}'
            if not os.path.exists(local_dir):
                os.mkdir(local_dir)
            outfile = f'{local_dir}/{file[:-4]}_ORFs.gtf'
            print(outfile)
            if not os.path.exists(outfile):
                cmd = f'cp {fullfile} {local_dir}/.'
                self.params.append(cmd)
                os.system(cmd)
                # self.__check_dirs([])
                print(f' {local_dir}/{file} ')
                # os.system(f'cat {local_dir}/{file}')
                cmd_gtf = f'{sys.path[0]}/dependencies/GTFtoFasta/GTFtoFasta {local_dir}/{file} {self.genome} met'
                self.params.append(cmd_gtf)
                os.system(cmd_gtf)

                cmd_gzip = f'gzip -d {local_dir}/*.gz'
                self.params.append(cmd_gzip)
                os.system(cmd_gzip)
