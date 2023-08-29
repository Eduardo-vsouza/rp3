import os
import sys
import inspect


class GTFtoFasta:
    def __init__(self, folder, genome, local_outdir):
        self.folder = folder
        self.genome = genome
        self.localOutdir = local_outdir

        self.__check_dirs()

        self.GTFFiles = os.listdir(self.folder)

        # params
        self.mode = 'ribocov'
        self.params = []

    def __check_dirs(self):
        if not os.path.exists(self.localOutdir):
            os.mkdir(self.localOutdir)

    def translate(self):
        genome_short_path = self.genome.split("/")[-1]
        cmd = f'cp {self.genome} {self.localOutdir}/.'
        self.params.append(cmd)
        os.system(cmd)
        for file in self.GTFFiles:

            local_dir = f'{self.localOutdir}/{file[:-4]}'
            os.mkdir(local_dir)

            cmd = f'cp {file} {self.localOutdir}/{local_dir}/.'
            self.params.append(cmd)
            os.system(cmd)

            cmd_gtf = f'{sys.path[0]}/dependencies/GTFtoFasta/GTFtoFasta {local_dir}/{file} {self.localOutdir}/{genome_short_path} met'
            self.params.append(cmd_gtf)
            os.system(cmd_gtf)

            cmd_gzip = f'gzip -d {self.localOutdir}/{local_dir}/*gz'
            self.params.append(cmd_gzip)
            os.system(cmd_gzip)
