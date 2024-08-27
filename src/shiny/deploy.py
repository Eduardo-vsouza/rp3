import os
import sys
from tokenize import group

from ..pipeline_config import PipelineStructure


class RPSDeployer:
    def __init__(self, outdir, shiny_r, paths_to_dfs, image_cols, image_dirs, df_alias, pdf_cols, pdf_dirs):
        self.pathsToDFs = paths_to_dfs
        self.imageCols = image_cols
        self.imageDirs = image_dirs
        self.dfAliases = df_alias
        self.pdfCols = pdf_cols
        self.pdfDirs = pdf_dirs

        # outdirs
        self.outdir = outdir
        self.shinyDir = f'{outdir}/RPS'
        self.shinyPDFDir = f'{self.shinyDir}/pdfs'
        self.shinyImageDir = f'{self.shinyDir}/images'
        self.shinyDFDir = f'{self.shinyDir}/data_frames'
        self.__check_dirs([self.outdir,
                           self.shinyDir,
                           self.shinyPDFDir,
                           self.shinyImageDir,
                           self.shinyDFDir])
        self.shinyPDFGroupDirs = ''
        self.shinyImageGroupDirs = ''
        self.shinyDFGroupDirs = ''

        # files
        self.shinyRScript = shiny_r
        self.customShinyRScript = f'{self.shinyDir}/app.R'


    def __check_dirs(self, folders):
        for folder in folders:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def copy_files(self):
        print(f"--Copying files to the ShinyApp folder.")

        pdf_outdirs = self.__generate_group_data(dirs=self.pdfDirs, app_outdir=self.shinyPDFDir, relative='pdfs')
        self.shinyPDFGroupDirs = pdf_outdirs
        print(self.shinyPDFGroupDirs)
        image_outdirs = self.__generate_group_data(dirs=self.imageDirs, app_outdir=self.shinyImageDir, relative='images')
        self.shinyImageGroupDirs = image_outdirs
        self.__get_data_frames()

    def __get_data_frames(self):
        files = self.pathsToDFs.split(',')
        for file in files:
            paths = file.split("/")
            if 'foldChange' in file:
                group = paths[-2]
                outdir = f'{self.shinyDFDir}/{group}'
                self.__check_dirs([outdir])
                df = paths[-1]
                app_path = f'{group}/{df}'
            else:
                app_path = f'{paths[-1]}'
            cmd = f'cp {file} {self.shinyDFDir}/{app_path}'
            os.system(cmd)
            self.shinyDFGroupDirs += f'data_frames/{app_path},'
        print(self.shinyDFGroupDirs)


    def __generate_group_data(self, dirs, app_outdir, relative):
        outdirs = ''
        folders = dirs.split(",")
        for folder in folders:
            groupdir = folder.split('/')[-4]
            outdir = f'{app_outdir}/{groupdir}'

            outdirs += f'{relative}/{groupdir},'
            self.__check_dirs([outdir])
            cmd = f'cp {folder}/* {outdir}/.'
            os.system(cmd)
        return outdirs[:-1]

    def edit_rscript(self):
        print(f"--Generating Shiny Script")
        outfile = []
        with open(self.shinyRScript, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if line.startswith("args"):
                    line = f'# {line}'
                if '(args' in line:
                    cols = line.split(" <- ")
                    end = ''
                    if 'args[1]' in cols[1]:
                        end = self.__get_value(self.shinyDFGroupDirs)
                    elif 'args[2]' in cols[1]:
                        end = self.__get_value(self.imageCols)
                    elif 'args[3]' in cols[1]:
                        end = self.__get_value(self.shinyImageGroupDirs)
                    elif 'args[4]' in cols[1]:
                        end = self.__get_value(self.dfAliases)
                    elif 'args[5]' in cols[1]:
                        end = self.__get_value(self.pdfCols)
                    elif 'args[6]' in cols[1]:
                        end = self.__get_value(self.shinyPDFGroupDirs)
                    line = f'{cols[0]} <- {end}'
                outfile.append(line)
        with open(self.customShinyRScript, 'w') as handler:
            handler.writelines(outfile)

    def __get_value(self, arg):
        value = f'strsplit(\"{arg}\", \",\")[[1]]\n'
        return value



