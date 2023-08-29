import os
import sys

from ..pipeline_config import PipelineStructure


class RiboSeqAlign(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.set_ribocov_params()
        self.adapterSequence = args.adapter

        self.index = f'{sys.path[0]}/data/STAR_indexes/hg19.star'
        self.indexCont = f'{sys.path[0]}/data/STAR_indexes/hg19_cont.star'
        self.__check_indexes()

    def trim_reads(self):
        files = os.listdir(self.args.fastq)
        for file in files:
            print(f'Performing read trimming for {file}.')
            full = f'{self.args.fastq}/{file}'
            cmd = (f'{self.args.fastx_clipper_path} -Q33 -l 20 -n -v -c -a {self.adapterSequence} -i {full} | '
                   f'{self.args.fastx_trimmer_path} -Q33 -f 1 2> {self.riboSeqTrimmedDir}/trim.log > '
                   f'{self.riboSeqTrimmedDir}/trimmed_{file}')
            os.system(cmd)
            print(f"Finished trimming {file}.")

    def __check_indexes(self):
        if self.args.genome_index is not None:
            if self.args.genome_index == 'hg38':
                self.index = f'{self.indexesDir}/hg38.star'
            else:
                self.index = self.args.genome_index
        else:
            print(f"Genome index not provided. Using default index at {self.index}.")

        if self.args.cont_index is not None:
            if self.args.cont_index == 'hg38':
                self.indexCont = f'{self.indexesDir}/hg38cont.star'
            else:
                self.indexCont = self.args.cont_index
        else:
            print(f"Contaminant index not provided. Using default index at {self.indexCont}.")

    def remove_ribosome(self):
        files = os.listdir(self.riboSeqTrimmedDir)
        for file in files:
            if file.endswith(".fastq"):
                print(f"Removing contaminants from {file}")
                cmd = (f'{self.toolPaths["STAR"]} --outSAMstrandField intronMotif --outReadsUnmapped Fastx --genomeDir '
                       f'{self.indexCont} --runThreadN {self.args.threads} --readFilesIn {self.riboSeqTrimmedDir}/{file} '
                       f'--outFileNamePrefix {self.riboSeqContaminantAlnDir}/no_contaminant_{file}')
                os.system(cmd)
                print(f"Finished removing contaminants from {file}")

    def align_rpfs(self):
        files = os.listdir(self.riboSeqContaminantAlnDir)
        for file in files:
            if 'Unmapped.out.mate1' in file:
                if file.endswith(".gz"):
                    gz = ' --readFilesCommand zcat'
                else:
                    gz = ''

                print(f"Aligning RPF reads from {file} to {self.index}")
                cmd = (f'{self.toolPaths["STAR"]} --outSAMstrandField intronMotif --genomeDir {self.index} --runThreadN '
                       f'{self.args.threads} --readFilesIn {self.riboSeqContaminantAlnDir}/{file} --outFileNamePrefix '
                       f'{self.riboSeqAlnDir}/aligned_to_genome_{file} --outFilterMismatchNmax 2 '
                       f'--outFilterMultimapNmax {self.args.multimappings} --chimScoreSeparation 10 --chimScoreMin '
                       f'20 --chimSegmentMin 15 --outSAMattributes All{gz}')
                os.system(cmd)
                print(f'Finished aligning {file}.')



