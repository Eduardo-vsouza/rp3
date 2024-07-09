import os

import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D


from ..pipeline_config import PipelineStructure


class Peptide:
    def __init__(self, start, seq, mod_pep, strand):
        self.start = start
        self.sequence = seq
        self.strand = strand
        self.modPep = mod_pep
        self.__define_end()
    def __define_end(self):
        if self.strand == '+':
            self.end = self.start + ((len(self.sequence)) * 3)
        else:
            self.end = self.start + ((len(self.sequence))*3)


class GTFInfo:
    def __init__(self, cols):
        self.cols = cols
        self.chrom = None
        self.feature = None
        self.start = None
        self.end = None
        self.transcript = None
        self.strand = None
        self.gene = None
        self.__get_info()

    def __get_info(self):
        cols = self.cols
        self.chrom = cols[0]
        self.feature = cols[2]
        self.start, self.end = int(cols[3]), int(cols[4])
        self.strand = f'{cols[6]}'
        if cols[6] == '+':
            self.strand = +1
        else:
            self.strand = -1
        attrs = cols[8].split(";")
        for a in attrs:
            if 'transcript_id' in a:
                if self.feature == 'exon' and cols[1] != 'GTF2FASTA':

                    transcript = a.split(" ")[-1].replace("\"", "")
                    self.transcript = f'{transcript}_transcript'
                else:
                    self.transcript = a.split(" ")[-1].replace("\"", "")
            if 'gene_name' in a:
                self.gene = a.split(" ")[-1].replace("\"", "")

class PGContext(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.gtf = self.select_gtf()
        self.gtf = self.rescoredMicroproteinsGTF

        self.pgContextDir = f'{self.outdir}/pg_context'
        self.intermediatePGCFiles = f'{self.pgContextDir}/intermediate_files'
        self.contextFiguresDir = f'{self.pgContextDir}/context_figures'
        self.check_dirs([self.pgContextDir, self.intermediatePGCFiles, self.contextFiguresDir])

        self.transcriptsGTF = f'{self.intermediatePGCFiles}/transcripts_only.gtf'
        self.expandedGTF = f'{self.intermediatePGCFiles}/expanded_rp3_smorfs.gtf'
        self.overlappedGTF = f'{self.intermediatePGCFiles}/overlapped_smorfs.gtf'

        self.overlappedSmorfs = {}
        self.smorfInfo = {}
        self.microproteinSequences = {}
        self.smorfLimits = {}



        self.mainFeaturesToInclude = ['CDS', 'exon']

        self.mainColor = '#ffd700'
        self.smorfColor = "#ccccff"
        self.peptideColor = "#cffccc"

        self.currentLevel = 10
        self.levelIncrease = 15

        self.currentIsoform = 1
        self.currentTranscripts = {}
        self.currentTranscriptLimits = {}
        self.addedName = {}

    def expand_genes(self):
        print(f"--Expanding smORF coordinates")
        cmd = f"grep '	transcript' {self.gtf} > {self.transcriptsGTF}"
        os.system(cmd)
        cmd = (f'bedtools slop -i {self.transcriptsGTF} -g {self.args.chromSizes} -b {self.args.neighLength} -s > '
               f'{self.expandedGTF}')
        os.system(cmd)

    def intersect(self):
        print(f"--Intersecting smORFs with provided GTF file")
        cmd = (f'{self.toolPaths["bedtools"]} intersect -wao -a {self.expandedGTF} -b {self.args.gtf} > '
               f'{self.overlappedGTF}')
        os.system(cmd)

    def gather_overlaps(self):
        with open(self.overlappedGTF, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split('\t')
                smorf = cols[:9]
                # print(smorf)
                main = cols[9:]
                # print(main)
                smorf_info = GTFInfo(smorf)
                main_info = GTFInfo(main)
                if main_info.feature in self.mainFeaturesToInclude:
                    if smorf_info.transcript not in self.overlappedSmorfs:
                        self.overlappedSmorfs[smorf_info.transcript] = []
                    self.overlappedSmorfs[smorf_info.transcript].append(main_info)
        print("smorf", self.overlappedSmorfs)

    def gather_microproteins_data(self):
        with open(self.gtf, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split('\t')
                info = GTFInfo(cols)
                if info.transcript not in self.smorfInfo:
                    self.smorfInfo[info.transcript] = []
                if info.feature == 'exon':
                    self.smorfInfo[info.transcript].append(info)

    def gather_ms_peptides(self):
        import re
        peptides = self.select_peptides_df()
        print(peptides)
        df = pd.read_csv(peptides, sep='\t')
        df = df[df["q-value"] != 'q-value']
        df = df[df["q-value"] <= 0.01]
        proteins, peptides = df["proteinIds"].tolist(), df["peptide"].tolist()
        for protlist, peptide in zip(proteins, peptides):
            splat = protlist.split(",")
            fixed_pep= re.sub(r'[^a-zA-Z]', '', peptide)


            if ',' in protlist:
                utp = False
            else:
                utp = True
            for protein in splat:
                if protein in self.microproteinSequences and protein in self.smorfLimits:
                    prot_coord = self.microproteinSequences[protein]['sequence'].find(fixed_pep)

                    self.microproteinSequences[protein]['peptides'].append(fixed_pep)
                    self.microproteinSequences[protein]['utp'].append(utp)
                    if '+chr' in protein:
                        self.microproteinSequences[protein]['pep_coords'].append((len(fixed_pep)*3) + self.smorfLimits[protein]['start'])
                    else:
                        self.microproteinSequences[protein]['pep_coords'].append(self.smorfLimits[protein]['end'] - (len(fixed_pep)*3))

                    self.microproteinSequences[protein]['mod_peptides'].append(peptide)
        # print(self.microproteinSequences)
        for protein in self.microproteinSequences:
            if len(self.microproteinSequences[protein]['peptides']) > 0:
                for pep in self.microproteinSequences[protein]['peptides']:
                    print(protein, pep, len(self.microproteinSequences[protein]['peptides']))


    def gather_microprotein_sequences(self):
        fasta = self.select_fasta()
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            if '_F:' in str(record.description):
                self.microproteinSequences[str(record.description)] = {'sequence': str(record.seq),
                                                                       'peptides': [],
                                                                       'utp': [],
                                                                       'pep_coords': [],
                                                                       'mod_peptides': []}


    # def analyze_context(self):
    #     for smorf in self.overlappedSmorfs:
    #         start, end = self.__define_coordinates(smorf)
    #
    #         features = []
    #         smorf_rf = self.determine_reading_frame(self.smorfInfo[smorf][0].start)
    #         for main_feature in self.overlappedSmorfs[smorf]:
    #             main_orf_name = main_feature.transcript
    #             if main_feature.gene is not None:
    #                 main_orf_name = f'{main_feature.gene}'
    #             features.append(CustomGraphicFeature(start=main_feature.start, end=main_feature.end,
    #                                            strand=main_feature.strand, color=self.mainColor,
    #                                            label=main_orf_name, height=0.1))
    #             print(main_feature.strand)
    #             print(main_feature.start, main_feature.end, main_feature.strand, main_feature.feature, main_feature.transcript)
    #         for feature in self.smorfInfo[smorf]:
    #             features.append(GraphicFeature(start=feature.start, end=feature.end, strand=feature.strand,
    #                                            color=self.smorfColor, label=f'smORF RF{smorf_rf}', height=0.8))
    #
    #
    #         record = GraphicRecord(sequence_length=1000000000, features=features)
    #
    #
    #         # ax, _ = record.plot(figure_width=5)
    #
    #         cropped_record = record.crop((start, end))
    #         cropped_record.plot()
    #
    #         plt.show()
    def define_smorf_limits(self):
        for smorf in self.overlappedSmorfs:
            if smorf not in self.smorfLimits:
                self.smorfLimits[smorf] = {'start': None, 'end': None}

            for feature in self.smorfInfo[smorf]:
                if self.smorfLimits[smorf]['start'] is None:
                    self.smorfLimits[smorf]['start'] = feature.start
                    self.smorfLimits[smorf]['end'] = feature.end
                else:
                    if feature.start < self.smorfLimits[smorf]['start']:
                        self.smorfLimits[smorf]['start'] = feature.start
                    if feature.end > self.smorfLimits[smorf]['end']:
                        self.smorfLimits[smorf]['end'] = feature.end


    def analyze_context(self):
        print(len(self.overlappedSmorfs))
        for smorf in self.overlappedSmorfs:
            print(smorf)
            self.currentLevel = 10
            self.currentIsoform = 1
            self.addedName = {}
            self.currentTranscripts = {}
            self.currentTranscriptLimits = {}
            plt.clf()
            start, end = self.__define_coordinates(smorf)
            self.fig, self.ax = plt.subplots()
            features = []
            smorf_rf = self.determine_reading_frame(self.smorfInfo[smorf][0].start)
            for main_feature in self.overlappedSmorfs[smorf]:
                main_orf_name = main_feature.transcript
                level = self.__define_isoform(transcript=main_feature.transcript)

                if main_feature.gene is not None:
                    gene = f'{main_feature.gene}'
                if main_feature.feature == 'exon':
                    self.__add_feature(feature=main_feature, name=main_orf_name, color="lightgreen", placement=self.currentLevel)
                else:
                    self.__add_feature(feature=main_feature, name=main_feature.transcript, color=self.mainColor, placement=level)

                self.__define_isoform_limits(transcript=main_feature.transcript, start=main_feature.start,
                                             end=main_feature.end, strand=main_feature.strand, feature=main_feature.feature)
                feature = self.currentTranscriptLimits[main_feature.transcript]['feature']

            for name in self.addedName:
                rf = self.determine_reading_frame(self.currentTranscriptLimits[name]["start"])
                if self.currentTranscriptLimits[name]['feature'] != 'exon':
                    rf_name = f'{name} | {gene} | RF{rf}'
                else:
                    rf_name = f'{name} | {gene} | Transcript'
                self.__add_text(name, rf_name)


            self.__add_intron_line()

            self.smorf_start, self.smorf_end = None, None
            self.currentLevel += self.levelIncrease
            for feature in self.smorfInfo[smorf]:
                self.__define_smorf_limits(feature)
                self.__add_feature(feature=feature, name='', placement=self.currentLevel, color=self.smorfColor)
                # width = feature.end - feature.start
                # rectangle = patches.Rectangle((feature.start, 30), width, 10, linewidth=1, edgecolor='black', facecolor=self.smorfColor)
                # ax.add_patch(rectangle)
            self.ax.text(self.smorf_start+(self.smorf_end-self.smorf_start)/2, self.currentLevel - 5, f"smORF | RF{smorf_rf}", color="black", fontsize=10, ha='center', va='center')

            self.__add_intron_line(smorf=True, start=self.smorf_start, end=self.smorf_end, placement=self.currentLevel,
                                   strand=feature.strand, rf=smorf_rf)

            x_increase = 200
            self.__integrate_peptide_data(smorf, start-x_increase, end+x_increase)

            self.ax.set_xlim(start-x_increase, end+x_increase)
            self.ax.set_ylim(0, self.currentLevel + 10)
            self.ax.yaxis.set_visible(False)

            self.ax.spines['top'].set_visible(False)
            self.ax.spines['right'].set_visible(False)
            self.ax.spines['left'].set_visible(False)
            self.fig.set_size_inches(18.5, 10.5)

            plt.title(f'{smorf}_RF_{smorf_rf}')
            plt.tight_layout()
            plt.savefig(f'{self.contextFiguresDir}/{smorf}_context.png')
            # plt.show()


    def __integrate_peptide_data(self, smorf, start, end):
        plt.plot([start, end], [self.currentLevel+10, self.currentLevel+10], color='black')
        for i, peptide in enumerate(self.microproteinSequences[smorf]['peptides']):
            start = self.microproteinSequences[smorf]['pep_coords'][i]
            seq = peptide
            if '+chr' in smorf:
                strand = '+'
            else:
                strand = '-'
            feature = Peptide(start=start, seq=seq, mod_pep=self.microproteinSequences[smorf]['mod_peptides'][i],
                              strand=strand)
            self.currentLevel += 10
            self.__add_feature(feature=feature, name=seq, placement=self.currentLevel, height=2, color=self.smorfColor)
            self.ax.text(start+(feature.end - feature.start)/2, self.currentLevel - 2.5, feature.modPep, color="black", fontsize=10, ha='center', va='center')

    def __define_smorf_limits(self, feature):
        if self.smorf_start is None:
            self.smorf_start, self.smorf_end = feature.start, feature.end
        else:
            if feature.start < self.smorf_start:
                self.smorf_start = feature.start
            if feature.end > self.smorf_end:
                self.smorf_end = feature.end
        self.smorfLimits[feature.transcript] = self.smorf_start

    def __add_feature(self, feature, name, color, placement=10, height=5):
        height = height
        width = feature.end - feature.start
        x, y = feature.start, placement
        bbox_props = dict(boxstyle="larrow,pad=0.0", edgecolor='black', facecolor=color, linewidth=1)
        # rectangle = patches.FancyBboxPatch((x, y), width, height, **bbox_props)

        rectangle = patches.Rectangle((x, y), width, height, linewidth=1, edgecolor='black', facecolor=color)
        arrow_width = 0.05 * width  # Adjust this value to make the arrow more pronounced
        # arrow = patches.FancyArrow(x + width, y + height / 2, arrow_width, 0, width=height,
        #                            length_includes_head=True, head_width=height, head_length=arrow_width,
        #                            edgecolor='black',
        #                            color=color)
        self.ax.add_patch(rectangle)


        if name not in self.addedName:
            self.addedName[name] = (feature.start, placement)
        else:
            if feature.start < self.addedName[name][0]:
                self.addedName[name] = (feature.start, placement)


        # features.append(CustomGraphicFeature(start=main_feature.start, end=main_feature.end,
        #                                strand=main_feature.strand, color=self.mainColor,
        #                                label=main_orf_name, height=0.1))
        self.ax.add_patch(rectangle)

    def __add_text(self, name, rf_name):
        # x, y = feature.start, placement
        x, y = self.addedName[name][0], self.addedName[name][1]
        text_x = x
        text_y = y

        self.ax.text(text_x, text_y - 2.5, rf_name, color="black", fontsize=10, ha='center', va='center')

    def __add_intron_line(self, rf=None, smorf=False, start=None, end=None, placement=None, strand=None):
        colors = {1: 'blue', 2: 'red', 3: 'green'}
        if not smorf:
            for transcript in self.currentTranscriptLimits:
                rf = self.determine_reading_frame(self.currentTranscriptLimits[transcript]["start"])
                marker = '>'
                if self.currentTranscriptLimits[transcript]['strand'] == -1:
                    marker = '<'
                # print(self.currentTranscriptLimits[transcript]['end'])
                start, end = self.currentTranscriptLimits[transcript]['start'], self.currentTranscriptLimits[transcript]['end']
                x = [start, end]
                y = [self.currentTranscripts[transcript]+2.5, self.currentTranscripts[transcript]+2.5]
                line_style = (0, (1, 1))  # Custom dash pattern: 1 point on, 1 point off
                if self.currentTranscriptLimits[transcript]['feature'] != 'exon':
                    color = colors[rf]
                else:
                    color = 'black'
                line = Line2D(x, y, marker=marker, color=color)
                self.ax.add_line(line)
                # plt.plot(x, y)
        else:
            if strand == -1:
                marker = '<'
            else:
                marker = '>'
            x = [start, end]
            y = [placement+2.5, placement+2.5]
            line = Line2D(x, y, marker=marker, color=colors[rf])
            self.ax.add_line(line)


    def __define_isoform_limits(self, transcript, start, end, strand, feature):
        # print(feature)
        if transcript not in self.currentTranscriptLimits:
            self.currentTranscriptLimits[transcript] = {'start': start, 'end': end, 'strand': strand, 'feature': feature}
        else:
            if self.currentTranscriptLimits[transcript]['start'] > start:
                self.currentTranscriptLimits[transcript]['start'] = start
            if self.currentTranscriptLimits[transcript]['end'] < end:
                self.currentTranscriptLimits[transcript]['end'] = end

    def __define_isoform(self, transcript):

        # print(transcript)
        # if '.' in transcript:
        #     isoform = int(transcript.split(".")[1])
        # else:
        #     isoform = 1
        if transcript not in self.currentTranscripts:
            self.currentLevel += 10
            self.currentTranscripts[transcript] = self.currentLevel
        return self.currentTranscripts[transcript]


    @staticmethod
    def determine_reading_frame(start):
        frame = (start % 3) + 1
        return frame

    def __define_coordinates(self, smorf):
        start = self.smorfInfo[smorf][0].start
        end = self.smorfInfo[smorf][-1].end
        for feature in self.overlappedSmorfs[smorf]:
            if feature.start < start:
                start = feature.start
            if feature.end > end:
                end = feature.end
        return start, end






