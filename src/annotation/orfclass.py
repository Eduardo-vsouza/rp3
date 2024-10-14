import os
import sys
import collections
from collections import OrderedDict

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import subprocess
import re

import pandas as pd

from ..pipeline_config import PipelineStructure


class ORFClassification(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.smorfsGTF = self.rescoredMicroproteinsGTF
        self.refGTF = self.args.gtf

        self.nonIntersectionFile = f'{self.orfClassDir}/non_intersection_annotation.gtf'
        self.intersectionFile = f'{self.orfClassDir}/intersection_annotation.gtf'
        self.annotationOutfile = f'{self.orfClassDir}/predicted_nonhomolog_smorfs_annotation'

    # @staticmethod
    # def check_dirs(folders):
    #     for folder in folders:
    #         if not os.path.exists(folder):
    #             os.mkdir(folder)

    def intersect(self):
        # Set the input and output file paths

        # Define the command for intersection
        intersection_command = ['bedtools', 'intersect', '-wo', '-f', '1', '-s', '-a', self.smorfsGTF, '-b', self.refGTF]
        # intersection_command = ['bedtools', 'intersect', '-wo', '-s', '-a', self.smorfsGTF, '-b', self.refGTF]

        # Execute the intersection command and redirect the output to a file
        with open(self.intersectionFile, 'w') as outfile:
            subprocess.run(intersection_command, stdout=outfile, check=True)

        # Define the command for non-intersection
        non_intersection_command = ['bedtools', 'intersect', '-v', '-s', '-a', self.smorfsGTF, '-b', self.refGTF]

        # Execute the non-intersection command and redirect the output to a file
        with open(self.nonIntersectionFile, 'w') as outfile:
            subprocess.run(non_intersection_command, stdout=outfile, check=True)

        print(f"Intersection output file '{self.intersectionFile}' created successfully.")
        print(f"Non-intersection output file '{self.nonIntersectionFile}' created successfully.")

    def annotate(self):

        gene_data = {}

        with open(self.intersectionFile, 'r') as file:
            for line in file:
                if line.startswith('chr'):
                    parts = line.strip().split('\t')
                    attributes = parts[8]
                    gene_id_match = re.search(r'gene_id "(.+?)"', attributes)
                    gene_id = gene_id_match.group(1) if gene_id_match else 'Unknown'

                    gene_info = parts[17]
                    gene_biotype_match = re.search(r'gene_biotype "([^"]+)"', gene_info)
                    gene_biotype = gene_biotype_match.group(1) if gene_biotype_match else 'Unnamed'

                    gene_name_match = re.search(r'gene_name "([^"]+)"', gene_info)
                    gene_name = gene_name_match.group(1) if gene_name_match else 'Unnamed'

                    # Remove quotes from gene name
                    gene_name = gene_name.replace('"', '')

                    transript_biotype_match = re.search(r'transcript_biotype "([^"]+)"', gene_info)
                    transcript_biotype = transript_biotype_match.group(1) if transript_biotype_match else 'Unknown'

                    annotation = 'UA'

                    # Let's check for the simple cases first
                    if gene_biotype == 'lncRNA' or gene_biotype == 'lincRNA' or gene_biotype == 'antisense' or gene_biotype == 'sense_intronic' or gene_biotype == 'sense_overlapping':
                        annotation = 'lncRNA'
                    elif parts[11] == 'five_prime_utr' or parts[11] == '5UTR':
                        annotation = 'uORF'
                    elif parts[11] == 'three_prime_utr' or parts[11] == '3UTR':
                        annotation = 'dORF'
                    elif parts[11] == 'exon':
                        annotation = 'eORF'
                    elif parts[11] == 'CDS':
                        annotation = 'oCDS'

                    # Now let's check for the more complicated cases
                    if 'ambiguous_orf' in transcript_biotype:
                        annotation = 'ndORF'
                    if 'nonsense_mediated_decay' in transcript_biotype:
                        annotation = 'ndORF'
                    if 'protein_coding_CDS_not_defined' in transcript_biotype:
                        annotation = 'ndORF'
                    if 'retained_intron' in transcript_biotype:
                        annotation = 'riORF'
                    if parts[11] == 'exon' and 'exon_number "1"' in gene_info:
                        annotation = 'aiORF'
                    if parts[
                        11] == 'exon' and 'exon_number "1"' in gene_info and 'protein_coding_CDS_not_defined' in transcript_biotype:
                        annotation = 'aindORF'
                    if 'processed_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'
                    if "transcribed_unitary_pseudogene" in transcript_biotype:
                        annotation = 'rtORF'
                    if 'unprocessed_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'
                    if 'transcribed_processed_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'
                    if 'transcribed_unprocessed_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'
                    if 'translated_unprocessed_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'
                    if 'translated_processed_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'
                    if 'polymorphic_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'
                    if 'unitary_pseudogene' in transcript_biotype:
                        annotation = 'rtORF'

                    # This order is arbitrary, but it's the order that I think is interesting...
                    # starting with "alternative initiation of a non-defined ORF (aindORF)"
                    # lncrna = long non-coding RNA
                    # aindORF = alternative initiation of a non-defined ORF
                    # uORF = upstream ORF
                    # aiORF = alternative initiation of a defined ORF
                    # ndORF = non-defined ORF
                    # rtORF = retrotransposed ORF
                    # riORF = retained intron ORF
                    # dORF = downstream ORF
                    # eORF = exon ORF
                    # oCDS = overlapping CDS
                    priority_order = ['rtORF', 'lncRNA', 'aindORF', 'uORF', 'aiORF', 'ndORF', 'riORF', 'dORF', 'eORF',
                                      'oCDS']

                    if gene_id not in gene_data:
                        gene_data[gene_id] = (annotation, gene_name)
                    else:
                        existing_annotation, _ = gene_data[gene_id]

                        if annotation in priority_order and existing_annotation not in priority_order:
                            gene_data[gene_id] = (annotation, gene_name)
                        elif annotation in priority_order and existing_annotation in priority_order:
                            if priority_order.index(annotation) < priority_order.index(existing_annotation):
                                gene_data[gene_id] = (annotation, gene_name)

            # How many gene_ids are there?
            print(len(gene_data))

        # Load non-intersection data and extract gene IDs
        non_intersection_gene_ids = set()

        with open(self.nonIntersectionFile, 'r') as file:
            for line in file:
                if line.startswith('chr'):
                    parts = line.strip().split('\t')
                    attributes = parts[8]
                    gene_id_match = re.search(r'gene_id "(.+?)"', attributes)
                    gene_id = gene_id_match.group(1) if gene_id_match else 'Unknown'
                    non_intersection_gene_ids.add(gene_id)

        # Delete duplicate gene IDs
        non_intersection_gene_ids = list(non_intersection_gene_ids)

        # Add gene IDs to gene_data
        for gene_id in non_intersection_gene_ids:

            if gene_id not in gene_data:
                gene_data[gene_id] = ('Intergenic', 'Intergenic')

        with open(self.annotationOutfile, 'w') as output:
            for gene_id, (annotation, gene_name) in gene_data.items():
                output.write(f'{gene_id}\t{annotation}\t{gene_name}\n')

        print(f"Output file '{self.annotationOutfile}' created successfully.")


class ORFClassVis(ORFClassification):
    def __init__(self, args):
        super().__init__(args=args)
        self.annotationDF = pd.read_csv(self.annotationOutfile, header=None,
                                        names=["smorf", "annotation", "gene"], sep='\t')

        self.annotations = {}

        self.smorfs = {}

        self.annotationByGroup = {}
        self.annotationPercentages = {}

    def classify_by_groups(self):
        df = pd.read_csv(self.microproteinMappingGroupsForPlotsUnion, sep='\t')
        smorfs, groups = df["smorf"].tolist(), df["group"].tolist()
        for smorf, group in zip(smorfs, groups):
            self.smorfs[smorf] = group

    def get_annotations(self):
        smorfs, annos = self.annotationDF["smorf"].tolist(), self.annotationDF["annotation"].tolist()
        for smorf, anno in zip(smorfs, annos):
            self.annotations[smorf] = anno
        print(self.annotations)

    def annotate_groups(self):
        for smorf in self.smorfs:
            groups = self.smorfs[smorf].split(",")
            for group in groups:
                if group not in self.annotationByGroup:
                    self.annotationByGroup[group] = []
                if smorf in self.annotations:
                    self.annotationByGroup[group].append(self.annotations[smorf])
        print(self.annotationByGroup)
        # groups = collections.Counter(self.annotationByGroup.keys())
        # print(groups)

    def plot_group_per_anno(self):
        annos = {}
        for group in self.annotationByGroup:
            # if group not in annos:
            #     annos[group] = {}
            for anno in self.annotationByGroup[group]:
                if anno not in annos:
                    annos[anno] = {}
                if group not in annos[anno]:
                    annos[anno][group] = 0
                annos[anno][group] += 1
                # annos[group][anno] = self.annotationByGroup[group]
        print(annos)
        norm_annos = {}
        for group in annos:
            if group not in norm_annos:
                norm_annos[group] = {}
            for anno in annos[group]:
                norm_annos[group][anno] = annos[group][anno]/len(annos[group])
        annos = norm_annos
        labels = set().union(*[v.keys() for v in annos.values()])

        # Extract the keys from the dictionary
        keys = list(annos.keys())

        # Initialize an empty dictionary for each label
        stacked_values = {label: [] for label in labels}

        # Populate the stacked_values dictionary
        for label in labels:
            for key in keys:
                value = annos[key].get(label, 0)
                stacked_values[label].append(value)

        # Plotting the stacked bar chart
        colors = ["#c45db9", "#7ab341", "#7b60cf", "#d09244", "#7e7fc5", "#87853a", "#c75980", "#51a876", "#ca5542",
                  "#45b0cf"]
        colors = ["#95c3ff",
                  "#ca6600",
                  "#016eee",
                  "#b7cc38",
                  "#e189ff",
                  "#007710",
                  "#7e3978",
                  "#a1cd97",
                  "#ff879f",
                  "#813e52"]
        bottom = None
        colors = sns.color_palette("husl", len(labels))

        for i, label in enumerate(labels):
            # plt.bar(keys, stacked_values[label], bottom=bottom, label=label, edgecolor='black',
            #         color=colors[i % len(colors)])
            plt.bar(keys, stacked_values[label], bottom=bottom, label=label, edgecolor='black',
                    color=colors[i % len(colors)])
            if bottom is None:
                bottom = stacked_values[label]
            else:
                bottom = [sum(x) for x in zip(bottom, stacked_values[label])]

        # Set the labels and title
        plt.xlabel('Groups')
        plt.ylabel('Annotation')
        # plt.title('Stacked Bar Plot')
        # ax.set_xticklabels(ax.get_xticklabels())
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # Add a legend
        # plt.legend()
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Update the ticks and tick labels
        ax.tick_params(top=False, right=False)
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # Show the plot
        plt.show()


    def get_annotation_percentages(self):
        totals = {}
        for group in self.annotationByGroup:
            if group not in self.annotationPercentages:
                self.annotationPercentages[group] = {}
                totals[group] = 0
            for anno in self.annotationByGroup[group]:
                totals[group] += 1

        for group in self.annotationByGroup:
            for anno in self.annotationByGroup[group]:
                if anno not in self.annotationPercentages[group]:
                    count = self.annotationByGroup[group].count(anno)
                    self.annotationPercentages[group][anno] = count/totals[group]
        print(self.annotationPercentages)

    def plot_stacked_bar(self):
        labels = set().union(*[v.keys() for v in self.annotationPercentages.values()])

        # Extract the keys from the dictionary
        keys = list(self.annotationPercentages.keys())

        # Initialize an empty dictionary for each label
        stacked_values = {label: [] for label in labels}
        # stacked_values = OrderedDict({label: [] for label in labels})

        # Populate the stacked_values dictionary
        for label in labels:
            for key in keys:
                value = self.annotationPercentages[key].get(label, 0)
                stacked_values[label].append(value)

        # Plotting the stacked bar chart
        colors = ["#c45db9", "#7ab341", "#7b60cf", "#d09244", "#7e7fc5", "#87853a", "#c75980", "#51a876", "#ca5542",
                  "#45b0cf"]
        colors = ["#003f5c", "#345176", "#5e638e", "#8775a3", "#b187b3", "#d99abf", "#ffb0c8"]
        palette = sns.color_palette("pastel", len(labels))
        colors = sns.color_palette("pastel", len(labels))

        # colors = ["#003f5c", "#034977", "#295090", '#5354a5', "#8053b2", "#ae4db7", "#db40b2"]
        # sns.color_palette("rocket")
        # palette = sns.color_palette("rocket", len(labels))
        # colors = ['#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd']
        bottom = None
        for i, label in enumerate(labels):
            plt.bar(keys, stacked_values[label], bottom=bottom, edgecolor='black', label=label, width=0.9,
                    color=colors[i])
            # plt.bar(keys, stacked_values[label], bottom=bottom, edgecolor='black', label=label, width=0.9)
            # plt.bar(keys, stacked_valReviewsues[label], bottom=bottom, edgecolor='black', label=label, width=0.9)

            if bottom is None:
                bottom = stacked_values[label]
            else:
                bottom = [sum(x) for x in zip(bottom, stacked_values[label])]

        # Set the labels and title
        plt.xlabel('Groups')
        plt.ylabel('Annotation')
        # plt.title('Stacked Bar Plot')
        # ax.set_xticklabels(ax.get_xticklabels())
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # Add a legend
        # plt.legend()
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # handles, labels = plt.gca().get_legend_handles_labels()
        # sorted_labels, sorted_handles = zip(*sorted(zip(labels, handles)))
        # plt.legend(sorted_handles, sorted_labels)

        # Update the ticks and tick labels
        ax.tick_params(top=False, right=False)
        handles, labels = plt.gca().get_legend_handles_labels()
        sorted_labels, sorted_handles = zip(*sorted(zip(labels, handles)))
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.xticks(rotation=45)

        plt.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1, 1))
        # plt.legend(sorted_handles, sorted_labels, loc='upper left', bbox_to_anchor=(1, 1))

        # plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # Show the plot
        plt.show()
