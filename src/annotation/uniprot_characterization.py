import os
import sys
import re

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm import tqdm

from ..pipeline_config import PipelineStructure


class UniprotAnnotation(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)
        self.annotatedDir = f'{self.outdir}/annotated_microproteins'
        self.annotatedPlotsDir = f'{self.annotatedDir}/plots'
        self.annotatedCountsDir = f'{self.annotatedDir}/counts'
        self.annotatedIntersectedFile = f'{self.annotatedDir}/intersected_annotation_uniprot.txt'

        self.check_dirs([self.annotatedDir, self.annotatedPlotsDir, self.annotatedCountsDir])
        self.annotation = {'seq': [], 'protein_evidence': [], 'reviewed': [], 'gene': [], 'annotation': []}
        self.filteredAnnotation = {}
        self.microproteinsSharedPeptides = {}
        self.identifiedUncharacterized = {}

        self.uniprotAnnotationMPIntersection = f'{self.annotatedCountsDir}/microproteins_sharing_peptides_with_annotated_uniprot_by_annotationLevel.txt'
        self.uniprotEvidenceMPIntersection = f'{self.annotatedCountsDir}/microproteins_sharing_peptides_with_annotated_uniprot_by_proteinEvidence.txt'

        self.uniprotAnnotationAll = f'{self.annotatedCountsDir}/annotated_uniprot_identified_by_annotationLevel.txt'
        self.uniprotEvidenceAll = f'{self.annotatedCountsDir}/annotated_uniprot_identified_by_proteinEvidence.txt'

        self.uniprotAnnotationMP = f'{self.annotatedCountsDir}/annotated_uniprot_microproteins_identified_by_annotationLevel.txt'
        self.uniprotEvidenceMP = f'{self.annotatedCountsDir}/annotated_uniprot_microproteins_identified_by_proteinEvidence.txt'

    def get_annotations(self, peps):
        df = pd.read_csv(self.args.uniprotTable, sep='\t')
        # reviewed = df["Reviewed"].tolist()
        # genes = df["Gene Names"].tolist()
        # seqs = df["Sequence"].tolist()
        # annotation = df["Annotation"].tolist()
        pattern = '|'.join(map(re.escape, peps))
        df = df[df["Sequence"].str.contains(pattern)]
        df.to_csv(self.annotatedIntersectedFile, sep='\t', index=False)

    #
    #
    # def filter_sequence(self, seq, peptides):
    #     if any(pep in seq for pep in peptides):
    #         return seq
    #
    # def parallel_filter_annotation(self, peptides):
    #     import concurrent.futures
    #
    #     filtered_annotation = {}
    #     with concurrent.futures.ThreadPoolExecutor() as executor:
    #         futures = [executor.submit(self.filter_sequence, seq, peptides) for seq in tqdm(self.annotation)]
    #         for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
    #             result = future.result()
    #             if result is not None:
    #                 filtered_annotation[result] = self.annotation[result]
    #     self.filteredAnnotation = filtered_annotation

    def filter_annotations(self):

        import concurrent.futures

        df = pd.read_csv(f'{self.rescorePostProcessDir}/group/peptides_fixed.txt', sep='\t')
        df = df[df["q-value"] <= 0.01]
        df = df[df["proteinIds"].str.contains("rev_") == False]
        peptides_unfixed = df["peptide"].tolist()
        peptides = []
        for pep in peptides_unfixed:
            letters_only = re.findall('[a-zA-Z]', pep)
            peptide = ''.join(letters_only)
            peptides.append(peptide)
        # peptides = list(set(peptides))
        # proteins = df["proteinIds"].tolist()
        self.get_annotations(peps=list(set(peptides)))
        # self.parallel_filter_annotation(peptides=peptides)

        # for i, prot in tqdm(enumerate(proteins)):
        #     pep = peptides[i]
        #     for seq in self.filteredAnnotation:
        #         if pep in seq:
        #             self.identifiedUncharacterized[seq] = self.filteredAnnotation[seq]
        #             if '_F:' in prot:
        #                 self.microproteinsSharedPeptides[seq] = self.filteredAnnotation[seq]
        #
        # for seq in self.microproteinsSharedPeptides:
        #     print(self.microproteinsSharedPeptides[seq])
        #
        # for seq in self.identifiedUncharacterized:
        #     print(self.identifiedUncharacterized[seq])

    def intersect_mass_spec(self):
        """
        uniprot annotation with mass spec evidence matching also a novel microprotein
        """
        df = pd.read_csv(f'{self.rescorePostProcessDir}/group/peptides_fixed.txt', sep='\t')
        df = df[df["q-value"] <= 0.01]
        df = df[df["proteinIds"].str.contains("rev_") == False]
        df = df[df["proteinIds"].str.contains("ANNO") == False]
        df = df[df["proteinIds"].str.contains('_F:')]  # gets only rows with a peptide that matched a novel MP
        peptides_unfixed = df["peptide"].tolist()
        proteins = df["proteinIds"].tolist()
        peptides = []
        for pep in peptides_unfixed:
            letters_only = re.findall('[a-zA-Z]', pep)
            peptide = ''.join(letters_only)
            peptides.append(peptide)


        adf = pd.read_csv(self.annotatedIntersectedFile, sep='\t')
        pattern = '|'.join(map(re.escape, peptides))
        adf = adf[adf["Sequence"].str.contains(pattern)]
        anno_df, proof_df = self.count_occurrences(df=adf)

        # annotations that intersected any novel microprotein identified by RP3
        anno_df.to_csv(self.uniprotAnnotationMPIntersection, sep='\t', index=False)
        proof_df.to_csv(self.uniprotEvidenceMPIntersection, sep='\t', index=False)

    def intersect_overall(self):
        """
        uniprot evidence with any mass spec evidence in the data, both regular and microproteins (annotated)
        """

        adf = pd.read_csv(self.annotatedIntersectedFile, sep='\t')
        anno_df, proof_df = self.count_occurrences(df=adf)
        anno_df.to_csv(self.uniprotAnnotationAll, sep='\t', index=False)
        proof_df.to_csv(self.uniprotEvidenceAll, sep='\t', index=False)

        adf = adf[adf["Sequence"].str.len() <= 150]
        anno_df_mp, proof_df_mp = self.count_occurrences(df=adf)
        anno_df_mp.to_csv(self.uniprotAnnotationMP, sep='\t', index=False)
        proof_df_mp.to_csv(self.uniprotEvidenceMP, sep='\t', index=False)


    @staticmethod
    def count_occurrences(df):
        annos = df["Annotation"].tolist()
        proofs = df["Protein existence"].tolist()
        occurrences_anno = {}
        occurrences_proof = {}

        for anno, proof in zip(annos, proofs):
            if anno not in occurrences_anno:
                occurrences_anno[anno] = 0
            if proof not in occurrences_proof:
                occurrences_proof[proof] = 0
            occurrences_anno[anno] += 1
            occurrences_proof[proof] += 1

        data_anno = {'annotation': [], 'count': []}

        for anno in occurrences_anno:
            data_anno['annotation'].append(anno)
            data_anno['count'].append(occurrences_anno[anno])

        data_proof = {'annotation': [], 'count': []}
        for proof in occurrences_proof:
            data_proof['annotation'].append(proof)
            data_proof['count'].append(occurrences_proof[proof])
        anno_df = pd.DataFrame(data=data_anno)
        proof_df = pd.DataFrame(data=data_proof)
        return anno_df, proof_df

    def generate_plots(self):
        files = os.listdir(self.annotatedCountsDir)
        for file in files:
            plt.clf()
            df = pd.read_csv(f'{self.annotatedCountsDir}/{file}', sep='\t')
            pal = sns.color_palette("YlOrBr", as_cmap=True)

            ax = sns.barplot(data=df, y="annotation", x="count", estimator=sum, edgecolor="black", orient='h',
                             color=pal)
            # plt.tight_layout()
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])

            title = file[:-4].replace("_", " ")
            plt.title(title)
            plt.savefig(f'{self.annotatedPlotsDir}/{file[:-4]}.png')
            plt.savefig(f'{self.annotatedPlotsDir}/{file[:-4]}.pdf')

