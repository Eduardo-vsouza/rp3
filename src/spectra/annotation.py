import os
import sys
import re

import pandas as pd
from pyteomics import pylab_aux as pa, usi
from pyteomics import mzml
from pyteomics import mgf, pepxml, mass, pylab_aux
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus


from Bio import SeqIO
import pylab
import matplotlib.pyplot as plt

from ..pipeline_config import PipelineStructure


class SpectrumAnnotator(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.check_dirs([self.spectrumDir])

        self.microproteins = []
        self.annotatedDir = f'{self.spectrumDir}/annotated_spectra'

    def prepare_input_files(self):
        records = SeqIO.parse(f'{self.summarizedResultsDir}/merged/microproteins_150.fasta_blast_filt.fasta',
                              'fasta')
        for record in records:
            self.microproteins.append(str(record.description))

        data = {'mzxml_file': [], 'scan': [], 'peptide': [], 'protein': [], 'precursor_charge': []}

        generator = self.get_content(self.postProcessDir)
        for content in generator:
            if content.file == 'peptides_filtered.txt':
                df = pd.read_csv(content.fullFile, sep='\t')
                df = df[df["proteinIds"].str.contains('|'.join(self.microproteins))]
                # smorfs = ['STRG.53636.1+chr4:3495610-3495711_F:1_P:5_M', 'STRG.53636.1+chr4:3495610-3495711_F:1_P:5_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.64145.1+chr6:7866472-7866516_F:2_P:9_M', 'MAPS.chr16.2348_-chr16:4416985-4420450-chr16:4420081-4420449_F:0_P:0', 'STRG.30218.1+chr17:56565394-56565474_F:1_P:174_M', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.43875.1+chr3:96414815-96414865_F:2_P:45_M', 'MAPS.chr3.2270_+chr3:96557924-96564079+chr3:96557959-96558087_F:1_P:1', 'STRG.42808.2+chr3:96557968-96558087_F:1_P:0', 'MAPS.chr3.1569_+chr3:96559589-96563900+chr3:96560673-96560720_F:2_P:16', 'STRG.51952.2+chr3:96561321-96561491_F:1_P:42', 'MAPS.chr3.2270_+chr3:96557924-96564079+chr3:96557959-96558087_F:1_P:1', 'STRG.42808.2+chr3:96557968-96558087_F:1_P:0', 'MAPS.chr3.1569_+chr3:96559589-96563900+chr3:96560673-96560720_F:2_P:16', 'STRG.51952.2+chr3:96561321-96561491_F:1_P:42', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.8179.2+chr11:87592769-87592969_F:2_P:12_M', 'STRG.8181.1+chr11:87618948-87618971_F:0_P:19_M', 'STRG.43875.1+chr3:96414815-96414865_F:2_P:45_M', 'MAPS.chr3.2270_+chr3:96557924-96564079+chr3:96557959-96558087_F:1_P:1', 'STRG.42808.2+chr3:96557968-96558087_F:1_P:0', 'MAPS.chr3.1569_+chr3:96559589-96563900+chr3:96560673-96560720_F:2_P:16', 'STRG.51952.2+chr3:96561321-96561491_F:1_P:42', 'MAPS.chr3.2270_+chr3:96557924-96564079+chr3:96557959-96558087_F:1_P:1', 'STRG.42808.2+chr3:96557968-96558087_F:1_P:0', 'MAPS.chr3.1569_+chr3:96559589-96563900+chr3:96560673-96560720_F:2_P:16', 'STRG.51952.2+chr3:96561321-96561491_F:1_P:42', 'STRG.75837.1-chr8:86719860-86719982_F:0_P:193_M', 'MAPS.chr8.905_-chr8:86708331-86746668-chr8:86745857-86746063_F:1_P:6', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'MAPS.chr17.1364_+chr17:36121620-36139566+chr17:36121830-36138765_F:0_P:1_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M', 'STRG.29399.2+chr17:36291728-36291859_F:2_P:61_M']
                # print(df)
                proteins = df["proteinIds"].tolist()
                psmids = df["PSMId"].tolist()
                peptides = df["peptide"].tolist()
                fixed_peptides = []
                for pep in peptides:
                    print(self.reformat_peptide(pep))
                    fixed_peptides.append(self.reformat_peptide(pep).replace("-", ""))

                for i, protein in enumerate(proteins):

                    prot_list = protein.split(",")
                    for prot in prot_list:
                        # if prot in smorfs:
                            data['protein'].append(prot)
                            print(psmids[i])

                            if '/' in psmids[i]:
                                mzml = '_'.join(psmids[i].split("/")[-1].split("_")[:-3])
                                # print(mzml)
                                scan = psmids[i].split("_")[-3]
                                data['precursor_charge'].append(psmids[i].split("_")[-2])

                            else:
                                mzml = psmids[i].split(".")[0]
                                scan = psmids[i].split(".")[1]
                                data['precursor_charge'].append(int(psmids[i].split(".")[-1].split("_")[0]))

                            data['mzxml_file'].append(f'{mzml}.mzML')
                            data['scan'].append(f'{scan}')
                            data['peptide'].append(fixed_peptides[i])

        # search = self.get_content(self.searchDir)
        # for content in search:
        #     if content.file.endswith(".tsv"):
        #
        #         scans = {}
        #         with open(content.fullFile, 'r') as handler:
        #             lines = handler.readlines()
        #             for line in lines:
        #                 cols = line.split("\t")
        #                 scan = cols[0]
        #                 charge = cols[3]
        #                 scans[scan] = charge
        #         # print(scans)
        #         for i, mzml in enumerate(data['mzxml_file']):
        #             if data['db'][i] == content.db:
        #                 if mzml.split(".")[0] == content.file.replace(".tsv", ""):
        #                     print(data['scan'][i], mzml)
        #                     data['precursor_charge'][i] = (scans[data['scan'][i]])
        df = pd.DataFrame(data=data)
        df.to_csv(f'{self.spectrumDir}/spectra_info.txt', sep='\t', index=False)

    @staticmethod
    def reformat_peptide(sequence):
        fixed = []
        mod = False
        for i, aa in enumerate(sequence):
            if aa == '[':
                mod = True
                fixed.append('[+')
            elif aa == ']':
                mod = False
                fixed.append(']')
            if aa == '.':
                if mod:
                    fixed.append(aa)
            else:
                if aa != '.' and aa != '[' and aa != ']' and aa != '-':
                    if i <= len(sequence):
                        # print(len(sequence))
                        # print(i)
                        # print(aa[i+1])
                        if not (i == 0 and sequence[i+1] == '.') and not (i == len(sequence) + 1 and sequence[i + 1] != '.'):
                            fixed.append(aa)
                    else:
                        if not (i == len(sequence) + 1 and sequence[len(sequence)] != '.'):
                            fixed.append(aa)
        fixed_pep = ''.join(fixed)
        print(sequence)
        print(fixed_pep, '\n')
        if not fixed_pep.endswith("]"):
            fixed_pep = fixed_pep[:-1]
        return fixed_pep

    def annotate_spectra(self):
        df = pd.read_csv(f'{self.spectrumDir}/spectra_info.txt', sep='\t')
        cols = {col: df[col].tolist() for col in df.columns}
        self.check_dirs([self.annotatedDir])
        peptides = []
        for i, file in enumerate(cols["mzxml_file"]):
            # print(f'{self.args.mzml}/{file}')
            peptide = cols["peptide"][i]
            if peptide not in peptides:
                # peptides.append(peptide)
                path = self.__search_mzml(file=file.replace(".mzML", "_uncalibrated.mzML"))
                if path != '':
                    self.__annotate_spectrum_with_peptide(mzxml_file=path,
                                                          spectrum_id=str(cols["scan"][i]),
                                                          peptide_sequence=cols["peptide"][i], protein=cols["protein"][i],
                                                          charge=cols["precursor_charge"][i])

    def __search_mzml(self, file):
        groups = os.listdir(self.args.mzml)
        path = ''
        for group in groups:
            files = os.listdir(f'{self.args.mzml}/{group}')
            # print(files)
            # print(file)
            if file in files:
                path = f'{self.args.mzml}/{group}/{file}'
        return path

    def __annotate_spectrum_with_peptide(self, mzxml_file, spectrum_id, peptide_sequence, charge, protein):
        peptide_sequence = peptide_sequence.replace("-", "").replace("n[+42.0106]", "")
        # print(peptide_sequence)
        print(mzxml_file)
        if not os.path.exists(
                f"{self.annotatedDir}/{peptide_sequence}_{protein.replace(':', '')}_{charge}_data.xls"):
            with mzml.read(mzxml_file, use_index=False) as spectra:
                # if spectrum_id != 60937:
                # print(vars(spectra))
                spectrum = spectra.get_by_id(f"controllerType=0 controllerNumber=1 scan={spectrum_id}")
                # spectrum = sus.MsmsSpectrum(identifier=)
                # spectrum = spectrum.annotate_proforma(peptide_sequence, 20, "ppm")
                # print(spectrum)
                plt.clf()

                mz = spectrum['m/z array'].tolist()
                intensities = spectrum['intensity array'].tolist()
                data = {'mz': mz, 'intensities': intensities}
                df = pd.DataFrame(data=data)
                df.to_csv(f"{self.annotatedDir}/{peptide_sequence}_{protein.replace(':', '')}_{charge}_data.xls",
                          sep='\t', index=False)
                fig, ax = plt.subplots(figsize=(12, 6))

                # fig = pylab.figure()
                pylab_aux.annotate_spectrum(spectrum, peptide_sequence,
                                            title=f'{peptide_sequence}',
                                            precursor_charge=charge, backend='spectrum_utils',
                                            ion_types='aby')
                ax.spines["right"].set_visible(False)
                ax.spines["top"].set_visible(False)
                # plt.show()

                fig.savefig(f"{self.annotatedDir}/{peptide_sequence}_{protein.replace(':', '')}_annotated_spectra.png",
                            dpi=600)

# if __name__ == '__main__':
#     data = SpectrumAnnotator(file='/media/eduardo/gold/Eduardo/smORFs_mice_cell_metabolism/cell_ms/20200731_LUM1_CPBA_EASY03_085_30_SA_1886-02_bRP1_07.mzXML')
#     data.annotate_spectrum_with_peptide(spectrum_id='34087', peptide_sequence='KAENGKLVINGEPITIFQERDPTNIKW')
