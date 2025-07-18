import os
import types

import pandas as pd
if not hasattr(pd, 'version'):
    pd.version = types.SimpleNamespace(version=pd.__version__)
# import numpy as np
from koinapy import Koina
from pyteomics import mzml
import matplotlib.pyplot as plt
from ..pipeline_config import PipelineStructure
from ..utils import print_state
import re


class ButterflyPlotter:
    """
    Given Percolator .pin files and an mzML file, locate a peptide, retrieve its
    experimental spectrum, predict its spectrum via Koina (Prosit), and produce a
    butterfly (mirror) plot using koinapy-style output (mz, intensities, annotation).
    """
    def __init__(self, pin_files, mzml_path, model_name, server_url, collision_energy=30):
        """
        Parameters:
        - pin_files: list of paths to .pin files
        - mzml_path: path to the mzML file containing spectra
        - model_name: name of the Koina model (e.g. "Prosit_2020_intensity_HCD")
        - server_url: Koina server URL (e.g. "koina.wilhelmlab.org:443")
        - collision_energy: normalized collision energy for Prosit predictions
        """
        self.pin_files = pin_files
        self.mzml_path = mzml_path
        self.collision_energy = collision_energy
        self.model = Koina(model_name, server_url)
        # load and concatenate only relevant columns from PIN files, skipping bad lines
        dfs = []
        for p in pin_files:
            pin = self._read_pin(p)
            if pin is not None:
                dfs.append(self._read_pin(p))
        self._df = pd.concat(dfs, ignore_index=True)

    def _read_pin(self, path):
        print_state(message="reading pin file", level=1, color='yellow')
        """
        Read a .pin file, selecting only essential columns (scan, charge, peptide, scores).
        Skips malformed lines automatically.
        """
        # sniff header to find available columns
        if os.path.getsize(path) == 0:
            # raise ValueError(f"Empty file: {path}")
            return None
        header = pd.read_csv(
            path, sep='\t', comment='#', nrows=0
        ).columns.tolist()
        # columns to extract
        keys = ['scan', 'charge', 'peptide', 'pep', 'score']
        usecols = []
        for key in keys:
            matches = [c for c in header if key in c.lower()]
            if matches:
                usecols.append(matches[0])
        # read only those columns, skip bad lines
        print_state(message=f"Read {path} succesfully.", level=1, color='green')
        return pd.read_csv(
            path,
            sep='\t',
            comment='#',
            usecols=usecols,
            engine='python',
            on_bad_lines='warn'
        )

    def find_peptide_entries(self, peptide_sequence):
        """
        Return DataFrame rows matching the target sequence exactly.
        """
        print_state(message=f"Finding entries for peptide: {peptide_sequence}", level=1, color='yellow')
        pep_cols = [c for c in self._df.columns if 'peptide' in c.lower()]
        if not pep_cols:
            print("No peptide column found in PIN data. Skipping...")
        return self._df[self._df[pep_cols[0]].str.contains(peptide_sequence)]

    def select_best_entry(self, entries):
        print_state(message="Selecting best entry based on PEP or score", level=1, color='yellow')
        """
        Given multiple PSMs, pick one by lowest PEP or highest score.
        """
        score_cols = [c for c in entries.columns if 'pep' in c.lower() or 'score' in c.lower()]
        if score_cols:
            sc = score_cols[0]
            if 'pep' in sc.lower():
                return entries.loc[entries[sc].idxmin()]
            else:
                return entries.loc[entries[sc].idxmax()]
        # print_state
        return entries.iloc[0]

    def load_experimental_spectrum(self, scan_number):
        """
        Pull m/z and intensity arrays for a given scan from the mzML.
        """
        print_state(message=f"Loading experimental spectrum for scan {scan_number}", level=1, color='yellow')
        for spec in mzml.read(self.mzml_path):
            if spec.get('id', '').endswith(str(scan_number)) or spec.get('scan number') == scan_number:
                return spec['m/z array'], spec['intensity array']
        raise ValueError(f"Scan {scan_number} not found in {self.mzml_path}")

    def predict_spectrum(self, peptide_sequence, charge):
        """
        Use Koina to predict fragments. Returns arrays of mz, intensity, and annotation.
        """
        print_state(message=f"Predicting spectrum for {peptide_sequence} with charge {charge}", level=1, color='yellow')
        df_in = pd.DataFrame({
            'peptide_sequences': [peptide_sequence],
            'precursor_charges': [int(charge)],
            'collision_energies': [float(self.collision_energy)]
        })
        pred = self.model.predict(df_in)
        # extract prosit-style output columns
        mzs = pred['mz'].values
        ints = pred['intensities'].values
        ann = pred['annotation'].values
        print_state(message="Prediction completed", level=1, color='green')
        return mzs, ints, ann

    def plot_butterfly(self, mz_exp, int_exp, mz_pred, int_pred,
                       annotations=None, title=None, annotate=False, out_file=None):
        """
        Mirror plot: experimental (up) vs predicted (down). Optionally annotate ions.
        """
        print_state(message="Plotting butterfly plot", level=1, color='yellow')
        plt.figure(figsize=(10, 4))
        plt.bar(mz_exp, int_exp, width=0.5, alpha=0.6, label='Experimental')
        plt.bar(mz_pred, -int_pred, width=0.5, alpha=0.6, label='Predicted')
        if annotate and annotations is not None:
            for x, y, label in zip(mz_pred, -int_pred, annotations):
                plt.text(x, y, label, rotation=90, va='top', fontsize=6)
        plt.xlabel('m/z')
        plt.ylabel('Intensity')
        top = max(int_exp.max() * 1.05, 1.0)
        bot = -int_pred.max() * 1.05
        plt.ylim(bot, top)
        if title:
            plt.title(title)
        plt.legend(loc='upper right')
        plt.tight_layout()
        if out_file:
            plt.savefig(out_file, dpi=300)
        else:
            plt.show()
        print_state(message="Plotting completed", level=1, color='green')

    def generate(self, peptide_sequence, out_file=None, annotate=False):
        """
        End-to-end: locate PSM, grab experimental, predict via Koina, and plot.
        Returns the DataFrame of predictions for downstream use.
        """
        
        if getattr(self, 'args', None) and getattr(self.args, 'proteinSeq', None):
            # Perform a simple tryptic digestion: cut after K or R (unless followed by P)
            peptides = re.split(r'(?<=[KR])(?!P)', self.args.proteinSeq)
            peptides = [pep for pep in peptides if pep]  # remove empty strings
            combined_entries = []
            for pep in peptides:
                pep_entries = self.find_peptide_entries(pep)
                if not pep_entries.empty:
                    combined_entries.append(pep_entries)
            if combined_entries:
                entries = pd.concat(combined_entries, ignore_index=True)
                # Optionally, choose one peptide (e.g., the one with the best score) for downstream processing
                peptide_sequence = entries[[c for c in entries.columns if 'peptide' in c.lower()][0]].iloc[0]
            else:
                return
                raise ValueError("No matching peptides found in PIN data for the given protein sequence")
        else:
            entries = self.find_peptide_entries(peptide_sequence)
        entries = self.find_peptide_entries(peptide_sequence)
        # if entries.empty:
        #     raise ValueError(f"Peptide {peptide_sequence} not found in PIN data")
        best = self.select_best_entry(entries)
        scan_col = [c for c in entries.columns if 'scan' in c.lower()][0]
        charge_col = [c for c in entries.columns if 'charge' in c.lower()][0]
        scan = best[scan_col]
        charge = best[charge_col]

        mz_exp, int_exp = self.load_experimental_spectrum(scan)
        mz_pred, int_pred, ann = self.predict_spectrum(peptide_sequence, charge)

        title = f"{peptide_sequence} (scan {scan}, z={charge})"
        self.plot_butterfly(
            mz_exp, int_exp, mz_pred, int_pred,
            annotations=ann, title=title, annotate=annotate,
            out_file=out_file
        )
        return pd.DataFrame({'mz': mz_pred, 'intensities': int_pred, 'annotation': ann})

class Butterfly(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)
    
    def fly(self):
        self.args.fileFormat = '.mzML'
        pin_files = []
        db = self.select_database(decoy=True).split("/")[-1]
        files = os.listdir(os.path.join(self.searchDir, 'group', db))
        for file in files:
            if file.endswith('.pin'):
                pin_files.append(os.path.join(self.searchDir, 'group', db, file))
        
        self.check_dirs([self.spectrumDir])
        mzml = []
        mzml_files = os.listdir(self.args.mzml)
        for file in mzml_files:
            if file.endswith(self.args.fileFormat):
                mzml_path = os.path.join(self.args.mzml, file)
                mzml.append(mzml_path)
            elif os.path.isdir(os.path.join(self.args.mzml, file)):
                subfiles = os.listdir(os.path.join(self.args.mzml, file))
                for subfile in subfiles:
                    if subfile.endswith(self.args.fileFormat):
                        mzml_path = os.path.join(self.args.mzml, file, subfile)
                        mzml.append(mzml_path)                
        for file in mzml:
            self.print_state(message=f"Processing {file}", color='yellow')
            plotter = ButterflyPlotter(pin_files, file, 'Prosit_2020_intensity_HCD', 'koina.wilhelmlab.org:443')
            plotter.generate(self.args.peptideSeq, out_file=f'{self.spectrumDir}/{self.args.peptideSeq}_{os.path.basename(file)}.png', annotate=True)

# Example usage:
# plotter = ButterflyPlotter(['a.pin'], 'run.mzML', 'Prosit_2020_intensity_HCD', 'koina.wilhelmlab.org:443')\# 
# df_pred = plotter.generate('AAAAAKAK', out_file='butterfly.png', annotate=True)
