import os
import pandas as pd
import numpy as np
from koinapy import Koina
from pyteomics import mzml
import matplotlib.pyplot as plt

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
        # concatenate all PIN data
        self._df = pd.concat(
            (pd.read_csv(p, sep='\t', comment='#') for p in pin_files),
            ignore_index=True
        )

    def find_peptide_entries(self, peptide_sequence):
        """
        Return DataFrame rows matching the target sequence exactly.
        """
        pep_cols = [c for c in self._df.columns if 'peptide' in c.lower()]
        if not pep_cols:
            raise ValueError("No peptide column found in PIN files")
        return self._df[self._df[pep_cols[0]] == peptide_sequence]

    def select_best_entry(self, entries):
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
        return entries.iloc[0]

    def load_experimental_spectrum(self, scan_number):
        """
        Pull m/z and intensity arrays for a given scan from the mzML.
        """
        for spec in mzml.read(self.mzml_path):
            if spec.get('id', '').endswith(str(scan_number)) or spec.get('scan number') == scan_number:
                return spec['m/z array'], spec['intensity array']
        raise ValueError(f"Scan {scan_number} not found in {self.mzml_path}")

    def predict_spectrum(self, peptide_sequence, charge):
        """
        Use Koina to predict fragments. Returns arrays of mz, intensity, and annotation.
        """
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
        return mzs, ints, ann

    def plot_butterfly(self, mz_exp, int_exp, mz_pred, int_pred,
                       annotations=None, title=None, annotate=False, out_file=None):
        """
        Mirror plot: experimental (up) vs predicted (down). Optionally annotate ions.
        """
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

    def generate(self, peptide_sequence, out_file=None, annotate=False):
        """
        End-to-end: locate PSM, grab experimental, predict via Koina, and plot.
        Returns the DataFrame of predictions for downstream use.
        """
        entries = self.find_peptide_entries(peptide_sequence)
        if entries.empty:
            raise ValueError(f"Peptide {peptide_sequence} not found in PIN files")
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

# Example usage:
# plotter = ButterflyPlotter(['a.pin'], 'run.mzML', 'Prosit_2020_intensity_HCD', 'koina.wilhelmlab.org:443')\# 
# df_pred = plotter.generate('AAAAAKAK', out_file='butterfly.png', annotate=True)
