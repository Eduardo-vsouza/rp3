import os
import sys
import re
#
from alphabase.peptide.fragment import get_charged_frag_types
from peptdeep.protein.fasta import PredictSpecLibFasta
from peptdeep.pretrained_models import ModelManager
from peptdeep.spec_lib.translate import translate_to_tsv

import pandas as pd

from ..pipeline_config import PipelineStructure


class RTPred(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.rtDir = f"{self.outdir}/RT_prediction"
        self.spectralLibraryHDF = f'{self.rtDir}/microproteins.speclib.hdf'
        self.check_dirs([self.rtDir])

        self.instrument = "Lumos"
        self.__configure_settings()
        self.rtPredictionComparison = f'{self.rtDir}/rt_predictions_comparison.txt'
        self.topKFrag = 16
        self.predictedRTDF = (
            f"{self.spectralLibraryHDF[:-len('.speclib.hdf')]}_frags={self.topKFrag}.speclib.tsv".replace("\'", "")
        )
    def __configure_settings(self):
        print("--Configuring settings to predict retention times with PeptDeep.")
        model_mgr = ModelManager(device='gpu')

        fasta = self.uniqueMicroproteins
        if os.path.exists(self.rescoredMicroproteinsFasta):
            fasta = self.rescoredMicroproteinsFasta

        self.fasta = [fasta]
        # output spectral library in hdf format
        hdf_path = self.spectralLibraryHDF

        protease = self.args.protease
        nce = 30
        instrument = 'Lumos'
        model_mgr.nce = nce
        model_mgr.instrument = instrument
        add_phos = False

        protease_dict = {
            "trypsin": "([KR])",  # this is in fact the "trypsin/P"
            "lysc": "([K])",
            "lysn": "\w(?=K)",
        }
        min_pep_len = self.args.minPepLen
        max_pep_len = self.args.maxPepLen
        max_miss_cleave = self.args.missCleavages
        max_var_mods = 1
        min_pep_mz = 400
        max_pep_mz = 1200
        precursor_charge_min = self.args.precursorCharge[0]
        precursor_charge_max = self.args.precursorCharge[1]

        var_mods = []
        var_mods += ['Oxidation@M']
        # var_mods += ['Phospho@S','Phospho@T','Phospho@Y']

        self.fragTypes = get_charged_frag_types(
            ['b', 'y'] +
            (['b_modloss', 'y_modloss'] if add_phos else []),
            2
        )
        digest = protease_dict[protease]  # Or digest = "trypsin/P"

        self.fastaLib = PredictSpecLibFasta(
            model_mgr,
            protease=digest,
            charged_frag_types=self.fragTypes,
            var_mods=var_mods,
            fix_mods=['Carbamidomethyl@C'],
            max_missed_cleavages=max_miss_cleave,
            max_var_mod_num=max_var_mods,
            peptide_length_max=max_pep_len,
            peptide_length_min=min_pep_len,
            precursor_charge_min=precursor_charge_min,
            precursor_charge_max=precursor_charge_max,
            precursor_mz_min=min_pep_mz,
            precursor_mz_max=max_pep_mz,
            decoy=None
        )

    def prepare_library(self):
        print("--Preparing spectra library.")
        self.fastaLib.get_peptides_from_fasta_list(self.fasta)
        self.fastaLib.append_decoy_sequence()
        self.fastaLib.add_modifications()
        self.fastaLib.precursor_df['nAA'] = self.fastaLib.precursor_df.sequence.str.len()
        self.fastaLib.precursor_df.sort_values('nAA', inplace=True)
        self.fastaLib.precursor_df.reset_index(drop=True, inplace=True)
        self.fastaLib.add_charge()
        self.fastaLib.hash_precursor_df()
        self.fastaLib.calc_precursor_mz()

        print(self.fastaLib.precursor_df)

    def predict_rts(self):
        print("Predicting retention times.")
        self.fastaLib.precursor_df["instrument"] = self.instrument
        res = self.fastaLib.model_manager.predict_all(
            self.fastaLib.precursor_df,
            predict_items=['rt', 'mobility', 'ms2'],
            frag_types=self.fragTypes,
        )
        self.fastaLib.set_precursor_and_fragment(**res)
        self.fastaLib.translate_rt_to_irt_pred()
        print(self.fastaLib.precursor_df.columns)
        print(self.fastaLib.precursor_df["rt_pred"])
        print(self.fastaLib.precursor_df["sequence"])


    def save_library(self):
        print(f"Saving library to {self.spectralLibraryHDF}.")
        self.fastaLib.save_hdf(self.spectralLibraryHDF)

        top_k_frag = 16
        frag_inten = 0.001

        min_frag_mz = 200
        min_frag_nAA = 0


        if 'decoy' in self.fastaLib.precursor_df.columns:
            self.fastaLib._precursor_df = self.fastaLib.precursor_df[self.fastaLib._precursor_df.decoy == 0]
        translate_to_tsv(
            self.fastaLib,
            self.predictedRTDF,
            keep_k_highest_fragments=top_k_frag,
            min_frag_nAA=min_frag_nAA,
            min_frag_mz=min_frag_mz,
            min_frag_intensity=frag_inten,
            batch_size=100000,
            translate_mod_dict=None,
            multiprocessing=True,
        )
        self.fastaLib.load_hdf(self.spectralLibraryHDF, load_mod_seq=True)
        print(f"Spectral library saved to {self.predictedRTDF}")

    def __get_predicted(self):
        predicted_rts = pd.read_csv(self.predictedRTDF, sep='\t')
        predicted = {}
        charges = predicted_rts["PrecursorCharge"].tolist()
        peptides = predicted_rts["ModifiedPeptide"].tolist()
        rts = predicted_rts["RT"].tolist()
        for i, pep in enumerate(peptides):
            pep = f'{pep.replace("_", "")}_{charges[i]}'
            if pep not in predicted:
                predicted[pep] = {'charge': charges[i], 'RT': rts[i]}
        return predicted

    def __fix_pin(self):
        with open(f'{self.outdir}/all_pins.pin', 'r') as handler, open(f'{self.rtDir}/all_pins_fixed.pin', 'w') as outfile:
            lines = handler.readlines()
            new_lines = []
            for line in lines:
                cols = '\t'.join(line.split('\t')[:20]).rstrip()
                new_lines.append(f'{cols}\n')
            outfile.writelines(new_lines)
    def compare_rts(self):
        self.__fix_pin()
        predicted = self.__get_predicted()

        peptides = pd.read_csv(f'{self.rtDir}/all_pins_fixed.pin', sep='\t')
        psms = peptides["SpecId"].tolist()
        # print(psms)
        mods = {'57.0215': 'Carbamidomethyl', "15.9949": "Oxidation"}
        peps = peptides["Peptide"].tolist()
        exp_rts = peptides["retentiontime"].tolist()

        data = {"peptide": [], "charge": [], "predicted_rt": [], "experimental_rt": []}

        for i, pep in enumerate(peps):
            # print(psms[i])
            charge = psms[i].split(".")[-1].split("_")[0]
            seq = pep
            # for i, string in enumerate(strings):
            for mod, replacement in mods.items():
                pattern = r'\[' + re.escape(mod) + r'\]'
                seq = re.sub(pattern, '[' + replacement + ']', seq)
            seq = seq.split('.')[1]
            seq = f'{seq.replace(".", "").replace("n", "").replace("-", "")}_{charge}'

            if seq in predicted:
                data['peptide'].append(seq)
                data['charge'].append(charge)
                data['predicted_rt'].append(predicted[seq]['RT'])
                data['experimental_rt'].append(exp_rts[i])
        df = pd.DataFrame(data=data)
        df.to_csv(self.rtPredictionComparison, sep='\t', index=False)

