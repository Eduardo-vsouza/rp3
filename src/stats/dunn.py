from sklearn.datasets import load_iris
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scikit_posthocs import posthoc_dunn, posthoc_tukey
import scikit_posthocs as sp
import statsmodels.stats.multicomp as multi


from statannotations.Annotator import Annotator

class DunnWithTukey:
    def __init__(self, data_frame):
        self.df = data_frame
        self.df = self.df.dropna()
        # iris_obj = load_iris()
        # iris_df = pd.DataFrame(iris_obj.data, columns=iris_obj.feature_names)
        # iris_df["species"] = [iris_obj.target_names[s] for s in iris_obj.target]
        # data_frame = iris_df
        # self.df = data_frame.select_dtypes(include=[np.number])
        # self.df["species"] = [.target_names[s] for s in self.df.target]
        # print(self.df)
        self.groups = None
        self.result_dunn = None
        self.result_posthoc = None

    def dunn_with_fdr_bh(self):
        # Perform Dunn test with FDR-Benjamini-Hochberg adjustment
        print(self.df["species"])
        self.result_dunn = sp.posthoc_dunn(self.df, p_adjust="fdr_bh", val_col="sepal length (cm)", group_col="species")
        self.groups = self.df.columns

    def posthoc_tukey(self, group, value):
        if self.result_dunn is None:
            raise ValueError("Please perform Dunn test first before running posthoc Tukey.")

        # Prepare data for Tukey's posthoc test
        stacked_data = self.df.stack().reset_index()
        stacked_data.columns = ['index', "group", "value"]

        # Perform Tukey's test
        tukey_result = multi.pairwise_tukeyhsd(stacked_data['value'], stacked_data['group'])
        self.result_posthoc = pd.DataFrame(data=tukey_result._results_table.data[1:],
                                           columns=tukey_result._results_table.data[0])

    def plot_with_pvalues(self):
        if self.result_posthoc is None:
            raise ValueError("Please perform posthoc Tukey test first before plotting.")

        # Plot boxplots
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=self.df, palette="Set3")
        sns.swarmplot(data=self.df, color=".25", alpha=0.7)
        plt.xticks(range(len(self.groups)), self.groups)

        # Add p-values to the plot
        p_values = self.result_posthoc['p-adj']
        annotations = [f"p = {p:.4f}" for p in p_values]
        annotator = Annotator(ax=plt.gca(), y=self.df.values.max(), text_annot_custom=annotations)
        annotator.set_custom_annotations(0.5)

        plt.xlabel('Groups')
        plt.ylabel('Values')
        plt.title('Boxplot with Posthoc Tukey p-values')
        plt.show()


if __name__ == '__main__':
    iris_obj = load_iris()
    iris_df = pd.DataFrame(iris_obj.data, columns=iris_obj.feature_names)
    iris_df["species"] = [iris_obj.target_names[s] for s in iris_obj.target]
    print(iris_df)
    data = DunnWithTukey(data_frame=iris_df)
    data.dunn_with_fdr_bh()
    data.posthoc_tukey(value='sepal length (cm)', group='species')
    data.plot_with_pvalues()

