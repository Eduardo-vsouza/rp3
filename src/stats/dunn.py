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
        print(self.df)
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

    def dunn_with_fdr_bh(self, val_col, group_col, test='tukey'):
        values = self.df[val_col].tolist()
        cols = self.df[group_col].tolist()
        data = {}
        for val, group in zip(values, cols):
            if group not in data:
                data[group] = []
            data[group].append(val)
        unpack = [data[group] for group in data]

        a = stats.kruskal(*unpack)
        print(a)
        # Perform Dunn test with FDR-Benjamini-Hochberg adjustment
        if test == 'dunn':
            self.result_posthoc = sp.posthoc_dunn(self.df, p_adjust="fdr_bh", val_col=val_col, group_col=group_col)
        elif test == 'tukey':
            print(stats.f_oneway(*unpack))
            self.result_posthoc = sp.posthoc_tukey(self.df, val_col=val_col, group_col=group_col)
        self.groups = self.df.columns
        self.valCol = val_col
        self.groupCol = group_col

    def posthoc_tukey(self):
        if self.result_dunn is None:
            raise ValueError("Please perform Dunn test first before running posthoc Tukey.")

        # Prepare data for Tukey's posthoc test
        # stacked_data = self.df.stack().reset_index()
        # stacked_data.columns = ['index', "group", "value"]
        # print(stacked_data)

        # Perform Tukey's test
        tukey_result = posthoc_tukey(self.df, val_col=self.valCol, group_col=self.groupCol)
        # self.result_posthoc = pd.DataFrame(data=tukey_result._results_table.data[1:],
        #                                    columns=tukey_result._results_table.data[0])
        self.result_posthoc = tukey_result
        # print(self.result_posthoc)


    def plot_with_pvalues(self, order, ylabel, xlabel, loc='outside'):
        # if self.result_dunn is None:
        #     raise ValueError("Please perform posthoc Tukey test first before plotting.")

        remove = np.tril(np.ones(self.result_posthoc.shape), k=0).astype("bool")
        self.result_posthoc[remove] = np.nan

        molten_df = self.result_posthoc.melt(ignore_index=False).reset_index().dropna()
        print(molten_df)
        molten_df = molten_df[molten_df["value"] <= 0.05]
        # print(molten_df)
        species = np.unique(self.df[self.groupCol])
        # print
        sns.set_palette(palette='coolwarm_r')
        # colors = ["#9E2A2B", "#3B604D", "#6E918C", "#639B6D", "#1D8BAC"]
         # = sns.color_palette("viridis", 6)
        # ax = sns.catplot(data=self.df, x=self.groupCol, y=self.valCol, kind='box', order=order)
        # plt.ylim(0,100)
        ax = sns.boxplot(data=self.df, x=self.groupCol, y=self.valCol, order=order, showfliers=False, palette=colors)
        # ax = sns.violinplot(data=self.df, x=self.groupCol, y=self.valCol, order=order, showfliers=False)
        # plt.show()
        pairs = [(i[1]["index"], i[1]["variable"]) for i in molten_df.iterrows()]
        p_values = [i[1]["value"] for i in molten_df.iterrows()]
        print(pairs)
        if pairs:
            annotator = Annotator(
                ax, pairs, data=self.df, x=self.groupCol, y=self.valCol, order=order
            )
            # print(molten_df)
            annotator.configure(text_format="star", loc=loc)
            annotator.set_pvalues_and_annotate(p_values)

        # molten_df
        #
        # # Plot boxplots
        # plt.figure(figsize=(10, 6))
        # sns.boxplot(data=self.df, palette="Set3")
        # sns.swarmplot(data=self.df, color=".25", alpha=0.7)
        # plt.xticks(range(len(self.groups)), self.groups)
        #
        # # Add p-values to the plot
        # p_values = self.result_posthoc['p-adj']
        # annotations = [f"p = {p:.4f}" for p in p_values]
        # annotator = Annotator(ax=plt.gca(), y=self.df.values.max(), text_annot_custom=annotations)
        # annotator.set_custom_annotations(0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.xticks(rotation=45, ha="center")

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        # plt.title('Boxplot with Posthoc Tukey p-values')
        plt.show()


if __name__ == '__main__':

    iris_obj = load_iris()
    iris_df = pd.DataFrame(iris_obj.data, columns=iris_obj.feature_names)
    iris_df["species"] = [iris_obj.target_names[s] for s in iris_obj.target]
    # print(iris_df)
    data = DunnWithTukey(data_frame=iris_df)
    data.dunn_with_fdr_bh(val_col='sepal length (cm)', group_col='species')
    # data.posthoc_tukey()
    data.plot_with_pvalues()

