from scipy import stats
# from statannot import add_stat_annotation
from statannotations.Annotator import Annotator

import itertools


class Stater:
    def __init__(self, ax, df, x, xlabel, ylabel, groups_x):
        self.ax = ax
        self.x = x
        self.df = df
        self.xlabel = xlabel
        self.ylabel = ylabel
        # self.x = self.df[self.xlabel].tolist()
        self.normal = self.__check_normality(groups_x)

    def __check_normality(self, groups_x):
        """

        :param groups_x: a dictionary containing group, x values
        :return:
        """
        for group in groups_x:
            stat, p = stats.shapiro(groups_x[group])
            print('Statistics=%.3f, p=%.3f' % (stat, p))
            alpha = 0.05
            if p > alpha:
                normal = True
                print('Sample looks Gaussian (fail to reject H0)')
            else:
                normal = False
                print('Sample does not look Gaussian (reject H0)')
        return normal

    def test(self, order, independent=True, overwrite_pairs=False, loc='inside'):
        """

        :param box_pairs: set
        :return:
        """
        correction = True

        list_combinations = []
        for n in range(len(order) + 1):
            list_combinations += list(itertools.combinations(order, 2))
        # print(list_combinations)
        list_combinations = list(set(list_combinations))
        if overwrite_pairs:
            list_combinations = overwrite_pairs
        if self.normal:
            stat_test = "t-test_ind"
        else:
            if independent:
                if len(order) > 2:
                    correction = True
                    stat_test = 'Kruskal'
                else:
                    correction = False
                    stat_test = 'Mann-Whitney'
        test_results = Annotator(self.ax, pairs=list_combinations, data=self.df, x=self.xlabel, y=self.ylabel, order=order)
        if correction:

            test_results.configure(test=stat_test, text_format='star', loc=loc, comparisons_correction="Bonferroni")
        else:
            test_results.configure(test=stat_test, text_format='star', loc=loc)

        test_results.apply_and_annotate()
        return test_results