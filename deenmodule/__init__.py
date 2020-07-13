# -*- coding: utf-8 -*-
"""
DeEnModule
================



"""
#   Copyright (C) 2020 by
#   Rion Brattig Correia <rionbr@gmail.com>
#   All rights reserved.
#   MIT license.
__package__ = 'deenmodule'
__title__ = u'DeEnModule: decomposition and entropy modules in complex networks'
__description__ = u'This package implements methods to identify modules in decompositions of adjacency matrixes of complex networks.'

__author__ = """\n""".join([
    'Rion Brattig Correia <rionbr@gmail.com>',
    'Marijn ten Thij <mctenthij@gmail.com >'
])

__copyright__ = u'2020, R.B. Correia, M. ten Thij'

__version__ = '0.1'
__release__ = '0.1-alpha'
#
__all__ = ['DeEnModule', 'datasets']


import numpy as np
import pandas as pd
from scipy import stats


def _in_ranges(x, bins):
    """Function for pandas.apply() that assigs values into bins
    """
    return [((x >= lower) & (x <= upper)) for lower, upper in bins]


def DeEnModule(df_dec,
               radius_window=1.0,
               radius_overlap=0.1,
               angle_window=30,
               angle_overlap=15,
               min_points=10,
               min_radius=1.0,
               max_cut_points=3,
               min_radius_separation=1.0,
               components=9,
               verbose=False):
    """ Computes DeEnModule

    Params:

    Returns:
        df_ent )()

    """
    #
    df_dec = df_dec.copy()
    #
    angle_items = int(angle_window / angle_overlap)
    a = np.arange(-180, (181), angle_overlap)
    angle_bins = [(i, j) for i, j in zip(a[0:-angle_items], a[angle_items:])]
    n_bins = len(angle_bins)
    max_entropy = stats.entropy((np.ones(shape=n_bins) / n_bins), base=2)
    list_df_ent = list()
    #
    for dim in range(1, (components + 1)):
        #
        cx = str(dim) + 'c'
        cy = str(dim + 1) + 'c'
        dist_label = '{cx:s}-{cy:s}-dist'.format(cx=cx, cy=cy)
        angle_label = '{cx:s}-{cy:s}-angle'.format(cx=cx, cy=cy)
        df_dec[dist_label] = np.hypot(df_dec[cx], df_dec[cy])
        df_dec[angle_label] = np.degrees(np.arctan2(df_dec[cy], df_dec[cx]))
        #
        df_dec.sort_values(dist_label, ascending=True, inplace=True)
        radius_max = df_dec[dist_label].max()
        #
        radius_items = int(radius_window / radius_overlap)
        #
        b = np.arange(0, (radius_max + radius_overlap), radius_overlap)
        radius_intervals = [(s, e) for s, e in zip(b[0:-radius_items], b[radius_items:])]

        # Loop radius intervals
        r = []
        for radius_start, radius_end in radius_intervals:

            df_dian_tmp = df_dec.loc[(df_dec[dist_label] >= radius_start) & (df_dec[dist_label] <= radius_end), :]

            dfc = df_dian_tmp[angle_label].apply(lambda x: pd.Series(_in_ranges(x, angle_bins), angle_bins))

            if len(dfc) > min_points:
                dfp = (dfc.sum(axis=0) / dfc.sum().sum()).rename('prob').to_frame()
                dfp['log2'] = dfp['prob'].apply(np.log2)
                #
                entropy = stats.entropy(dfp['prob'], base=2)

            else:
                entropy = np.nan

            entropy_norm = entropy / max_entropy
            r.append((dim, radius_start, radius_end, entropy, entropy_norm))

        #
        df_ent_tmp = pd.DataFrame(r, columns=['dim', 'radius-start', 'radius-end', 'entropy', 'entropy-norm'])
        #
        # Interpolation
        df_ent_tmp['entropy-norm-smooth'] = df_ent_tmp['entropy-norm'].interpolate(method='linear', limit_direction='both')
        # Rank
        df_ent_tmp['radius-rank'] = df_ent_tmp['radius-start'].rank(method='min')
        df_ent_tmp['entropy-rank'] = df_ent_tmp['entropy-norm'].rank(method='min')
        # Rank Sum
        df_ent_tmp['rank-sum'] = ((df_ent_tmp['radius-rank']) + (df_ent_tmp['entropy-rank']))
        #
        # Define cut points
        cut_points = []
        # Index % Sort
        df_cp = df_ent_tmp.sort_values('rank-sum').loc[(df_ent_tmp['radius-start'] > min_radius), :]
        possible_rank = 1
        for possible_id, row in df_cp.iterrows():
            possible_value = row['radius-start']
            if not any([True if abs(possible_value - existing_value) <= min_radius_separation else False for existing_id, existing_value, existing_rank in cut_points]):
                cut_points.append((possible_id, possible_value, possible_rank))
                possible_rank += 1

            if len(cut_points) >= max_cut_points:
                break
        #
        dict_cut_points = {idx: rank for idx, value, rank in cut_points}
        df_ent_tmp['circle-rank'] = df_ent_tmp.index.map(dict_cut_points)

        # Add to list
        list_df_ent.append(df_ent_tmp)

    df_en = pd.concat(list_df_ent, axis='index')
    #
    return df_en, df_dec
