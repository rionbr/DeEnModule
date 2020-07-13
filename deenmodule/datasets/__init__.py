# -*- coding: utf-8 -*-
"""
Biological Boolean Networks
=================================

Datasets used in tutorials


"""
#   Copyright (C) 2020 by
#   Rion Brattig Correia <rionbr@gmail.com>
#   All rights reserved.
#   MIT license.
import os
import pandas as pd
import networkx as nx


_path = os.path.dirname(os.path.realpath(__file__))
""" Make sure we know what the current directory is """


def load_drosophila_pca():
    """PCA files for the Drosophila melanogaster enterocyte gene interaction network

    Returns:
        tuple (pandas.DataFrame)
    """
    df_pca = pd.read_csv(_path + '/drosophila/pca-enterocyte-thr-0p5-DM-dim.csv.gz', index_col=0)
    df_s = pd.read_csv(_path + '/drosophila/pca-enterocyte-thr-0p5-DM-s.csv.gz', index_col=0)
    return df_pca, df_s


def load_sociopatterns_graph():
    """Network file for the SocioPatterns primary school network

    Returns:
        (networkx.Graph)
    """
    G = nx.read_gpickle(_path + 'primary-school-social.gpickle')
    return G
