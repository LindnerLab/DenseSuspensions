# -*- coding: utf-8 -*-
"""
@author: Lars Kool
"""

def renumber_particles(tracked, rename, start_idx=0):
    """
    Renumber the particle ID's of tracked particles sequentially, starting
    from start_idx.

    Parameters
    ----------
    tracked : n-by-m Pandas Dataframe
        Pandas Dataframe with n elements and m keys to be renumbered. Names of
        keys are arbitrary, except the key that has to be renumbered. That
        key's name should be given by the variable 'rename'.
    rename : str
        Name of the column of which the indices should be renumbered
    start_idx : int, optional
        Start sequential numbering from start_idx. The default is 0.

    Returns
    -------
    tracked : n-by-m Pandas Dataframe
        Pandas Dataframe where the values of key 'rename' are renumbered
        sequentially, starting from start_idx.

    """
    idx = tracked[rename].unique().astype('int')
    nParticles = len(idx)
    i_list = list(range(start_idx, nParticles+start_idx))
    tracked[rename] = tracked[rename].replace(dict(zip(idx, i_list)))
    return tracked