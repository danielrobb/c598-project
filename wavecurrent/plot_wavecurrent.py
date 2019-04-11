"""Plot linear dispersive wave results."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

Fr_mins = [0.3, 0.4, 0.5, 0.6]
Fr_maxs = [0.7, 0.8, 0.9, 1.0]
for Fr_min, Fr_max in zip(Fr_mins, Fr_maxs):
    fileroot = f'wavecurrent_2048_{Fr_min:.2f}_{Fr_max:.2f}'
    filename = f'../data/{fileroot}.csv'
    df = pd.read_csv(filename)
    ts = pd.unique(df.t)
    x = pd.unique(df.x)
    v = (0.5*(Fr_max - Fr_min) * np.tanh(24*np.cos(2*np.pi*(x-0.25))+23.) - 0.5*(Fr_max + Fr_min))
    sns.set_context('paper', font_scale=1.)
    figsize=(5, 3)
    for i, t in enumerate(ts):
        fig, ax = plt.subplots(figsize=figsize)
        df.loc[df.t == t].plot(ax=ax, x='x', y='phi_r', legend=False)
        y1 = (-v-Fr_min)/25 - 0.045
        ax.plot(x, y1, color='lightgray')
        ax.set_xlabel('$x$')
        ax.set_ylabel('$\phi$')
        ax.text(0.02, 0.98, f'$t =$' +  f' {t:.2f}', va='top', transform=ax.transAxes)
        label = '$Fr_{min} =$' +  f'{Fr_min:.1f}, ' + '$Fr_{max} =$' + f'{Fr_max:.1f}'
        ax.text(0.98, 0.98, label, ha='right', va='top', transform=ax.transAxes, color='lightgray')
        ax.set_ylim(-0.05, 0.05)
        ax.set_xlim(0, 1.)
        filename = f'../fig/{fileroot}_{i:04d}.png'
        plt.tight_layout()
        plt.savefig(filename, dpi=400)
        print(f'saving {filename}')
        plt.close()
