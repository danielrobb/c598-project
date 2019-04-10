"""Plot linear dispersive wave results."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

fileroot = 'wavecurrent_4096'
filename = f'../data/{fileroot}.csv'
df = pd.read_csv(filename)
ts = pd.unique(df.t)
sns.set_context('paper', font_scale=1.)
figsize=(5, 3)
for i, t in enumerate(ts):
    fig, ax = plt.subplots(figsize=figsize)
    df.loc[df.t == t].plot(ax=ax, x='x', y='phi_r', legend=False)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$\phi$')
    ax.text(0.02, 0.98, f'$t =$' +  f' {t:.2f}', va='top', transform=ax.transAxes)
    ax.set_ylim(-0.11, 0.11)
    ax.set_xlim(0, 1.)
    filename = f'../fig/{fileroot}_{i:04d}.png'
    plt.tight_layout()
    plt.savefig(filename, dpi=400)
    print(f'saving {filename}')
    plt.close()
