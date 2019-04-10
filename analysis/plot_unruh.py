"""Plot weakly nonlinear internal waves and compare to horn et al. 2002."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

filename = '../data/unruh_2048_1.0.csv'
df = pd.read_csv(filename)
df.phi = df.phi_r
ts = pd.unique(df.t)
sns.set_context('talk', font_scale=1.5)
figsize=(9, 6)
for i, t in enumerate(ts):
    fig, ax = plt.subplots(figsize=figsize)
    df.loc[df.t == t].plot(ax=ax, x='x', y='phi_r', legend=False)
    ax.axhline(xmin=0, xmax=6, y=0, color='black', linestyle=':', zorder=0)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$\phi$')
    ax.set_ylim(-0.15, 0.15)
    ax.set_xlim(df.x.min(), df.x.max())
    ax.text(0.02, 0.98, f'$t =$' +  f' {t:.2f}', va='top', transform=ax.transAxes)
    filename = f'../fig/unruh2048_1.0_{i:04d}.png'
    plt.tight_layout()
    plt.savefig(filename, dpi=400)
    print(f'saving {filename}')
    plt.close()
