"""Plot weakly nonlinear internal waves and compare to horn et al. 2002."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



filename = '../data/melville_1.2.csv'
df = pd.read_csv(filename)
df.u = df.u 
ts = pd.unique(df.t)
sns.set_context('talk', font_scale=1.5)
figsize=(9, 6)
for i, t in enumerate(ts):
    fig, ax = plt.subplots(figsize=figsize)
    df.loc[df.t == t].plot(ax=ax, x='x', y='u', legend=False)
    ax.axhline(xmin=-40, xmax=40, y=0, color='black', linestyle=':', zorder=0)
    ax.set_xlabel('$x$')
    #ax.set_ylabel('$\alpha \eta$')
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlim(-40, 40)
    ax.text(0.02, 0.98, f'$t =$' +  f' {t:.0f}', va='top', transform=ax.transAxes)
    filename = f'../fig/melville_1.2_{i:04d}.png'
    plt.tight_layout()
    plt.savefig(filename, dpi=400)
    print(f'saving {filename}')
    plt.close()
