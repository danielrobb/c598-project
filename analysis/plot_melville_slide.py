"""Plot weakly nonlinear internal waves and compare to horn et al. 2002."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns





Fs = [0.9, 1.1, 1.2]
sns.set_context('talk', font_scale=0.8)
figsize=(6, 4)
fig, ax = plt.subplots(figsize=figsize)
for idx, F in enumerate(Fs):
    filename = f'../data/melville_{F}.csv'
    df = pd.read_csv(filename)
    df.u = df.u * 3 + (idx - 1) * 1.5
    ts = pd.unique(df.t)
    df.loc[df.t == ts[-1]].plot(ax=ax, x='x', y='u', legend=False, color='C0')
    ax.axhline(xmin=-40, xmax=40, y=(idx - 1) * 1.5, color='black', linestyle=':', zorder=0)
    ax.set_xlabel('$x/L$')
    ax.set_ylabel('$\eta/\eta_{max}$')
ax.axvline(ymin=-2, ymax=2, x=0, color='black', linestyle='--', zorder=0)
ax.set_ylim(-2.5, 2.5)
ax.set_xlim(-30, 50)
filename = f'../fig/melville_slide.png'
plt.tight_layout()
plt.savefig(filename, dpi=400)


