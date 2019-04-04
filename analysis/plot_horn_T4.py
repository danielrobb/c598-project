"""Plot weakly nonlinear internal waves and compare to horn et al. 2002."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

filename = '../data/horn_58_5_T4.csv'
df = pd.read_csv(filename)
df.u = df.u * 100.
ts = pd.unique(df.t)
ts = ts[:8]
print(ts)
sns.set_context('talk', font_scale=0.8)
figsize=(9, 6)
fig, ax = plt.subplots(nrows=4, ncols=2, figsize=figsize, sharex=True, sharey=True)

labels = ['0', '$T_i/4$', '$T_i/2$', '$3/4 \, T_i$', '$T_i$', '$5/4 T_i$', '$3/2 T_i$', '$7/4 T_i$']
for a, t, label in zip(ax.transpose().flat, ts, labels):
    df.loc[df.t == t].plot(ax=a, x='x', y='u', legend=False)
    a.axhline(xmin=0, xmax=6, y=0, color='black', linestyle=':', zorder=0)
    #ax.set_xlabel('$x \; \mathregular{(m)}$')
    #ax.set_ylabel('$\eta \; \mathregular{(cm)}$')
    a.set_ylim(-5, 5)
    a.set_xlim(0, 6)
    a.text(0.02, 0.98, f'$t =$' +  f' {label}', va='top', transform=a.transAxes)

for a in ax[:, 0].flat:
    a.set_ylabel('$\eta \; \mathregular{(cm)}$')
for a in ax[-1, :].flat:
    a.set_xlabel('$x \; \mathregular{(m)}$')
filename = f'../fig/horn_T42.png'
plt.tight_layout(h_pad=0)
plt.savefig(filename, dpi=400)
print(f'saving {filename}')

