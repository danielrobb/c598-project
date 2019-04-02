"""Plot weakly nonlinear internal waves at two gauges."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

filename = '../data/horn_58_5.csv'
df = pd.read_csv(filename)
df.u = df.u * 100.
ts = pd.unique(df.t)
xs = pd.unique(df.x)
sns.set_context('talk', font_scale=1.25)
gauges = np.array([3., 4.5])
ids = list()
for i, gauge in enumerate(gauges):
    ids.append(np.argmin(np.abs(gauge - xs)))
print(ids)
x_gauges = xs[ids]
figsize=(9, 7.5)
fig, ax = plt.subplots(nrows=2, figsize=figsize, sharex=True, sharey=True)
labels = ['a. Wave gauge at $x = 3.0 \, \mathregular{m}$', 'b. Wave gauge at $x = 4.5 \, \mathregular{m}$']
for a, x_gauge, label in zip(ax, x_gauges, labels):
    print(x_gauge)
    df.loc[df.x == x_gauge].plot(ax=a, x='t', y='u', legend=False)
    a.set_ylabel('$\eta \; \mathregular{(cm)}$')
    a.set_ylim(-2, 2)
    a.set_yticks(np.arange(-2, 2.1, 1))
    a.set_xlim(0, 300)
    a.text(0., 1.01, label, va='bottom', transform=a.transAxes)
filename = f'../fig/horn_gauge.png'
ax[-1].set_xlabel('$t \; \mathregular{(s)}$')
plt.tight_layout()
plt.savefig(filename, dpi=400)
print(f'saving {filename}')
