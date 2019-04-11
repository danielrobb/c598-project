import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

cmaps = ['jet', 'viridis', 'plasma', 'magma', 'inferno', 'RdBu', 'gray', 'RdYlBu', 'Blues', 'bone']

dF = 0.4
F_maxs = np.array([0.5, 0.6, 0.7, 1.0, 1.1, 1.2])
F_mins = F_maxs - dF
sns.set_context('paper', font_scale=0.8)

for cmap in cmaps:
    fig, ax = plt.subplots(nrows=2, ncols=3, sharex='row', sharey='row', figsize=(6.5, 4))
    for a, F_min, F_max in zip(ax.flat, F_mins, F_maxs):
        run = f'wavecurrent_2048_{F_min:0.2f}_{F_max:0.2f}'
        filename = f'../data/{run}.csv'
        df = pd.read_csv(filename)
        x = pd.unique(df.x)
        t = pd.unique(df.t)
        u = df.phi_r.values
        nx, nt = (len(x), len(t))
        u = np.reshape(u, (nt, nx), order='F')
        threshold = 0.0005
        ind = np.where(np.abs(u) < threshold)
        u[ind] = 0.
        X, T = np.meshgrid(x, t)
        a.contourf(X, T, u, levels=np.arange(-0.025, 0.025, 0.001), cmap=cmap, extend='both')
        a.set_xlim(0.5, 1.)
        if F_max < 1.:
            a.set_ylim(0, 1.4)
        else:
            a.set_ylim(0, 4.)
        fig_title = '$Fr_{min} =$' +  f'{F_min:.1f}, ' + '$Fr_{max} =$' + f'{F_max:.1f}'
        a.set_title(fig_title, fontsize=8)
        a.axvline(x=0.65, linewidth=0.5, linestyle='--', color='black')
        a.axvline(x=0.85, linewidth=0.5, linestyle='--', color='black')
    for a in ax[:, 0].flat:
        a.set_ylabel('$t$')
    for a in ax[1, :].flat:
        a.set_xlabel('$x$')
    plt.tight_layout()
    #plt.show()
    plt.savefig(f'../fig/spacetime_{cmap}.png', dpi=400)
