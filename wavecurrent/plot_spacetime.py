import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

run = 'wavecurrent_1024'
filename = f'../data/{run}.csv'
df = pd.read_csv(filename)


x = pd.unique(df.x)
t = pd.unique(df.t)
u = df.phi_r.values
print(u.shape, x.shape, t.shape)
nx, nt = (len(x), len(t))
u = np.reshape(u, (nt, nx), order='F')
sns.set_context('paper')
fig, ax = plt.subplots(figsize=(6.5, 5))
X, T = np.meshgrid(x, t)
cb = plt.contourf(X, T, u, levels=np.arange(-0.1, 0.1+1e-6, 0.002), cmap='gist_rainbow', extend='both')
ax.set_xlabel('$x$ (-)')
ax.set_ylabel('$t$ (-)')
ax.set_xlim(0, 1)
ax.set_ylim(0, 4.)
ax.axvline(x=0.75, linewidth=1, linestyle='--', color='red')
fig.colorbar(cb, label='$\phi$ (-)', extendrect=True)
plt.tight_layout()
plt.show()
#plt.savefig(f'../fig/{run}_contour.png', dpi=400)
