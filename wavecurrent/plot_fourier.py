import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def calc_k(x, L=1):
    nx = len(x)
    kl = np.linspace(0, nx/2-1, nx/2)
    kr = np.linspace(1, nx/2, nx/2)
    kr = -1*kr[::-1]
    k = (2.*np.pi/L)*np.concatenate((kl, kr))
    return k

run = 'wavecurrent_2048_0.4'
filename = f'../data/{run}.csv'
df = pd.read_csv(filename)

x = pd.unique(df.x)
t = pd.unique(df.t)
u = df.phi_r.values + 1j * df.phi_i.values

print(u.shape, x.shape, t.shape)
nx, nt = (len(x), len(t))
u = np.reshape(u, (nt, nx), order='F')
print(u.shape, x.shape, t.shape)
k = calc_k(x)
omega = calc_k(t)
sns.set_context('paper')
fig, ax = plt.subplots(nrows=2, figsize=(8, 5))
#iss = np.arange(0, 498, 70)
iss = [0, 200]
for i in iss:
    ax[0].plot(x, np.real(u[i, :]))
    ax[1].plot(k, np.absolute(np.fft.fft(np.real(u[i, :]))))
ax[0].set_xlabel('$x$ (-)')
ax[1].set_xlabel('$k$ (-)')
ax[0].set_ylabel('$\phi$ (-)')
ax[1].set_ylabel('$\hat{\phi}$ (-)')
plt.tight_layout()
plt.show()

#plt.savefig(f'../fig/{run}_contour.png', dpi=400)


