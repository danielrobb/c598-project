"""Solve the kdv equation for a weakly nonlinear long internal wave."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.fftpack import diff as psdiff

def c0(gp, h1, h2):
    return np.sqrt(gp * h1 * h2 / (h1 + h2))

def alpha(gp, h1, h2):
    return 1.5 * c0(gp, h1, h2) * (h1 - h2) / (h1 * h2)

def beta(gp, h1, h2):
    return c0(gp, h1, h2) * h1 * h2 / 6.

def kdv(t, u):
    """Differential equations for the KdV equation, discretized in x."""
    # Compute the x derivatives using the pseudo-spectral method.
    L = 12.
    H = 0.29
    h2 = 0.085
    h1 = H - h2
    gp = 10. * 20 / 1000.
    ux = psdiff(u, period=L)
    uxxx = psdiff(u, period=L, order=3)
    dudt = -c0(gp, h1, h2)*ux - alpha(gp, h1, h2)*u*ux - beta(gp, h1, h2)*uxxx
    return dudt

if __name__ == "__main__":
    # Set the size of the domain, and create the discretized grid.
    L = 12.
    N = 256
    H = 0.29
    h2 = 0.085
    h1 = H - h2
    dx = L / (N - 1.0)
    x = np.linspace(0, (1-1.0/N)*L, N)
    u0 = np.empty_like(x)
    idx = np.arange(0, int(N/2))
    theta = 0.75 * np.pi / 180
    eta0 = np.tan(theta) * (L/4)
    m = eta0 / (L/4)
    x0 = L/4
    u0[idx] = m * (x[idx] - x0) 
    idx = np.arange(int(N/2), N) 
    u0[idx] = -m * (x[idx] - 3*x0)
    t0, tf = (0., 300.)
    nout = 1000.
    t_eval = np.linspace(t0, tf, nout)
    print("Computing the solution.")
    sol = solve_ivp(kdv, t_span=(t0, tf), y0=0.5*u0, method='RK45', t_eval=t_eval, rtol=1e-10, atol=1e-10)
    f = sol.y[:int(N/2), :]
    g = np.flip(sol.y[int(N/2):, :], axis=0)
    fg = f + g
    #fig, ax = plt.subplots(nrows=nout, figsize=(10, 6))
    #for idx, (a, t, y) in enumerate(zip(ax, sol.t, sol.y.transpose())):
    #    a.plot(x, y)
    #    a.plot(x[:int(N/2)], fg[:, idx])
    #plt.show()

    fig, ax = plt.subplots(figsize=(10, 6))
    #ax.plot(sol.t, sol.y[int(N/4), :])
    ax.plot(sol.t, fg[int(N/4), :])
    plt.show()

