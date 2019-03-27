"""Solve the kdv equation for a weakly nonlinear long internal wave."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.fftpack import diff as psdiff

def calc_c0(gp, h1, h2):
    return np.sqrt(gp * h1 * h2 / (h1 + h2))

def calc_alpha(gp, h1, h2):
    return 1.5 * calc_c0(gp, h1, h2) * (h1 - h2) / (h1 * h2)

def calc_beta(gp, h1, h2):
    return calc_c0(gp, h1, h2) * h1 * h2 / 6.

def kdv(t, u, params):
    """Differential equations for the KdV equation, discretized in x."""
    p = params
    k, L, c0, alpha, beta, diss_coeff = (p['k'], p['L'], p['c0'], p['alpha'], p['beta'], p['diss_coeff'])
    uhat = np.fft.fft(u)
    diss = np.real(diss_coeff * np.fft.ifft(np.abs(k)**0.5 * (-1 + np.sign(k)) * uhat))
    ux = psdiff(u, period=L)
    uxxx = psdiff(u, period=L, order=3)
    dudt = -c0*ux -alpha*u*ux - beta*uxxx + diss
    return dudt

def rk4(f, t, y, dt, c):
    k1 = f(t, y, c) 
    k2 = f(t + dt/2., y + dt/2. * k1, c)
    k3 = f(t + dt/2., y + dt/2. * k2, c)
    k4 = f(t + dt, y + dt * k3, c)
    return y + dt/6. * (k1 + 2*k2 + 2*k3 + k4)

if __name__ == "__main__":
    # Set the size of the domain, and create the discretized grid.
    L = 12.
    H = 0.29
    h2 = 0.058
    h1 = H - h2
    nu = 1.e-6
    theta = 0.5 * np.pi / 180.
    g = 9.8
    gp = g * 20 / 1000.
    c0 = calc_c0(gp, h1, h2)
    alpha = calc_alpha(gp, h1, h2)
    beta = calc_beta(gp, h1, h2)
    diss_coeff = 0.5 * np.sqrt(0.5 * nu * c0) * (h1 + h2) / (h1 * h2)

    nx = 256
    x = np.linspace(0.,L, nx+1)
    x = x[:nx]  
    kx1 = np.linspace(0, nx/2-1, nx/2)
    kx2 = np.linspace(1, nx/2, nx/2)
    kx2 = -1*kx2[::-1]
    k = (2.*np.pi/L)*np.concatenate((kx1, kx2))

    u0 = np.empty_like(x)
    idx = np.arange(0, int(nx/2))
    eta0 = np.tan(theta) * (L/4)
    m = eta0 / (L/4)
    x0 = L/4
    u0[idx] = m * (x[idx] - x0) 
    idx = np.arange(int(nx/2), nx) 
    u0[idx] = -m * (x[idx] - 3*x0)

    t0, tf = (0., 300.)
    nout = 300.
    t_eval = np.linspace(t0, tf, nout)
    print("Computing the solution.")
    params = {'k': k, 'L': L, 'c0': c0, 'alpha': alpha, 'beta': beta, 'diss_coeff': diss_coeff}
    dt = 0.01 
    t = np.arange(t0, tf, dt)
    nt = len(t)
    y = np.zeros((nt, nx))
    y[0, :] = 0.5 * u0
    for i in range(nt-1):
        y[i+1, :] = rk4(kdv, t[i], y[i, :], dt, params)
    f = y[:, :int(nx/2)]
    g = np.flip(y[:, int(nx/2):], axis=1)
    fg = f + g
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, fg[:, int(nx/4)])
    plt.show()
    #sol = solve_ivp(kdv, t_span=(t0, tf), y0=0.5*u0, method='RK45', t_eval=t_eval, rtol=1e-10, atol=1e-10)

