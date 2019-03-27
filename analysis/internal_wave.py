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
    """The KdV equation discretized in x."""
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

def solve_kdv(y0, dt, t_span, t_eval, params):
    t0, tf = t_span
    t, y = (t0, y0)
    nt, nx = (len(t_eval), len(y0))
    sol = np.zeros((nt, nx))
    id_out = idt = 0
    while t < tf:
        y = rk4(kdv, t, y, dt, params)
        if np.isclose(t, t_eval[id_out], atol=dt/2):
            print_diagnostics(idt, id_out, dt, t, y)
            sol[id_out, :] = y
            id_out += 1
        t += dt
        idt += 1
    sol_wrapped = wrap_solution(sol)
    return sol_wrapped

def wrap_solution(sol):
    return sol[:, :int(nx/2)] + np.flip(sol[:, int(nx/2):], axis=1)

def print_diagnostics(idt, id_out, dt, t, u):
    message = f'iteration: {idt}, output: {id_out}, dt: {dt:.5f}, time: {t:.2f}, min_u: {np.min(u):.5f}, max_u: {np.max(u):.5f}'
    print(message)    

if __name__ == "__main__":
    # Constants
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
    print(f'c0: {c0}, alpha: {alpha}, beta: {beta}, diss_coeff: {diss_coeff}')
    # Grid
    nx = 512
    dx = L / nx
    x = np.linspace(0.,L, nx+1)
    x = x[:nx]  
    kx1 = np.linspace(0, nx/2-1, nx/2)
    kx2 = np.linspace(1, nx/2, nx/2)
    kx2 = -1*kx2[::-1]
    k = (2.*np.pi/L)*np.concatenate((kx1, kx2))
    # Initial conditions
    u0 = np.empty_like(x)
    idx = np.arange(0, int(nx/2))
    eta0 = np.tan(theta) * (L/4)
    m = eta0 / (L/4)
    x0 = L/4
    u0[idx] = m * (x[idx] - x0) 
    idx = np.arange(int(nx/2), nx) 
    u0[idx] = -m * (x[idx] - 3*x0)
    params = {'k': k, 'L': L, 'c0': c0, 'alpha': alpha, 'beta': beta, 'diss_coeff': diss_coeff}
    # Time stepping
    t_span = (0., 300.)
    (t0, tf) = t_span
    nout = 300.
    t_eval = np.linspace(t0, tf, nout+1)    
    CFL = 0.3
    dt = (dx)**3 / (c0) ** 3 * CFL
    y0 = 0.5 * u0
    sol = solve_kdv(y0, dt, t_span, t_eval, params)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_eval, sol[:, int(nx/4)])
    plt.show()

