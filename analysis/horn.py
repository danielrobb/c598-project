"""
Korteweg–de Vries equation solver for internal waves.

"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns

class KdVModel(object):
    """Dispersive wave model."""

    def __init__(self,
                 nx=None,
                 x_span=None,
                 t_span=None,
                 t_out=None,
                 gp=None,
                 h1=None,
                 h2=None,
                 theta=None,
                 cfl=1.,):
        self.nx = nx
        self.xmin = x_span[0]
        self.xmax = x_span[1]
        self.tmin = t_span[0]
        self.tmax = t_span[1]
        self.t_out = t_out
        self.id_out = 0
        self.nout = len(self.t_out)
        self.cfl = cfl
        self.L = self.xmax - self.xmin
        self.dx = self.L / nx
        self.x = np.arange(0., self.L, self.dx)
        # Coefficients
        self.gp = gp
        self.h1 = h1
        self.h2 = h2
        self.theta = theta
        self._calc_c0()
        self._calc_alpha1()
        self._calc_alpha2()
        self._calc_beta()
        self._calc_diss_coeff()
        # Time
        self.dt = 0.001
        #self.dt = self.dx**3 / self.c0**3 * self.cfl
        self.nt = int(np.ceil(self.tmax/self.dt))
        self.idt = 0
        # Wave number
        kl = np.linspace(0, self.nx/2 - 1, self.nx/2)
        kr = np.linspace(1, self.nx/2, self.nx/2)
        kr = -1*kr[::-1]
        self.k = (2.*np.pi/self.L)*np.concatenate((kl, kr))
        # Initial conditions
        u0 = np.empty_like(self.x)
        idx = np.arange(0, int(self.nx/2))
        eta0 = np.tan(theta) * (self.L/4)
        m = eta0 / (self.L/4)
        x0 = self.L/4
        u0[idx] = m * (self.x[idx] - x0) 
        idx = np.arange(int(self.nx/2), self.nx) 
        u0[idx] = -m * (self.x[idx] - 3*x0)
        self.u = 0.5 * u0
        # Initialize u_out array
        dtype = np.float64
        shape = (self.nout, self.nx)
        self.u_out = np.zeros(shape, dtype)

    def run(self):
        """Main loop."""
        t = self.tmin
        try:
            while t < self.tmax:
                if np.isclose(t, self.t_out[self.id_out], atol=self.dt/2.):
                    self._print_diagnostics()
                    self.u_out[self.id_out, :] = self.u 
                    self.id_out += 1
                self._stepforward()
                self.idt += 1
                t += self.dt
        except:
            print('foo')
        self.wrap_solution()

    def wrap_solution(self):
        self.sol = self.u_out[:, :int(nx/2)] + np.flip(self.u_out[:, int(nx/2):], axis=1)

    def _kdv(self, u):
        """KdV equation discretized in x."""
        uhat = np.fft.fft(u)
        ux = np.real(np.fft.ifft(1j * self.k * uhat))
        uxxx = np.real(np.fft.ifft(-1j * self.k**3 * uhat))
        diss = np.real(self.diss_coeff * np.fft.ifft(np.abs(self.k)**0.5 * (-1 + np.sign(self.k)) * uhat))
        dudt = -1.*(self.c0 + self.alpha1*u + self.alpha2*u**2)*ux - self.beta*uxxx + diss
        return dudt
    
    def _rk4(self, u):
        """Fourth Order Runge-Kutta."""
        dt = self.dt
        k1 = self._kdv(u) 
        k2 = self._kdv(u + dt/2.*k1)
        k3 = self._kdv(u + dt/2.*k2)
        k4 = self._kdv(u + dt*k3)
        return u + dt/6.*(k1 + 2*k2 + 2*k3 + k4)

    def _stepforward(self):
        """Step forward."""
        self.u = self._rk4(self.u)

    def _calc_c0(self):
         gp, h1, h2 = (self.gp, self.h1, self.h2)
         self.c0 = np.sqrt(gp * h1 * h2 / (h1 + h2))
     
    def _calc_alpha1(self):
        c0, h1, h2 = (self.c0, self.h1, self.h2)
        self.alpha1 = 1.5 * c0 * (h1 - h2) / (h1 * h2)

    def _calc_alpha2(self):
        self.alpha2 = 0.
    
    def _calc_beta(self):
         c0, h1, h2 = (self.c0, self.h1, self.h2)
         self.beta = c0 * h1 * h2 / 6.

    def _calc_diss_coeff(self):
        nu = 1.e-6
        c0, h1, h2 = (self.c0, self.h1, self.h2)
        self.diss_coeff = 0.5 * np.sqrt(0.5 * nu * c0) * (h1 + h2) / (h1 * h2)

    def to_csv(self, filename='out.csv'):
        """Save model results to a csv."""
        u = self.sol.flatten('F')
        t = np.tile(self.t_out, (int(self.nx/2)))
        x = np.repeat(self.x[:int(self.nx/2)], (len(self.t_out)))
        data = {'t': t, 'x': x, 'u': u}
        df = pd.DataFrame(data)
        df.to_csv(filename, index=False, float_format='%.5f')

    def _print_diagnostics(self):
        idt = self.idt
        id_out = self.id_out
        dt = self.dt
        umin = np.min(self.u)
        umax = np.max(self.u)
        t = self.t_out[self.id_out]
        message = f'iteration: {idt}, output: {id_out}, dt: {dt:.5f},' \
                  + f'time: {t:.2f}, min_u: {umin:.5f}, max_u: {umax:.5f}'
        print(message)    

if __name__ == "__main__":
    # Constants
    H = 0.29
    h2 = 0.058
    h1 = H - h2
    theta = 0.5 * np.pi / 180.
    g = 9.8
    gp = g * 20. / 1000.
    c = np.sqrt(gp * h1 * h2 / (h1 + h2))
    L = 6
    T = 2 * L / c
    print(f'T/4 = {T/4} seconds')
    # Grid
    xmin, xmax = (0., 12.)
    tmin, tmax = (0., 300.)
    nx = 512
    #nout = 600
    #t_out = np.linspace(tmin, tmax, nout)
    t_out = np.arange(tmin, tmax + T/4, T/4)
    m = KdVModel(nx=nx, x_span=(xmin, xmax), t_span=(tmin, tmax), t_out=t_out,
                 gp=gp, h1=h1, h2=h2, theta=theta, cfl=0.2)
    m.run()
    m.to_csv('../data/horn_58_5_T4.csv')
