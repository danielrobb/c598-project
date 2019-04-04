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
                 dt=None,
                 F=None,
                 R=0.35,
                 alpha=0.472,
                 beta=0.0154,
                 gamma=0.055):

        self.F = F
        self.R = R
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.B = 1. #self.gamma / (self.alpha * self.beta)
        self._calc_h1()
        self._calc_h2()
        # Space
        self.nx = nx
        self.xmin = x_span[0]
        self.xmax = x_span[1]
        self.L = self.xmax - self.xmin
        self.dx = self.L / nx
        self.x = np.arange(-self.L/2., self.L/2., self.dx)
        # Wave number
        kl = np.linspace(0, self.nx/2 - 1, self.nx/2)
        kr = np.linspace(1, self.nx/2, self.nx/2)
        kr = -1*kr[::-1]
        self.k = (2.*np.pi/self.L)*np.concatenate((kl, kr))
        # Initial conditions
        self.u = np.zeros_like(self.x)
        # Coefficients
        self._calc_H()
        self._calc_Hx()
        self.alpha0 = self.F - 1.
        self.alpha1 = -3./2. * self.alpha * self._calc_dn(-2)
        self.alpha2 = 3. * self.alpha**2 * self._calc_dn(-3)
        self.beta1 = 1./6. * self.beta * self._calc_dn(-1)
        self.topo = self.beta * self.B  / (2 * self.h2) * self.Hx
        # Time
        self.tmin = t_span[0]
        self.tmax = t_span[1]  
        self.t_out = t_out 
        self.id_out = 0
        self.nout = len(self.t_out)
        self.dt = dt
        self.nt = int(np.ceil(self.tmax/self.dt))
        self.idt = 0
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

    def _kdv(self, u):
        """KdV equation discretized in x."""
        uhat = np.fft.fft(u)
        ux = np.real(np.fft.ifft(1j * self.k * uhat))
        uxxx = np.real(np.fft.ifft(-1j * self.k**3 * uhat))
        diss = np.real(0.005 * np.fft.ifft(np.abs(self.k)**0.5 * (-1 + np.sign(self.k)) * uhat))
        dudt = -1.*(self.alpha0 + self.alpha1*u + self.alpha2*u**2)*ux + self.beta1*uxxx + self.topo + diss
        #print(f'diss: {np.max(diss)} topo: {np.max(self.topo)} disp: {np.max(self.beta1*uxxx)}')
        return dudt

    def _calc_h1(self):
        self.h1 = 1. - self.R

    def _calc_h2(self):
        self.h2 = self.R

    def _calc_H(self):
        self.H = np.ones_like(self.x)
        idx = np.where((self.x >= -0.5) & (self.x <= 0.5))[0]
        self.H[idx] = 1. - 4. * (self.x[idx])**2
        #self.H[idx] = np.cos(np.pi * self.x[idx])**2

    def _calc_Hx(self):
        self.Hx = np.zeros_like(self.x)
        idx = np.where((self.x >= -0.5) & (self.x <= 0.5))[0]
        #self.Hx[idx] = -np.sin(2*np.pi*self.x[idx])
        self.Hx[idx] = -8. * self.x[idx]

    def _calc_dn(self, n):
        return self.h2**n + (-1.)**(n-1)*self.h1**n
    
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

    def to_csv(self, filename='out.csv'):
        """Save model results to a csv."""
        u = self.u_out.flatten('F')
        t = np.tile(self.t_out, self.nx)
        x = np.repeat(self.x, (len(self.t_out)))
        print(u.shape, t.shape, x.shape)
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
    # Grid
    xmin, xmax = (-100., 100.)
    tmin, tmax = (0., 97.4)
    nx = 2048
    dt = 5.e-3
    nout = 6
    t_out = np.linspace(tmin, tmax, nout)

    Fs = [1.2, 1.3, 1.4]
    for F in Fs:
        m = KdVModel(nx=nx, x_span=(xmin, xmax), t_span=(tmin, tmax), t_out=t_out, dt=dt, F=F)
        m.run()
        m.to_csv(f'../data/melville_{F:.1f}.csv')
