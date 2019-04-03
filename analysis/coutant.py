"""
Korteweg–de Vries equation solver surface waves over an obstacle in a 
counterpropagating current.

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
                 ul=None,
                 ur=None,
                 a=None,
                 x0=None,
                 sigma=None,
                 npeaks=None,
                 cfl=1.):
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
        self.x = np.arange(-self.L/2, self.L/2, self.dx)
        # Coefficients
        self.ul = ul
        self.ur = ur
        self.a = a
        self.x0 = x0
        self.sigma = sigma
        self.npeaks = npeaks
        self._calc_v()
        # Initial conditions
        self._calc_u0()
        self.u = self.u0
        # Time
        self.dt = 0.0025
        #self.dt = self.dx**3 * self.cfl
        self.nt = int(np.ceil(self.tmax/self.dt))
        self.idt = 0
        # Wave number
        kl = np.linspace(0, self.nx/2 - 1, self.nx/2)
        kr = np.linspace(1, self.nx/2, self.nx/2)
        kr = -1*kr[::-1]
        self.k = (2.*np.pi/self.L)*np.concatenate((kl, kr))
        # Initial conditions
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
        dudt = -1.*(self.v - 1.)*ux 
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

    def _calc_v(self):
         self.v = (ur + ul)/2. + (ur - ul)/2. * np.tanh(self.a * self.x)
         plt.plot(self.x, self.v)
         plt.show()
     
    def _calc_u0(self):
        x = self.x
        x0 = self.x0
        sigma = self.sigma
        n = self.npeaks
        self.u0 = np.cos(n*np.pi*x/sigma)*np.exp(-np.pi*((x-x0)/sigma)**2)
        plt.plot(self.x, self.u0)
        plt.show()

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
    # Constants
    ul = 1.5
    ur = 0.5
    x0 = 200.
    sigma = 100.
    a = 1./sigma*2
    npeaks = 3.
    # Grid
    xmin, xmax = (-600., 600.)
    tmin, tmax = (0., 1000.)
    nx = 2028
    nout = 10
    t_out = np.linspace(tmin, tmax, nout)
    m = KdVModel(nx=nx, x_span=(xmin, xmax), t_span=(tmin, tmax), t_out=t_out,
                 ul=ul, ur=ur, a=a, x0=x0, sigma=sigma, npeaks=npeaks, cfl=1.)
    m.run()
    m.to_csv('../data/coutant.csv')
