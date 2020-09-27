"""
Korteweg–de Vries equation solver.

"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class KdVModel(object):
    """Dispersive wave model."""

    def __init__(self,
                 nx=None,
                 x_span=None,
                 t_span=None,
                 t_out=None,
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
        self.x = np.arange(0., self.L, self.dx)
        # Coefficients
        self.c0 = 0.
        self.alpha1 = 1. 
        self.alpha2 = 0. 
        self.beta = 0.022 ** 2 
        #self.dt = self.cfl * (self.dx)**3 / (self.c0) ** 3 
        self.dt = 9.e-6
        self.nt = int(np.ceil(self.tmax/self.dt))
        self.idt = 0
        # Wave number
        kl = np.linspace(0, self.nx/2 - 1, self.nx/2)
        kr = np.linspace(1, self.nx/2, self.nx/2)
        kr = -1*kr[::-1]
        self.k = (2.*np.pi/self.L)*np.concatenate((kl, kr))
        # Initial conditions
        self.u = np.cos(np.pi * self.x)

        dtype = np.float64
        shape = (self.nout, self.nx)
        self.u_out = np.zeros(shape, dtype)

    def run(self):
        """Main loop."""
        t = self.tmin
        while(t < self.tmax):
            if np.isclose(t, self.t_out[self.id_out], atol=self.dt/2.):
                self._print_diagnostics()
                self.u_out[self.id_out, :] = self.u 
                self.id_out += 1
            self._stepforward()
            self.idt += 1
            t += self.dt

    def _print_diagnostics(self):
        idt = self.idt
        id_out = self.id_out
        dt = self.idt
        umin = np.min(self.u)
        umax = np.max(self.u)
        t = self.t_out[self.id_out]
        message = f'iteration: {idt}, output: {id_out}, dt: {dt:.5f},' \
                  + f'time: {t:.2f}, min_u: {umin:.5f}, max_u: {umax:.5f}'
        print(message)    

    def _kdv(self, u):
        """KdV equation discretized in x."""
        uhat = np.fft.fft(u)
        ux = np.real(np.fft.ifft(1j * self.k * uhat))
        uxxx = np.real(np.fft.ifft(-1j * self.k**3 * uhat))
        dudt = -1.*(self.c0 + self.alpha1*u + self.alpha2*u**2)*ux - self.beta*uxxx
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

    def to_csv(self, filename='kdv_out.csv'):
        """Save model results to a csv."""


if __name__ == "__main__":
    xmin, xmax = (0., 2.)
    tmin, tmax = (0., 3.6/np.pi)
    nx = 512
    #t_out = np.array([tmin, np.pi**-1, tmax])
    t_out = [0, 1./np.pi, 3.6/np.pi]
    m = KdVModel(nx=nx, x_span=(xmin, xmax), t_span=(tmin, tmax), t_out=t_out)
    m.run()
    #m.save(filename=filename)
    nt, nx = m.u_out.shape
    sns.set_context('paper', font_scale=0.8)
    width = 3.25
    figsize=(width, width/1.5)
    fig, ax = plt.subplots(figsize=figsize)
    linestyles = [':', '--', '-']
    labels = ['0', '$\pi^{-1}$', '$3.6\pi^{-1}$']
    for i, label, linestyle in zip(range(nt), labels, linestyles):
        ax.plot(m.x, m.u_out[i, :], color='k', linestyle=linestyle, label=f't={label}')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u(x, t)$')
    ax.set_ylim(-1., 3.)
    ax.set_xlim(0, 2)
    ax.legend(loc='upper right')
    filename = f'../fig/zabusky_report.png'
    plt.tight_layout()
    plt.savefig(filename, dpi=800)
    plt.close()


    
