"""
Linear dispersive wave equation solver for surface waves over 
topography in the presence of a current.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class WaveModel(object):
    """Linear dispersive wave model."""

    def __init__(self,
                 nx=None,
                 dt=None,
                 tout=None,
                 x_span=None,
                 t_span=None,
                 Fr_span=None,
                 plot_initial_conditions=True):

        self.nx = nx
        self.xmin = x_span[0]
        self.xmax = x_span[1]
        self.tmin = t_span[0]
        self.tmax = t_span[1]
        self.Fr_min = Fr_span[0]
        self.Fr_max = Fr_span[1]
        self.tout = tout
        self.dt = dt
        self.id_out = 0
        self.nout = len(self.tout)
        self.L = self.xmax - self.xmin
        self.dx = self.L / nx
        self.x = np.arange(self.xmin, self.xmax, self.dx)
        self.nt = int(np.ceil(self.tmax/self.dt))
        self.idt = 0
        self._calc_k()
        self._calc_v()
        self._calc_disp2()
        self._init_variables()
        self._init_conditions()
        self._calc_v()
        if plot_initial_conditions is True:
            self._plot_check()

    def run(self):
        """Main loop."""
        self.t = self.tmin
        try:
            while self.t < self.tmax:
                if np.isclose(self.t, self.tout[self.id_out], atol=self.dt/2.):
                    self._print_diagnostics()
                    self.phi_out[self.id_out, :] = self.phi[self.idt, :]
                    self.id_out += 1
                self._stepforward()
                self.idt += 1
                self.t += self.dt
        except:
            print('Simulation finished: t > tmax')

    def _stepforward(self):
        """Step forward."""
        # Velocity potential phi
        self._build_coeff_b()
        self._build_matrix_B()
        self._update_phi()
        # Dispersion term f2
        self._update_f2()
        # Momentum mom
        self._build_coeff_c()
        self._build_matrix_C()
        self._update_mom()

    def _update_phi(self):
        """Update velocity potential phi."""
        m = self.idt
        self.phi[m+1, :] = np.linalg.solve(self.B, self.b[m, :])

    def _update_mom(self):
        """Update momentum mom."""
        m = self.idt
        self.mom[m+1, :] = np.linalg.solve(self.C, self.c[m, :])

    def _update_f2(self):
        """Update dispersion term F^2(-i \partial_x)phi(t,x)."""
        m = self.idt
        phi_hat = np.fft.fft(self.phi[m+1, :])
        f2 = self.disp2
        self.f2[m+1, :] = np.real(np.fft.ifft(f2 * phi_hat))

    def _build_coeff_b(self):
        """Build coefficient b."""
        a, m = (self.a, self.idt)
        for i in range(self.nx):
            (i_prev, i_next) = self._get_prev_next(i)
            self.b[m, i] = (self.phi[m, i]
                            - a*self.v[i]*(self.phi[m, i_next] - self.phi[m, i_prev])
                            + self.dt * self.mom[m, i])

    def _build_coeff_c(self):
        """Build coefficient c."""
        a, m = (self.a, self.idt)
        for i in range(self.nx):
            (i_prev, i_next) = self._get_prev_next(i)
            self.c[m, i] = (self.mom[m, i]
                            - a*(self.v[i_next]*self.mom[m, i_next] - self.v[i_prev]*self.mom[m, i_prev])
                            - self.dt * self.f2[m+1, i])

    def _build_matrix_B(self):
        """Build cyclic tridiagonal matrices B."""
        self.B[self.nx-1, 0] = self.a * self.v[self.nx-1]
        self.B[0, self.nx-1] = self.a * self.v[0] * -1.
        for i in range(self.nx):
            self.B[i, i] = 1.
            if i > 0:
                self.B[i, i-1] = self.a * self.v[i] * -1.
            if i < self.nx-1:
                self.B[i, i+1] = self.a * self.v[i]

    def _build_matrix_C(self):
        """Build cyclic tridiagonal matrices C."""
        self.C[self.nx-1, 0] = self.a * self.v[0]
        self.C[0, self.nx-1] = self.a * self.v[self.nx-1] * -1.
        for i in range(self.nx):
            self.C[i, i] = 1.
            if i > 0:
                self.C[i, i-1] = self.a * self.v[i-1] * -1.
            if i < self.nx-1:
                self.C[i, i+1] = self.a * self.v[i+1]

    def _get_prev_next(self, i):
        """Return left and right spatial points on a periodic domain."""
        i_prev = i - 1
        if i < self.nx-1:
            i_next = i + 1
        else:
            i_next = 0
        return (i_prev, i_next)

    def _init_variables(self):
        """Initialize variables."""
        dtype = np.dtype(np.complex128)
        shape = (self.nt, self.nx)
        # Velocity potential
        self.phi = np.zeros(shape, dtype)
        self.phi0 = np.zeros(self.nx, dtype)
        self.mom = np.zeros(shape, dtype)
        self.mom0 = np.zeros(self.nx, dtype)
        self.B = np.zeros((self.nx, self.nx), dtype)
        # Momentum
        self.mom  = np.zeros(shape, dtype)
        self.C = np.zeros((self.nx, self.nx), dtype)
        # Coefficients
        self.a = self.dt / (4 * self.dx)
        self.b  = np.zeros(shape, dtype)
        self.c  = np.zeros(shape, dtype)
        # Action on phi
        self.f2 = np.zeros(shape, dtype)
        # Background velocity
        self.v = np.zeros(self.nx, dtype)
        # Output
        shape = (self.nout, self.nx)
        self.phi_out = np.zeros(shape, dtype)

    def _init_conditions(self):
        self._calc_phi0()
        self.phi[0, :] = self.phi0
        self._calc_mom0()
        self.mom[0, :] = self.mom0

    def _calc_v(self):
        """Calculate velocity of the background current."""
        x = self.x
        self.v = (0.5*(self.Fr_max - self.Fr_min)
                  * np.tanh(24*np.cos(2*np.pi*(x-0.25))+23.)
                  - 0.5*(self.Fr_max + self.Fr_min))

    def _calc_k(self):
        """Calculate wave numbers."""
        kl = np.linspace(0, self.nx/2-1, self.nx/2)
        kr = np.linspace(1, self.nx/2, self.nx/2)
        kr = -1*kr[::-1]
        self.k = (2.*np.pi/self.L)*np.concatenate((kl, kr))
     
    def _calc_phi0(self):
        """Calculate the initial velocity potential."""
        x = self.x
        x0 = 0.3
        sigma = 0.025
        amp = 0.1
        self.phi0 = amp * np.exp(-0.5*((x-x0)/sigma)**2 + 8j*(x-x0)/sigma)

    def _calc_mom0(self):
        """Calculate the initial momentum."""
        x = self.x
        phi_hat = np.fft.fft(self.phi0)
        f2 = self.disp2
        f = f2 ** 0.5 
        mom0 = -1j * np.fft.ifft(f * phi_hat)
        dispersive_term = np.real(np.fft.ifft(f2 * phi_hat)) 
        for i in range(self.nx):
            (i_prev, i_next) = self._get_prev_next(i)
            self.mom0[i] = (mom0[i]
                            - self.a*(self.v[i_next]*mom0[i_next]
                                      - self.v[i_prev]*mom0[i_prev])
                            - self.dt/2.*dispersive_term[i])

    def _calc_disp2(self):
        k0 = 512.
        self.disp2 = (k0 * np.tanh(self.k / k0))**2

    def to_csv(self, filename='out.csv'):
        """Save model results to a csv."""
        phi = self.phi_out.flatten('F')
        t = np.tile(self.tout, self.nx)
        x = np.repeat(self.x, (len(self.tout)))
        print(phi.shape, t.shape, x.shape)
        data = {'t': t, 'x': x, 'phi_r': np.real(phi), 'phi_i': np.imag(phi)}
        df = pd.DataFrame(data)
        df.to_csv(filename, index=False, float_format='%.5f')

    def _print_diagnostics(self):
        idt = self.idt
        id_out = self.id_out
        dt = self.dt
        umin = np.real(np.min(self.phi))
        umax = np.real(np.max(self.phi))
        mmin = np.real(np.min(self.mom))
        mmax = np.real(np.max(self.mom))
        t = self.t
        message = f'iteration: {idt}, output: {id_out}, dt: {dt:.5f},' \
                  + f'time: {t:.2f}, min_u: {umin:.5f}, max_u: {umax:.5f}, min_m: {mmin:.5f}, max_m: {mmax:.5f}'
        print(message)

    def _plot_check(self):
        plt.plot(self.x, self.phi[0, :])
        plt.plot(self.x, 0.1*self.mom[0, :])
        plt.plot(self.x, self.v)
        plt.show()

if __name__ == "__main__":
    xmin, xmax = (0., 1.)
    tmin, tmax = (0., 4.)
    Fr_min, Fr_max = (0.5, 0.7)
    nx = 1024
    dx = (xmax - xmin)/nx
    dt = dx
    nout = 100
    tout = np.linspace(tmin, tmax, nout)
    m = WaveModel(nx=nx,
                  dt=dt,
                  tout=tout,
                  x_span=(xmin, xmax),
                  t_span=(tmin, tmax),
                  Fr_span=(Fr_min, Fr_max),
                  plot_initial_conditions=True)
    m.run()
    m.to_csv(f'../data/wavecurrent_{nx}.csv')
