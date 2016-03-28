"""This module implements the shallow water equations with rotation.

References:
    [1] Griffith, A. and Nichols, N. (2000). Adjoint methods in data
        assimilation for estimating model error. Flow, Turbulence, and
        Combustion, 65: 469--488.

    [2] Pedlosky, J. (1987). Geophysical Fluid Dynamics. Second
        edition.

    [3] Parrett, C. A., & Cullen, M. J. P. (1984). Simulation of
        hydraulic jumps in the presence of rotation and mountains.
        Quarterly Journal of the Royal Meteorological Society,
        147--165.

Author:             Jeremy Shaw
Institution:        Portland State University
Date Created:       23 March 2016
Last Modified Date: 28 March 2016
"""
import numpy as np
import math
import matplotlib.pyplot as plt

import swe_fortran as fortran

class swe(object):
    """Two-dimensional shallow water equations solved by FTCS.

    Consider the shallow water equations with rotation and bottom
    topography

                u_t + u * u_x + phi_x = f * v - g * H_x         (1a)
                        v_t + u * v_x = -f * u                  (1b)
        phi_t + u * phi_x + phi * u_x = 0                       (1c)

    where u and v are the eastward and northward wind of a fluid, and
    phi is the geopotential, i.e. phi(x, t) = g * h(x, t), where g is
    the acceleration due to gravity. These equations are the
    two-dimensional shallow water equations derived in [2] with the
    additional assumption that all functions are independent of y.
    Writing equations (1) in flux form, the equations become

    m_t + [m * u + (phi ** 2) / 2]_x = f * n - g * phi * H_x    (2a)
                     n_t + (m * v)_x = -f * m                   (2b)
                         phi_t + m_x = 0                        (2c)

    by letting m = phi * u and n = phi * v. Both [1] and [3] add
    diffusion terms to (2).

    Attributes:
        g: Acceleration due to gravity (m / s**2).
    
        J: Number of nodes in the spatial domain.
    
        K: Diffusion constant (m**2 / s).

        lat: Latitude in degrees.

        topo: String to determine whether the bottom topography will
              contain a ridge or be flat. The choices are "ridge" and
              "zero".

        f: Coriolis parameter for rotation (1 / s), computed from the
           latitude input.

        L: Circumference of latitude circle on the Earth (km),
           computed from the latitude input.

        n: Total number of variables.
    
        dx: Uniform spatial distance between variables.
    
        dt: Uniform time-step.

        H: Discretized bottom topography.

    Functionality:
        forecast: Numerical integration of the PDE system using the
            forward-time and centered-space finite difference scheme
            (FTCS).
    
        forecast_tlm: Tangent linear model of the FTCS integration.
    
        forecast_adj: Adjoint model of the FTCS integration.

        animate: Matplotlib animation of the shallow water equations
            using the initial conditions of [1].

        _topography: Bottom topography with a ridge [1].
    """
    def __init__(self, g=9.8, J=100, K=500.0, lat=30.0, topo="zero"):
        """Initializes model parameters.

        Griffith and Nichols [1] use the value of L corresponding to
        a latitude of 60 degrees and not 30 degrees. The possible
        values for topo are "zero", which means flat bottom
        topography, and "ridge", which uses _topography.
        """
        self.g = g
        self.J = J
        self.K = K

        lat_radians = math.radians(lat)
        self.f = 1.45842e-4 * math.sin(lat_radians)
        self.L = 40074.0 * math.cos(lat_radians)

        self.n = 3 * J
        self.dx = self.L / float(J)
        self.dt = 0.1 * self.dx

        if topo == "zero":
            self.H = np.zeros(J)
        elif topo == "ridge":
            x = np.arange(0.0, self.L, self.dx)
            self.H = np.empty(x.shape)
            
            for (index, x_value) in enumerate(x):
                self.H[index] = self._topography(x_value)
        else:
            raise ValueError("Choose a valid option.")

    def forecast(self, x, ndt=1):
        """Uses forward-time centered-space (FTCS) finite differences.

        Arguments:
            x: Initial state of the model trajectory
            ndt: Number of time-steps to integrate

        Returns:
            Model forecast of x.
        """
        J = self.J

        m0 = x[0 : J]
        n0 = x[J : 2 * J]
        phi0 = x[2 * J : 3 * J]
        
        (m, n, phi) = fortran.ftcs(m0, n0, phi0, self.H, self.f, \
                                   self.g, self.dx, self.dt, self.K)

        for i in xrange(1, ndt):
            (m, n, phi) = fortran.ftcs(m, n, phi, self.H, self.f, \
                                self.g, self.dx, self.dt, self.K)

        return np.hstack((m, n, phi))

    def forecast_tlm(self, x, xd, ndt=1):
        """Tangent linear model of FTCS.

        Arguments:
            x: The state at which the TLM is initialized
            xd: Direction vector
            ndt: Number of time-steps to integrate

        Returns:
            yd: The derivative of FTCS in the direction of xd.
        """
        J = self.J

        m0 = x[0 : J]
        n0 = x[J : 2 * J]
        phi0 = x[2 * J : 3 * J]
        m0_d = xd[0 : J]
        n0_d = xd[J : 2 * J]
        phi0_d = xd[2 * J : 3 * J]

        (md, nd, phid) = fortran.ftcs_tlm(m0, m0_d, n0, n0_d, \
            phi0, phi0_d, self.H, self.f, self.g, self.dx, \
            self.dt, self.K)

        for i in xrange(1, ndt):
            (m0, n0, phi0) = fortran.ftcs(m0, n0, phi0, self.H, self.f, \
                                   self.g, self.dx, self.dt, self.K)
            (md, nd, phid) = fortran.ftcs_tlm(m0, md, n0, nd, \
                phi0, phid, self.H, self.f, self.g, self.dx, \
                self.dt, self.K)

        return np.hstack((md, nd, phid))

    def forecast_adj(self, x, yb, ndt=1):
        """Adjoint model of FTCS.

        Arguments:
            x: The state at which the adjoint model is formed.
            yb: Direction vector
            ndt: Number of time-steps

        Returns:
            xb: The backward integration of yb.
        """
        xlist = [[]] * ndt
        xlist[0] = x
        J = self.J

        m0b = yb[0 : J]
        n0b = yb[J : 2 * J]
        phi0b = yb[2 * J : 3 * J]

        # Forward sweep
        for i in xrange(1, ndt):
            xlist[i] = self.forecast(xlist[i - 1])

        # Backward sweep
        for i in xrange(ndt - 1, -1, -1):
            m0 = xlist[i][0 : J]
            n0 = xlist[i][J : 2 * J]
            phi0 = xlist[i][2 * J : 3 * J]

            (m0b, n0b, phi0b) = fortran.ftcs_adj(m0, n0, phi0, m0b, \
                n0b, phi0b, self.H, self.f, self.g, self.dx, \
                self.dt, self.K)
        
        return np.hstack((m0b, n0b, phi0b))

    def animate(self):
        """Display animation of model variables.

        The initial state for the animation is the same that Griffith
        and Nichols [1] used.
        """        
        phi = 10.0 * np.ones(self.J)
        m = np.zeros(self.J)
        n = np.zeros(self.J)
        xaxis = np.arange(self.J)
        
        plt.clf()
        plt.axis([0, self.J - 1, -10.0, 15.0])
        plt.plot(xaxis, phi, ".-", label="$\phi$")
        plt.plot(xaxis, n, ".-", label="$n$")
        plt.plot(xaxis, m, ".-", label="$m$")
        plt.legend(loc="upper left", numpoints=1)
        plt.pause(0.0001)

        for k in xrange(1, 500 + 1):
            [m, n, phi] = fortran.ftcs(m, n, phi, self.H, self.f, \
                                   self.g, self.dx, self.dt, self.K)
        
            plt.clf()
            plt.axis([0, self.J - 1, -10.0, 15.0])
            plt.plot(xaxis, phi, ".-", label="$\phi$")
            plt.plot(xaxis, n, ".-", label="$n$")
            plt.plot(xaxis, m, ".-", label="$m$")
            plt.legend(loc="upper left", numpoints=1)
            plt.pause(0.0001)

        plt.show()
        
    def _topography(self, x):
        """A function describing the bottom topography.

        This function is the bottom topography as implemented in [1].
        This generates a ridge in the middle of the domain of width
        20 * dx. The function is
    
        H(x) = 0.5 * [1 - (x - L / 2) ** 2 / a ** 2],

        if |x - L / 2| <= a, with a chosen to be 10 as in [1].

        Argument:
            x: Spatial location in the domain.

        Returns:
            output: Topography at location x.
        """
        check = abs((x - self.L / 2.0) / (10.0 * self.dx))

        if check <= 1.0:
            output = 0.5 * (1.0 - check ** 2)
        else:
            output = 0.0
        
        return output
