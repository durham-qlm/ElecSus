import arc
import numpy as np
from scipy import constants as c
import scipy as sp
import sympy as sy
from sympy.physics.wigner import wigner_6j
from tqdm import tqdm


class state:
    def __init__(self, n, l, j, f=None):
        self.n = n
        self.l = l
        self.j = j
        self.f = f

    def __call__(self, precision):
        if precision == 'nlj':
            return (self.n, self.l, self.j)
        elif precision == 'nljf':
            return self.n, self.l, self.j, self.f

    def F(self, f):
        self.f = f
        return self


class atomicSystem:
    def __init__(self, element, states, T=20+273.15, beam_diameter=5e-3):
        if element.lower() in ['rb85', 'rubidium85']:
            self.atom = arc.Rubidium85()
            self.isotopeShift = 21.679e6
        elif element.lower() in ['rb87', 'rubidium87']:
            self.atom = arc.Rubidium87()
            self.isotopeShift = -56.219e6
        elif element.lower() in ['na', 'sodium', 'na23']:
            self.atom = arc.Sodium()
            self.isotopeShift = 0
        elif element.lower() in ['cs', 'caesium', 'cesium', 'cs23']:
            self.atom = arc.Caesium()
            self.isotopeShift = 0
        else:
            print('Could not parse desired element! Default to Rb85')
            self.atom = arc.Rubidium85()
            self.isotopeShift = 21.679e6

        import copy
        self.states = copy.deepcopy(states)
        self.max_allowed_states = 3
        self.n_states = len(states)
        if self.n_states > self.max_allowed_states:
            raise RuntimeError('Exceeded allowed number of states!')

        padding = [None] * (self.max_allowed_states - len(states))
        states.extend(padding)
        self.groundState = states[0]
        self.excitedState = states[1]
        self.rydbergState = states[2]

        self.T = T
        self.beam_diameter = beam_diameter

        self.initSystemProperties()
        self.generateSymbols()
        self.generateMatrices()
        dict = self.generateDictionaries()

        self.system_matrix = self.master_equation.subs(dict)
        # Add constrain that total population has to be 1
        self.system_matrix = self.system_matrix.as_mutable()
        self.system_matrix[0] = -1 + self.r.trace()
        self.A, self.b = self.generate_linear_system()

    def getSFF(self, state1, state2):
        # Ref: Steck, Daniel A. "Rubidium 85 D Line Data" (2009)
        f1 = np.atleast_1d(state1.f)
        f2 = np.atleast_1d(state2.f)
        SFF = np.zeros((f1.size, f2.size))
        for i, fi in enumerate(f1):
            for j, fj in enumerate(f2):
                SFF[i, j] = (2 * fj + 1) * (2 * state1.j + 1) \
                    * wigner_6j(state1.j, state2.j, 1, fj, fi, self.atom.I)**2
        return SFF

    def getF(self, state):
        # Number of Hyperfine states F:
        # |J-I| <= F <= J+I
        return np.arange(np.abs(self.atom.I - state.j),
                         self.atom.I + state.j + 1)

    def getBranchingRatio(self, state1, state2):
        # Ref: Wenting, Chen "Two-Photon spectroscopy of rubidium in the
        # vicinity of silicon devices" (2019)

        # State is assumed to be the lower energy state
        if self.atom.getEnergy(*state1('nlj')) > self.atom.getEnergy(*state2('nlj')):  # noqa
            state1, state2 = state2, state1
        f1 = np.atleast_1d(state1.f)
        f2 = np.atleast_1d(state2.f)
        B = np.zeros((f1.size, f2.size))
        for i, fi in enumerate(f1):
            for j, fj in enumerate(f2):
                B[i, j] = (2 * state2.j + 1) * (2 * state1.f + 1) \
                    * wigner_6j(state2.j, fj, self.atom.I, fi, state1.j, 1)**2
        return B

    def getTransitTime(self):
        # Ref: ARC-Alkali-Rydberg-Calculator Web interface (click 'View code')
        # in s
        mean_speed = self.atom.getAverageSpeed(self.T)
        # 2w(1/e2) = 1.699 * FWHM
        beam_fwhm = self.beam_diameter / 1.699
        return np.sqrt(c.pi)*beam_fwhm/(2.*mean_speed)

    def getEnergySeparation(self, state):
        sublevels = self.getF(state)
        energies = np.zeros_like(sublevels)
        try:
            A, B = self.atom.getHFSCoefficients(*state('nlj'))
            for i, f in enumerate(sublevels):
                # Get HFS magnetic dipole/quadrupole constant
                energies[i] = self.atom.getHFSEnergyShift(state.j, f, A, B)
        except:
            # we simply return zeros here
            print('Detected high energy state!')
            energies = np.arange(sublevels.size) * 1e2
        return energies

    def initSystemProperties(self):
        self.f_resonance = np.array([self.atom.getTransitionFrequency(
            *self.states[i]('nlj'), *self.states[i+1]('nlj'))
            for i in range(self.n_states-1)])

        self.sublevels = [self.getF(state) for state in self.states]
        self.n = np.array([len(level) for level in self.sublevels])
        self.total_levels = self.n.sum()

        self.transit_time = self.getTransitTime()
        self.energySeparation = [self.getEnergySeparation(state) for state in self.states]
        self.energySeparation[1] -= self.isotopeShift
        self.slices = [slice(self.n[0:i].sum(), self.n[0:i+1].sum()) for i in range(self.n_states)]
        # dipole matrix moment: |<J||er||J'>|
        # SFF'(F->F') = (2F'+1)(2J+1){J, J', 1, F', F, I}^2
        # d^2 = 1/3 * SFF' * |<J||er||J'>|^2
        # Ω = d * E / hbar = sqrt(SFF/3) * |<J||er||J'>| * E / hbar
        # Ref: Steck, Daniel A. "Rubidium 85 D Line Data" (2009), p. 10
        DME = [self.atom.getReducedMatrixElementJ_asymmetric(
            *self.states[i]('nlj'), *self.states[i+1]('nlj'))
            * c.e * c.physical_constants['Bohr radius'][0]
            for i in range(self.n_states-1)]
        SFF = [self.getSFF(self.states[i].F(self.sublevels[i]),
                           self.states[i+1].F(self.sublevels[i+1]))
               for i in range(self.n_states-1)]
        self.dipole_moments = [np.sqrt(SFF[i].ravel()/3) * DME[i]
                               for i in range(self.n_states-1)]
        self.Gammas = [self.atom.getTransitionRate(
            *self.states[i+1]('nlj'), *self.states[i]('nlj')) / 2 / c.pi
            for i in range(self.n_states-1)]

    def generateSymbols(self):
        #######################################################################
        # Generate symbols and variables
        #######################################################################
        # Symbols for both laser frequencies
        self.wL = np.array([sy.symbols(f'w_{i}{i+1}') for i in range(self.n_states-1)])
        self.E = [np.array(sy.symbols(f'E_{self.n[0:i].sum()}:{self.n[0:i+1].sum()}')) for i in range(self.n_states)]

        self.O = [np.array(sy.symbols(
            f'O_({self.n[0:i].sum()}:{self.n[0:i+1].sum()})({self.n[0:i+1].sum()}:{self.n[0:i+2].sum()})'
            )).reshape((self.n[i], self.n[i+1]))
            for i in range(self.n_states-1)]

        self.G_01, self.G_10 = sy.symbols('G_01 G_10')
        self.G = [np.array(sy.symbols(
            f'G_({self.n[0:i+1].sum()}:{self.n[0:i+2].sum()})({self.n[0:i].sum()}:{self.n[0:i+1].sum()})'
            )).reshape((self.n[i+1], self.n[i]))
            for i in range(self.n_states-1)]

        self.r_individual = sy.symbols(
            f'\\rho_{{(0:{self.total_levels})(0:{self.total_levels})}}')
        self.r = sy.Matrix(self.total_levels,
                           self.total_levels, self.r_individual)

    def generateMatrices(self):
        #######################################################################
        # Generate matrices
        #######################################################################
        self.H_rabi = sy.zeros(self.total_levels, self.total_levels)
        for i in range(self.n_states-1):
            self.H_rabi[self.slices[i], self.slices[i+1]] = 0.5 * self.O[i]
        self.H_rabi = self.H_rabi + sy.conjugate(self.H_rabi.T)

        detunings = np.concatenate([-self.E[i]+self.wL[0:i].sum()
                                    for i in range(self.n_states)], axis=None)
        self.H_energylevels = sy.diag(*detunings)
        self.H = self.H_rabi + self.H_energylevels

        G = sy.zeros(self.total_levels, self.total_levels)
        G[0, 1] = self.G_01
        G[1, 0] = self.G_10
        for i in range(self.n_states-1):
            G[self.slices[i+1], self.slices[i]] = self.G[i]

        # Ref: Weller, PhD thesis, p. 14, eq. 1.13
        L = sy.zeros(self.total_levels, self.total_levels)
        for i in range(self.total_levels):
            for j in range(self.total_levels):
                L[i, i] += G[j, i] * self.r[j, j] - G[i, j] * self.r[i, i]
                if (i != j):
                    for k in range(self.total_levels):
                        L[i, j] -= 0.5 * (G[i, k] + G[j, k]) * self.r[i, j]

        self.master_equation = -sy.I * (self.H*self.r - self.r*self.H) - L

    def generate_linear_system(self):
        self.r_list = self.matrix2list(self.r)
        # Create list of off-diagonal elements relevant for i->j transition
        self.transition_list = []
        for i in range(self.n_states-1):
            mask = np.full((self.total_levels, self.total_levels), False)
            mask[self.slices[i], self.slices[i+1]] = True
            self.transition_list.append(self.matrix2list(mask))

        A, b = sy.linear_eq_to_matrix(self.system_matrix.expand(), self.r_list)
        A = sy.lambdify([*self.wL, *np.concatenate(self.O, axis=None)], A, 'numpy')
        # b = sy.lambdify([self.wL[0],
        #                  self.wL[1],
        #                  *self.O[0].flatten(),
        #                  *self.O[1].flatten()], b, 'numpy')
        b = np.zeros((self.total_levels**2, 1))
        b[0] = 1
        return A, b

    def generateDictionaries(self):
        n_mf1 = 2 * self.sublevels[0][0] + 1
        n_mf2 = 2 * self.sublevels[0][1] + 1
        dict_transit = {
            self.G_01: n_mf2 / (n_mf1 + n_mf2) / self.transit_time,
            self.G_10: n_mf1 / (n_mf1 + n_mf2) / self.transit_time,
        }

        dict_G = {}
        for h in range(self.n_states-1):
            for i, f1 in enumerate(self.sublevels[h]):
                for j, f2 in enumerate(self.sublevels[h+1]):
                    B = self.getBranchingRatio(
                        self.states[h+1].F(f2), self.states[h].F(f1))
                    dict_G[self.G[h][j, i]] = float(B * self.Gammas[h])

        dict_E = [{key: val
            for (key, val) in zip(self.E[i], self.energySeparation[i])}
            for i in range(self.n_states)]

        dict_all = {**dict_transit, **dict_G}
        for i in range(self.n_states):
            dict_all = {**dict_all, **dict_E[i]}
        return dict_all

    def printStats(self):
        print(f'============= {self.atom.elementName} =============')
        for i, f in enumerate(self.f_resonance):
            print(f'{i}->{i+1}: {f/1e12:.3f} THz / {c.c / f * 1e9:.2f} nm')
        for i, SL in enumerate(self.sublevels):
            print(f'F of state {i}: {SL}')
        print(f'Transit time: {self.transit_time * 1e6:.3f} μs')
        for i, G in enumerate(self.Gammas):
            print(f'Γ_{i}{i+1}: 2π * {G/1e6:.3f} MHz')

        # print(f'Γ_ge: 2π * {self.Gammas_ge/1e6:.3f} MHz')
        print(f'Branching ratios e->g:')
        for F_p in self.sublevels[1]:
            for F in self.sublevels[0]:
                B = self.getBranchingRatio(
                    self.excitedState.F(F_p), self.groundState.F(F))
                print(f"F'={F_p} -> F={F}: B={float(B):.3f}")
        print(f'==================================================')

    def v_dist(self, v):
        return np.sqrt(self.atom.mass / (2 * c.pi * c.k * self.T)) \
            * np.exp(-self.atom.mass * v**2 / (2 * c.k * self.T))

    def cdf(self, v):
        o = np.sqrt(c.k * self.T / self.atom.mass) * np.sqrt(2)
        return 0.5 * (1 + sp.special.erf(v/o))

    def cdfinv(self, p):
        o = np.sqrt(c.k * self.T / self.atom.mass) * np.sqrt(2)
        return o * sp.special.erfinv(2 * p - 1)

    def matrix2list(self, mat):
        # Generate list containing all entries of density matrix
        # First diagonal entries: r[0,0], r[1,1], ...
        # Then upper half: r[0,1], r[0,2], ...
        # Then lower half: r[1,0], r[2,0], ...
        l1 = []
        l2 = []
        for i in range(mat.shape[0]):
            for j in range(i+1, mat.shape[1]):
                l1.append(mat[i, j])
                l2.append(mat[j, i])
        return list(mat.diagonal()) + l1 + l2

    def unpack_beam(self, beam):
        if len(beam) == 3:
            w, P, D = beam
            sgn = 1
        elif len(beam) == 4:
            w, P, D, sgn = beam
        return w, P, D, sgn

    def beam2rabi_Efield(self, beam, type='ge'):
        # sympy matrix is constructed with 2*pi*w "space"
        w, P, D, _ = self.unpack_beam(beam)
        w = np.atleast_1d(w)
        A = c.pi * (D / 2) ** 2
        # I = 1/2 * c * epsilon_0 * E0**2
        I = 2 * P / A
        E = np.atleast_1d(np.sqrt(2 * I / c.c / c.epsilon_0))
        if type == 'ge':
            dipole_moments = self.dipole_moments[0]
        elif type == 'er':
            dipole_moments = self.dipole_moments[1]
        rabi = np.outer(E, dipole_moments) / c.hbar
        return rabi, E

    def solve(self, beam_ge, beam_er=None, v=None):
        #######################################################################
        # Calculate Rabi Frequencies
        #######################################################################
        w_ge, _, _, sgn_ge = self.unpack_beam(beam_ge)
        w_ge = np.atleast_1d(w_ge)
        rabi_ge, E_ge = self.beam2rabi_Efield(beam_ge)
        wavenumber_ge = self.f_resonance[0] / c.c

        if beam_er is not None:
            w_er, _, _, sgn_er = self.unpack_beam(beam_er)
            w_er = np.atleast_1d(w_er)
            rabi_er, E_er = self.beam2rabi_Efield(beam_er)
            wavenumber_er = self.f_resonance[1] / c.c

        #######################################################################
        # Solve linear system
        #######################################################################
        if self.rydbergState is None:
            A = [[self.A(w, *o) for o in rabi_ge] for w in w_ge]
            b = np.zeros((self.total_levels**2, 1))
            b[0] = 1
        else:
            if v is None:
                A = [[[[self.A(w1, w2, *o1, *o2) for o2 in rabi_er]
                    for o1 in rabi_ge] for w2 in w_er] for w1 in w_ge]
                b = [[[[self.b(w1, w2, *o1, *o2) for o2 in rabi_er]
                    for o1 in rabi_ge] for w2 in w_er] for w1 in w_ge]
            else:
                w_ge = w_ge[0]
                w_er = w_er[0]
                A = [[[self.A(w_ge + sgn_ge * wavenumber_ge * vi,
                              w_er + sgn_er * wavenumber_er * vi,
                              *o1, *o2)
                       for o2 in rabi_er]
                      for o1 in rabi_ge]
                     for vi in v]
                b = [[[self.b(w_ge + sgn_ge * wavenumber_ge * vi,
                              w_er + sgn_er * wavenumber_er * vi,
                              *o1, *o2)
                       for o2 in rabi_er]
                      for o1 in rabi_ge]
                     for vi in v]
        res = np.linalg.solve(A, b)

        #######################################################################
        # Extract relevant information
        #######################################################################
        # Move density matrix dimension to the front
        res = np.moveaxis(res.squeeze(), -1, 0)
        k_ge = np.divide.outer(self.dipole_moments[0], E_ge) / c.epsilon_0
        # - Return sum of excited states. Indexing order is given by order
        #   of arguments of 'sy.linear_eq_to_matrix' above
        # - Population of excited states is given by diagonal entries
        # - Complex-valued susceptibility is given by off-axis entries
        #   (only one side, they are complex conjugated anyway)
        exci_state = np.sum(res[self.slices[1]], axis=0).real
        chi_ge = self.atom.abundance * 2 \
            * np.sum(res[self.transition_list[0]] * k_ge, axis=0)
        # I don't know why, but currently this fixes to correct absorption
        chi_ge = chi_ge / 2 / c.pi
        return exci_state.squeeze(), chi_ge.squeeze()

    def solve_w_doppler(self, beam_ge, beam_er=None):
        # chi_dopp(∆) = \int p(v,T)chi(∆-kv)dv = (p*chi)(∆)
        # k = 1 / lam2bda = w/c
        w_ge, P_ge, D_ge, sgn_ge = self.unpack_beam(beam_ge)
        w_ge = np.atleast_1d(w_ge)
        P_ge = np.atleast_1d(P_ge)
        k_ge = self.f_resonance[0] / c.c

        if beam_er is not None:
            w_er, P_er, D_er, sgn_er = self.unpack_beam(beam_er)
            w_er = np.atleast_1d(w_er)
            P_er = np.atleast_1d(P_er)
            k_er = self.f_resonance[1] / c.c

        if beam_er is None:
            exci_state = np.ones((len(w_ge), len(P_ge)), dtype='float64')
        else:
            exci_state = np.zeros((len(w_ge), len(w_er), len(P_ge)), dtype='float64')
        chi = np.zeros_like(exci_state, dtype='complex128')

        for i, P in enumerate(P_ge):
            if beam_er is None:
                resolution = 2  # MHz
                v = np.linspace(
                    w_ge.min()/1e6 - 1000,
                    w_ge.max()/1e6 + 1000,
                    int((w_ge.ptp()/1e6 + 2000) / resolution))
                dv = v[1] - v[0]
                v_distribution = self.v_dist(np.subtract.outer(w_ge / k_ge, v))
                # Use symmetry of convolution to calculate population_number
                # once and instead calculate Maxwell Boltzmann distribution
                # more often (which is computationally cheaper)
                E, C = self.solve((k_ge * v, P, D_ge, 1))
                exci_state[:, i] = np.sum(v_distribution * E, axis=1) * dv
                chi[:, i] = np.sum(v_distribution * C, axis=1) * dv
            else:
                # Create nonequidistant velocity distribution
                # with equal bin probabilities. See nonequidistant_sampling.py
                # for more.
                binsize = 10000
                epsilon = 1e-6
                p_bounds = np.linspace(epsilon, 1-epsilon, binsize + 1)
                p = (p_bounds[1:] + p_bounds[:-1]) / 2
                v_bounds = self.cdfinv(p_bounds)
                v = self.cdfinv(p)
                dv = np.diff(v_bounds)
                v_distribution = self.v_dist(v)

                for m, wm in enumerate(tqdm(w_er, position=1)):
                    for n, wn in enumerate(tqdm(w_ge, position=0, leave=False)):
                        E, C = self.solve((wn, P, D_ge, sgn_ge),
                                          (wm, P_er, D_er, sgn_er), v)
                        exci_state[n, m, i] = np.sum(E * dv * v_distribution)
                        chi[n, m, i] = np.sum(C * dv * v_distribution)
        return exci_state.squeeze(), chi.squeeze()

    def transmission(self, beam_ge, beam_er=None, z=50e-3, doppler=True):
        alpha = self.optical_depth(beam_ge, beam_er, doppler)
        return np.exp(alpha * z)

    def absorption(self, beam_ge, beam_er=None, z=50e-3, doppler=True):
        return 1 - self.transmission(beam_ge, beam_er, z, doppler)

    def optical_depth(self, beam_ge, beam_er=None, doppler=True):
        n = self.atom.getNumberDensity(self.T)
        if doppler:
            _, chi = self.solve_w_doppler(beam_ge, beam_er)
        else:
            _, chi = self.solve(beam_ge, beam_er)
        n_imag = np.sqrt(1.0 + chi * n).imag
        return 4 * c.pi * self.f_resonance[0] / c.c * n_imag


if __name__ == '__main__':
    groundState = state(5, 0, 1/2)     # 5S1/2
    excitedState = state(5, 1, 1/2)    # 5P1/2
    rydbergState = state(17, 0, 1/2)  # 17S1/2
    rb85 = atomicSystem(
        'Rb85',
        [groundState,
        excitedState,
        rydbergState],
        T=20 + 273.15,
        beam_diameter=3e-3)

    rb85.printStats()

    # # Calculation of C_p via C_mf1mf2
    # # C_mf1mf2 is given by combination of
    # # Steck, eq. (35) and (36) as well as
    # # Chen, eq. 38.
    # J = 1/2
    # Jp = 3/2
    # F = 2
    # Fp = 3
    # q = +1

    # C_q = 0
    # for mf1 in np.arange(-F, F+1):
    #     for mf2 in np.arange(-Fp, Fp+1):
    #         C_mf1mf2 = np.sqrt(2*F+1)\
    #             * arc.wigner.Wigner3j(Fp, 1, F, mf2, q, mf1)\
    #             * np.sqrt((2*Fp+1)*(2*J+1))\
    #             * arc.wigner.Wigner6j(J, Jp, 1, Fp, F, 5/2)
    #         if mf2 == mf1 + q:
    #             print(mf1, mf2, C_mf1mf2)
    #         C_q += C_mf1mf2**2
    # C_q /= (2*F + 1)
    # print(C_q)
