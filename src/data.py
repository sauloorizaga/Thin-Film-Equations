import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, ifft2
from src.utils import ConfigReader

class ThinFilmData(ConfigReader):
    def __init__(self, config_path: str) -> None:
        ConfigReader.__init__(self, config_path)
        self.data = self._read_yaml()
        self._set_simulation_attrs()
        self._set_more_attrs()
        self._define_lhs()
        self._energy_computations()

    def gaussian_bump(self, center_x, center_y, coefficient=0.5):
        """
        This function defines a 2D Gaussian bump centered at (center_x, center_y).

        Args:
            center_x: x-coordinate of the center.
            center_y: y-coordinate of the center.
            coefficient: Scales the amplitude of the bump (default is 0.5).

        Returns:
            A 2D numpy array representing the Gaussian bump.
        """
        return coefficient * np.exp(-0.5 * (np.power(self.X - center_x, 2) + np.power(self.Y - center_y, 2)))
    
    def set_initial_condition(self):
        """
        Defines the initial condition of the wave simulation by creating a wave-like profile 
        using mathematical expressions.
        """
        z1 =  self.gaussian_bump(2*np.pi, 1.5*np.pi)
        z2 =  self.gaussian_bump(4*np.pi, 1.5*np.pi)
        z3 = -self.gaussian_bump(2*np.pi, 4*np.pi)
        z4 =  self.gaussian_bump(4*np.pi, 4*np.pi)
        z5 = -self.gaussian_bump(3*np.pi, 4.5*np.pi)
        
        # TODO: dont hardcode the .6
        self.U = self.init_U_z_shift - z1 - z2 + z3 - z4 + z5
        
    def plot_initial_condition(self):
        """
        Plots the initial condition of the wave using the matplotlib library, 
        if the initial condition (self.U) has already been set.
        """
        if self.U is not None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            ax.plot_surface(self.X, self.Y, self.U, cmap=plt.cm.YlGnBu_r)
            plt.show()
            
    def _set_more_attrs(self):
        """
        Sets various attributes likely related to the simulation grid and wave properties, 
        including:
        - b (related to domain size)
        - h (grid spacing)
        - N (number of grid points)
        - x (grid positions in x-direction)
        - X, Y (meshed grids in x and y directions)
        - k (wavenumbers)
        - k1x, k1y, kx, ky (derived wavenumber related quantities)
        - k2, k4 (derived wavenumber related quantities)
        - eps2 (possibly related to wave properties)
        """
        # self.b
        self.b = self.M * np.pi
        
        # %uniform mesh thickness
        self.h=(self.b-self.a)/self.N
        
        # (Periodic bdy conditions)
        self.n = self.N
        # xgrid formtation (a b] and eventually (a b]^2
        self.x = np.arange(self.a + self.h, self.b + self.h, self.h)
        
        self.X, self.Y = np.meshgrid(self.x, self.x)
        
        # wave number generation (same as in 1d)
        self.k = np.concatenate((np.arange(self.N // 2+1), -np.arange(self.N // 2-1, 0, -1)))/(self.M//2)
        
        self.k1x, self.k1y = np.meshgrid(np.power(self.k, 1), np.power(self.k, 1))
        self.kx, self.ky = np.meshgrid(np.power(self.k, 2), np.power(self.k, 2))
        
        self.k2 = self.kx + self.ky
        self.k4 = np.power(self.k2, 2)
        
        self.set_initial_condition()
        # self.U
        
        self.eps2 = self.epsilon**2
        
        # self.Y.size
        self.Mvar = [np.max(np.max(self.U ** 3))]    
    
    def _define_lhs(self):
        # The LHS, left hand side of problem----------
        self.lhs = 1 + self.dt*self.M1*self.k4
        
    def _energy_computations(self):
        self.hat_U = fft2(self.U)    
        self.Ue1 = np.real(ifft2(-1j*self.k1x * fft2(self.U)))
        self.Ue2 = np.real(ifft2(-1j*self.k1y * fft2(self.U)))
        
        self.energy=-(self.eps2/(self.U**2)) * ((1/2)-self.epsilon/(3*self.U))+(1/2)*0.1*self.U**2+(1/2)*( self.Ue1**2+self.Ue2**2)
        self.Energy = [self.h*self.h*sum(sum(self.energy))]
        self.time = [0]