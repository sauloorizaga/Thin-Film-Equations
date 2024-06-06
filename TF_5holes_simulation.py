import yaml
import numpy as np
from src.utils import ConfigReader
import matplotlib.pyplot as plt
import scipy

class Data(ConfigReader):
    def __init__(self, config_path: str) -> None:
        ConfigReader.__init__(self, config_path)
        self.data = self._read_yaml()
        self._set_simulation_attrs()
        self._set_more_attrs()
        self.define_lhs()

    
    def set_initial_condition(self):
        z1 = .5 * np.exp(-.5 * (np.power(self.X - 2 * np.pi, 2) + 1 * np.power(self.Y - 1.5 * np.pi, 2)))
        z2 = .5 * np.exp(-.5 * (np.power(self.X - 4 * np.pi, 2) + 1 * np.power(self.Y - 1.5 * np.pi, 2)))
        z3 = -.5 * np.exp(-.5 * (np.power(self.X - 2 * np.pi, 2) + 1 * np.power(self.Y - 4 * np.pi, 2)))
        z4 = .5 * np.exp(-.5 * (np.power(self.X - 4 * np.pi, 2) + 1 * np.power(self.Y - 4 * np.pi, 2)))
        z5 = -.5 * np.exp(-.5 * (np.power(self.X - 3 * np.pi, 2) + 1 * np.power(self.Y - 4.5 * np.pi, 2)))
        
        self.U = .6 - z1 - z2 + z3 - z4 + z5
        
    def plot_initial_condition(self):
        if self.U is not None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            ax.plot_surface(self.X, self.Y, self.U, cmap=plt.cm.YlGnBu_r)
            plt.show()
            
    
    def _set_more_attrs(self):
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
    def define_lhs(self):
        
        # The LHS, left hand side of problem----------
        self.lhs = 1 + self.dt*self.M1*self.k4
        self.lhs.shape
        scipy.fft.fft2
        
        pass
data = self = Data("simulations/TF_5holes_simulation.yaml")

type(data.N)
data.lhs



