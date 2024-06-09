import numpy as np
from src.data import ThinFilmData
from scipy.fft import ifft2, fft2
import matplotlib.pyplot as plt
        
d = ThinFilmData("simulations/TF_5holes_simulation.yaml")
d.plot_initial_condition()

class ThinFilmModel(object):
    def __init__(self) -> None:
        pass

    def phiU(self):
        pass
    
    def rhs1(self):
        pass
    
    def rhs2(self):
        pass
    
    
while d.t <  d.tfinal - d.dt*1*0-.0000001:
    print(f"time: {d.t}")
    # Just Eyres scheme - no initial guess      
    d.U1 = d.U.copy()
    
    for i in range(d.iter):
        
        d.MU = d.U1 * d.U1 * d.U1
        
        d.MU.shape
        d.MU
        phiU = (-d.eps2 / (d.U1**4)) * (3-(4*d.epsilon)/d.U1)+0.1
        gU = d.MU * phiU

        # |> RHS1 -------------------------
        
        LapU1=ifft2(-1*d.k2 * fft2(d.U1))
        lapU1_x = ifft2(1j * d.k1x * fft2(LapU1))
        lapU1_y = ifft2(1j * d.k1y * fft2(LapU1))
        
        lapU1_x = (d.M1 - d.MU) * lapU1_x
        lapU1_y = (d.M1 - d.MU) * lapU1_y

        rhs1 = ifft2(1j * d.k1x * fft2(lapU1_x)) + ifft2(1j * d.k1y * fft2(lapU1_y))
        rhs1 = np.real(rhs1)
        # |> RHS2 -------------------------

        f1 = ifft2(1j * d.k1x * fft2(d.U1))
        f2 = ifft2(1j * d.k1y * fft2(d.U1))
        # factor-in the M(u)
        f1=f1 * gU
        f2=f2 * gU
        
        rhs2 = ifft2(1j * d.k1x * fft2(f1)) + ifft2(1j * d.k1y * fft2(f2))
        rhs2 = np.real(rhs2)
      
        # |> -------------
        RHS = rhs1 + rhs2
        hat_rhs = d.hat_U + d.dt * fft2(RHS)
        hat_U1 = hat_rhs / d.lhs
        d.U1 = ifft2(hat_U1)
    
    d.U = d.U1
    d.hat_U=hat_U1
    d.it += 1
    d.t += d.dt
    
    # Energy Computations
    d.Ue1 = np.real(ifft2(-1j * d.k1x * fft2(d.U)))
    d.Ue2 = np.real(ifft2(-1j * d.k1y * fft2(d.U)))
    d.energy = -(d.eps2 / (d.U **2)) * ((1/2) - d.epsilon / (3 * d.U))+(1/2) * 0.1*d.U**2 + (1/2) * (d.Ue1**2 + d.Ue2**2)
    
    d.Energy.append(d.h * d.h * sum(sum(d.energy)))
    d.time.append(d.t)
    d.Mvar.append(np.max(np.max(d.U ** 3)))
    
    if d.Energy[-1] > d.Energy[int(d.it)]:
        Stability = 0
        break

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(d.X, d.Y, d.U, cmap=plt.cm.YlGnBu_r)
plt.title(f"t_F={np.round(d.t)}")
plt.show()

fig, ax = plt.subplots()
ax.plot(d.time, d.Mvar)
# ax = fig.add_subplot(projection='3d')
ax.set_xlabel("time")
ax.set_ylabel("Mvar")
plt.show()