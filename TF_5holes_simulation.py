import yaml
import numpy as np
from src.data import ThinFilmData

        
d = ThinFilmData("simulations/TF_5holes_simulation.yaml")
d.plot_initial_condition()

d.tfinal

error = 10
while error>.1:
    for i in range(1,100):
        error /= i
        print(error)
