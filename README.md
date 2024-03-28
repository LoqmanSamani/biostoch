
# BioStoch
<p><b>Simulation in Biology<p></p>

![ssa1.png](https://github.com/LoqmanSamani/biostoch/blob/systembiology/examples/plots/ssa1.png)
<p><i>Simulated with GillespySimulator: Concentration trajectories of species A and B over time. Ten independent simulations were performed, each represented by a different line color<i></p>


### Overview

- `biostoch` is a Python library for simulating chemical and biological models using various deterministic and stochastic methods. It provides implementations of methods such as Euler's method, the Runge-Kutta algorithm, the Stochastic Simulation Algorithm (SSA), Tau-Leaping, and the Chemical Langevin Equation (CLE). These methods can be used to model and analyze the dynamics of biochemical reactions, gene regulatory networks, and other biological systems.

### Article

- [Introducing BioStoch: A Python Library for Simulating Biological Systems](https://medium.com/@samaniloqman91/introducing-biostoch-a-python-library-for-simulating-biological-systems-5d38df9b6991) is the  article, which was published on Towards Data Science, provides an exploration of the mathematical bases behind various simulation methods used by BioStoch. From deterministic to stochastic algorithms, learn how to use the library to model biological systems.


### Installation

- Install `biostoch` and its dependencies using pip:

```bash
pip install numpy matplotlib 
pip install biostoch
```

### Usage
```python
import numpy
import matplotlib
import time
from biostoch.model import Model
from biostoch.ode import EulerSimulator, RungeKuttaSimulator
from biostoch.ssa import GillespieSimulator
from biostoch.tau_leaping import TauLeaping
from biostoch.cle import ChemicalLangevin
from biostoch.visualization import Visualization

# Define a system of reactions, in this case, a simple system with two reactions: A <-> B; with rate constants K1 = 0.1, K2 = 0.05
model = Model() 

# Add rate constants for each reaction in the system
model.parameters({
    "K1": 0.1, 
    "K2": 0.05
})

# Add species and rate of change for each of them
model.species(
    components={
        "A": 100.0, 
        "B": 0.0
    }, 
    rate_change={
        "A": "K2 * B - K1 * A", 
        "B": "K1 * A - K2 * B"
    }
)

# Add reactions and rate for each of them
model.reactions(
    reacts={
        "reaction1": "A -> B", 
        "reaction2": "B -> A"
    },
    rates={
        "reaction1": "K1 * A", 
        "reaction2": "K2 * B"
    }
)

# Simulate the model using ordinary differential equations (ODE) with the Euler method
euler_model = EulerSimulator(
    model=model, 
    start=0, 
    stop=100, 
    epochs=1000
)
euler_model.reset() # Reset the model to initialize if it has been used before
euler_model.simulate() # Simulate the model 
euler_model.species # Print the model species after the simulation, a dictionary containing the change in species concentration during the simulation time
euler_model.time # Show how long the simulation took to complete

# Simulate the model using ordinary differential equations (ODE) with the Runge-Kutta method
runge_model = RungeKuttaSimulator(
    model=model, 
    start=0, 
    stop=100, 
    epochs=1000
)

runge_model.simulate()
runge_model.species
runge_model.time

# Simulate the model using the Stochastic Simulation Algorithm (SSA)
ssa_model = GillespieSimulator(
    model=model, 
    start=0, 
    stop=100, 
    max_epochs=1000
)

ssa_model.simulate()
ssa_model.species
ssa_model.time

# Simulate the model using the Tau-Leaping method
tau_model = TauLeaping(
    model=model, 
    start=0, 
    stop=100, 
    max_epochs=100
)

tau_model.simulate()
tau_model.species
tau_model.time

# Simulate the model using the Chemical Langevin Equation method
cle_model = ChemicalLangevin(
    model=model,
    start=0,
    stop=100, 
    max_epochs=1000
)

cle_model.simulate()
cle_model.species
cle_model.time

# Visualize the simulated models using integrated matplotlib.pyplot in biostoch
euler_plot = Visualization(euler_model)
runge_plot = Visualization(runge_model)
ssa_plot = Visualization(ssa_model)
tau_plot = Visualization(tau_model)
cle_plot = Visualization(cle_model)

euler_plot.plot()
runge_plot.plot()
ssa_plot.plot()
tau_plot.plot()
cle_plot.plot()

```

### License

- This project is licensed under the [MIT License](https://github.com/LoqmanSamani/biostoch/blob/systembiology/LICENSE)
