import gillespy2 as g
import numpy as np
import matplotlib.pyplot as plt


def reaction_system(X_0, algorithm, num_steps):
    
    model = g.Model(name='reaction_network')
    
    # rate constants
    k1 = g.Parameter(name='k1', expression=0.0000003)
    k2 = g.Parameter(name='k2', expression=0.0001)
    k3 = g.Parameter(name='k3', expression=0.001)
    k4 = g.Parameter(name='k4', expression=3.50)
    
    # species
    A = g.Species(name='A', initial_value=100000)
    B = g.Species(name='B', initial_value=200000)
    X = g.Species(name='X', initial_value=X_0)

    # Propensity functions 
    a1 = "k1 * A * X * (X - 1) / 2"
    a2 = "k2 * X * (X - 1) * (X - 2) / 6"
    a3 = "k3 * B"
    a4 = "k4 * X"

    # reactions
    r1 = g.Reaction(name='r1', propensity_function=a1,  reactants={A:1, X:2}, products={X:3, A:1})
    r2 = g.Reaction(name='r2', propensity_function=a2,  reactants={X:3, A:1}, products={A:1, X:2})
    r3 = g.Reaction(name='r3', propensity_function=a3,  reactants={B:1}, products={X:1, B:1})
    r4 = g.Reaction(name='r4', propensity_function=a4,  reactants={X:1, B:1}, products={B:1})

    # add rate constants, species and reactions to the model
    model.add_parameter([k1, k2, k3, k4])
    model.add_species([A, B, X])
    model.add_reaction([r1, r2, r3, r4])

    # define simulations time
    tspan = g.TimeSpan.linspace(t=100, num_points=num_steps)
    model.timespan(tspan)
    # simulate the model
    result = model.run(algorithm = algorithm)
    
    return result

algorithms = ['ODE','ODE','ODE','SSA','SSA','SSA','SSA','SSA','SSA','SSA','SSA','SSA',
              'SSA','SSA','SSA','SSA','SSA','SSA','SSA','SSA','SSA','SSA','SSA','SSA',
              'SSA','SSA','SSA','SSA','SSA','SSA']
x_zeros = [100,250, 600,100,250,200,250,250,250,250,250,250,250,250,250,400,600,250,200,
           200,100,400,600,250,250,200,400,600,250,250]

plt.figure(figsize=(13,7))
for i in range(len(algorithms)):
    result = reaction_system(x_zeros[i], algorithms[i],1000)
    plt.plot(result['time'], result['X'])
plt.xlabel('Time', fontsize=12)
plt.ylabel('X Concentration', fontsize=12)
plt.title('Simulation of The Schlögl Model', fontsize=14)
plt.xlim(0, 100)
plt.show()

