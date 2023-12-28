#!/usr/bin/env python
# coding: utf-8

# ODE Simulation of Simple Chemical Reaction System:
# 
#  - simulate a simple system using ordinary differential equations(ODE) with just two   species 
# 
# The system:
#      - Forward_reaction =  A --> B,  rate_constant = k1, rate_equation = k1 * A
#      - Reverse_reaction =  B --> A,  rate_constant = k2, rate_equation = k2 * B

# In[30]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import tellurium as te
import gillespy2 



# Using numpy, matplotlib 
def ode_simulator(sp1_0, sp2_0, k_1, k_2, t_0, t_end, num_steps):
     
    # Initializing Arrays for Concentration and Time Tracking
    simulation_time = np.linspace(t_0, t_end, num_steps)
    sp1_conc = np.zeros(num_steps)
    sp2_conc = np.zeros(num_steps)

    # Initialize the system
    sp1_conc[0] = sp1_0
    sp2_conc[0] = sp2_0

    # Simulate the system
    for i in range(1, num_steps):
        Rate_of_change_sp1 = k_2 * sp2_conc[i-1] - k_1 * sp1_conc[i-1]
        Rate_of_change_sp2 = k_1 * sp1_conc[i-1] - k_2 * sp2_conc[i-1]
        
        # Update the concentrations
        sp1_conc[i] = sp1_conc[i-1] + Rate_of_change_sp1 * (simulation_time[i] - simulation_time[i-1])
        sp2_conc[i] = sp2_conc[i-1] + Rate_of_change_sp2 * (simulation_time[i] - simulation_time[i-1])
    
    plt.figure(figsize=(8,6))
    plt.plot(simulation_time, sp1_conc, label = 'sp1')
    plt.plot(simulation_time, sp2_conc, label = 'sp2')
    plt.xlabel('Time', fontsize=10)
    plt.ylabel('Concentration of Species', fontsize=10)
    plt.title('Deterministic Simulation of a System using ODEs', fontsize=10)
    plt.legend()
    
    return plt.show()
    
    
k_1 = 0.1
k_2 = 0.05
sp1_0 = 100
sp2_0 = 0
t_0 = 0
t_end = 100
num_steps = 1000

f = ode_simulator(sp1_0, sp2_0, k_1, k_2, t_0, t_end, num_steps)    


# In[27]:


# Using numpy, matplotlib and scipy
def ode_simulator(sp1_0, sp2_0, k_1, k_2, t_0, t_end, num_steps):
    
    
    # Generate vectors for both species 
    def reaction_system(y, t):
        sp1, sp2 = y
        dsp1_dt = -k_1 * sp1 + k_2 * sp2
        dsp2_dt =  k_1 * sp1 - k_2 * sp2
        return [dsp1_dt, dsp2_dt]

    
    # Initialization
    y_0 = [sp1_0, sp2_0]
    simulation_time = np.linspace(t_0, t_end, num_steps) 
    
    # Using odeint from scipy to generate the populations
    result = odeint(reaction_system, y_0, simulation_time)
    
    # Separate the generated populations from each other
    pop_sp1 = result[:, 0]
    pop_sp2 = result[:, 1]
    
    plt.figure(figsize=(8,6))
    plt.plot(simulation_time, pop_sp1, label = 'sp1')
    plt.plot(simulation_time, pop_sp2, label = 'sp2')
    plt.xlabel('Time', fontsize=10)
    plt.ylabel('Concentration of Species', fontsize=10)
    plt.title('Deterministic Simulation of a System using ODEs', fontsize=10)
    plt.legend()
    
    return plt.show()

model = ode_simulator(sp1_0, sp2_0, k_1, k_2, t_0, t_end, num_steps)    


# In[28]:


# Using matplotlib and tellurium
def ode_simulator(sp1_0, sp2_0, k_1, k_2, t_0, t_end, num_steps):
    
    # Using antimony generate the system of reactions
    model= te.loada('''
      k1=0.1; k2=0.05; A=100; B=0;
      A -> B; k1*A;
      B -> A; k2*B;
      ''')
    # Make sure to reset the model
    model.reset()
    # Run the simulation
    result = model.simulate(t_0, t_end, num_steps)


    # Separate the generated populations from each other
    simulation_time = result[:, 0]
    pop_sp1 = result[:, 1]
    pop_sp2 = result[:, 2]
    
    plt.figure(figsize=(8,6))
    plt.plot(simulation_time, pop_sp1, label = 'sp1')
    plt.plot(simulation_time, pop_sp2, label = 'sp2')
    plt.xlabel('Time', fontsize=10)
    plt.ylabel('Concentration of Species', fontsize=10)
    plt.title('Deterministic Simulation of a System using ODEs', fontsize=10)
    plt.legend()
    
    return plt.show()

model = ode_simulator(sp1_0, sp2_0, k_1, k_2, t_0, t_end, num_steps)    


# In[36]:


# Using gillespy2
def ode_simulator(k_1, k_2, sp1_0, sp2_0, t_0, t_end, num_steps):
    
    model = gillespy2.Model()

    # Define parameters 
    k_1 = gillespy2.Parameter(name='k_1', expression=k_1)
    k_2 = gillespy2.Parameter(name='k_2', expression=k_2)
    model.add_parameter([k_1, k_2])

    # Define species
    sp1 = gillespy2.Species(name='reactant', initial_value=sp1_0)
    sp2 = gillespy2.Species(name='product',   initial_value=sp2_0)
    model.add_species([sp1, sp2])

    # Define the reactions
    reaction_1 = gillespy2.Reaction(name="r_forward", rate=k_1, reactants={sp1:1}, products={sp2:1})
    reaction_2 = gillespy2.Reaction(name="r_backward", rate=k_2, reactants={sp2:1}, products={sp1:1})
    model.add_reaction([reaction_1, reaction_2])

    # Set the timespan for the simulation.
    tspan = gillespy2.TimeSpan.linspace(t=t_end, num_points=num_steps)
    model.timespan(tspan)
    
    # Choose the model to run(ode)
    results = model.run(algorithm ='ODE') 
    
    p = results['product']
    t = results['time']
    r = results['reactant']
    plt.plot(t, r, label='A')
    plt.plot(t, p, label='B')
    plt.title('Deterministic Simulation of a System using ODE(with gillespy2)', fontsize=10)
    plt.xlabel('Time', fontsize=10)
    plt.ylabel('Concentration os Species', fontsize=10)
    plt.legend()
    
    return plt.show()

model = ode_simulator(0.1, 0.05, 100, 0, 0, 100, 1000)



