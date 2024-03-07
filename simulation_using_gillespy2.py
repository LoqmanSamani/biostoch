
"""
Simulate a simple system using different simulation methods in gillespy2 package with just two   species 

The system:
    - Forward_reaction =  A --> B,  rate_constant = k1, rate_equation = k1 * A
    - Reverse_reaction =  B --> A,  rate_constant = k2, rate_equation = k2 * B
"""


import matplotlib.pyplot as plt
import gillespy2 as gi





class GillespySimulator:
    
    def __init__(self, k_1, k_2, sp1_0, sp2_0, t_0, t_end, num_steps, method):
        
        self.k_1 = k_1
        self.k_2 = k_2
        self.sp1_0 = sp1_0
        self.sp2_0 = sp2_0
        self.t_0 = t_0
        self.t_end = t_end
        self.num_steps = num_steps
        self.method = method


        
    def simulation(self):
            
        # Define the simulation system
        model = gi.Model()
    
        # Define the rate constant for each reaction
        rate_constant_1 = gi.Parameter(name='k_1', expression = self.k_1)
        rate_constant_2 = gi.Parameter(name='k_2', expression = self.k_2)
   
        # Define the species in each reaction
        species_1 = gi.Species(name = 'reactant', initial_value = self.sp1_0)
        species_2 = gi.Species(name = 'product', initial_value = self.sp2_0)
   
        # Define the reactions
        reaction_1 = gi.Reaction(name = "forward", rate = rate_constant_1, reactants = {species_1:1}, products = {species_2:1})
        reaction_2 = gi.Reaction(name = "backward", rate = rate_constant_2, reactants = {species_2:1}, products = {species_2:1})
    
        # Define the time of the simulation
        simulation_time = gi.TimeSpan.linspace(t = self.t_end, num_points = self.num_steps)
    
    
        # Add rate constants, species, reactions and the simulation time to the system
        model.add_parameter([rate_constant_1, rate_constant_2])
        model.add_species([species_1, species_2])
        model.add_reaction([reaction_1, reaction_2])
        model.timespan(simulation_time)
    
        # Run the model with the appropriate simulation method
        # There are five methods in gillespy2 package:
        #    1) ODE (using ordinary differential equations, it is a deterministic approach)
        #    2) SSA (using stochastic simulation algorithm), it is a stochastic approach
        #    3) Tau-leaping algorithm, another stochastic approach
        #    4) CLE (chemical langevin equation), another stochastic approach
        #    5) Tau-Hybrid method: a combination between stochastic and deterministic approach 
        result = model.run(algorithm = self.method) 
            
            
        return result
        
        
methods = ['ODE', 'SSA', 'CLE', 'Tau-Leaping', 'Tau-Hybrid']  
num_steps = [1000, 1000, 500, 250, 1000]

results = []
for i in range(len(methods)):
    model = GillespySimulator(0.1, 0.05, 100, 0, 0, 100, num_steps[i], methods[i])
    result = model.simulation()
    results.append(result)
            
    
    
# Visualize the simulations
plt.figure(figsize=(17,10))

# ODE approach
sp1_population = results[0]['reactant']
sp2_population = results[0]['product']
simulation_time = results[0]['time']   
plt.subplot(2,3,1)
plt.plot(simulation_time, sp1_population, label = 'Sp1', linewidth = 0.2)
plt.plot(simulation_time, sp2_population, label = 'Sp2', linewidth = 0.2)
plt.title('Deterministic Simulation of a System using ODE', fontsize = 10)
    
# SSA approach
sp1_population1 = results[1]['reactant']
sp2_population1 = results[1]['product']
simulation_time1 = results[1]['time'] 
plt.subplot(2,3,2)
plt.plot(simulation_time1, sp1_population1, label = 'Sp1', linewidth = 0.2)
plt.plot(simulation_time1, sp2_population1, label = 'Sp2', linewidth = 0.2)
plt.title('Stochastic Simulation of a System using SSA', fontsize = 10)
    
# CLE approach
sp1_population2 = results[2]['reactant']
sp2_population2 = results[2]['product']
simulation_time2 = results[2]['time'] 
plt.subplot(2,3,4)
plt.scatter(simulation_time2, sp1_population2, s = 1)
plt.scatter(simulation_time2, sp2_population2, s = 1)
plt.plot(simulation_time2, sp1_population2, label = 'Sp1', linewidth = 0.2)
plt.plot(simulation_time2, sp2_population2, label = 'Sp2', linewidth = 0.2)
plt.title('Stochastic Simulation of a System using CLE', fontsize = 10)
    
# Tau-Leaping approach
sp1_population3 = results[3]['reactant']
sp2_population3 = results[3]['product']
simulation_time3 = results[3]['time'] 
plt.subplot(2,3,3)
plt.scatter(simulation_time3, sp1_population3, s = 1)
plt.scatter(simulation_time3, sp2_population3, s = 1)
plt.plot(simulation_time3, sp1_population3, linestyle = 'dashed', label = 'Sp1', linewidth = 0.4)
plt.plot(simulation_time3, sp2_population3, linestyle = 'dashed', label = 'Sp2', linewidth = 0.4)
plt.title('Stochastic Simulation of a System using t-leaping', fontsize = 10)
    
# Tau-Hybrid approach
sp1_population4 = results[4]['reactant']
sp2_population4 = results[4]['product']
simulation_time4 = results[4]['time'] 
plt.subplot(2,3,(5,6))
plt.plot(simulation_time4, sp1_population4, linewidth = 0.2)
plt.plot(simulation_time4, sp2_population4, linewidth = 0.2)
plt.title('Stochastic Simulation of a System using t-hybrid', fontsize = 10)

plt.legend()

plt.plot()
    

