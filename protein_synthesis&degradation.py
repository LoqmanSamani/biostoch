#!/usr/bin/env python
# coding: utf-8

# In[48]:


import numpy as np
import matplotlib.pyplot as plt

k0 = 0.01
k1 = 1
k2 = 5
mRNA = 3
a_as = 3
pro_conc = 0



t_start = 0
t_end = 10
n_points = 100

Time = np.linspace(t_start, t_end, n_points)
mRNA_conc = np.linspace(mRNA, 0, n_points)
steady_state_sol = (k0+k1*mRNA_conc) / k2
steady_state_sol
mRNA_conc

    
plt.figure(figsize=(10, 7))    
plt.plot(mRNA_conc, steady_state_sol)
plt.xlabel('Signal Strength (S)', fontsize=12)
plt.ylabel('The Steady-Rtate Response (R_ss)', fontsize=12)
plt.title('The Signal-Response Curve of Protein Synthesis & Degradation', fontsize=16)
plt.grid(True)
plt.show()





# In[47]:


import numpy as np
import matplotlib.pyplot as plt

k0 = 0.01
k1 = 1
k2 = 5
mRNA = 5
a_as = 3
pro_conc = 0



t_start = 0
t_end = 10
n_points = 100

Time = np.linspace(t_start, t_end, n_points)
mRNA_conc = np.linspace(mRNA, 0, n_points)
steady_state_sol = (k0+k1*mRNA_conc) / k2
steady_state_sol
mRNA_conc

Rate = k0+(k1*mRNA_conc) - (k2* pro_conc)
plt.figure(figsize=(10, 7))
plt.plot(steady_state_sol, Rate, label='Protein Responce')
points_y = [1, 2, 3]
for point_y in points_y:
    plt.scatter(steady_state_sol[np.abs(Rate - point_y).argmin()], point_y, color='red', marker='o')
    plt.axhline(y=point_y, color='red', linestyle='--', label= 'mRNA Concentrations')

plt.xlim(0, 1)  
plt.ylim(0, 5) 
plt.xlabel('Response (Concentration of Protein)',fontsize=12)
plt.ylabel('Rate Equation(dR/dt)',fontsize=12)
plt.title('The Rate Curve of Protein Synthesis & Degradation',fontsize=16)
plt.legend()
plt.grid(True)
plt.show()

