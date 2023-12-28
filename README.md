
### BioStoch

BioStoch is a library that contains various stochastic and deterministic simulation methods used in computational biology.

The library is under development and cannot yet be used.

How to install BioStoch:
```bash
    pip  install numpy matplotlib seaborn pandas  # install required libraries
    pip install BioStoch # install BioStoch library
```

How to use BioStoch:
```python
    from BioStoch import SSA, TauLeaping, CLE, RRE, Visualization



    model1 = SSA.StochSim()  # This model will simulate a biological system with The Stochastic Simulation Algorithm (SSA)
    
    model2 = TauLeaping.TauLeap()  # This model will simulate a biological system with The Tau-Leaping Algorithm
    
    model3 = RRE.RateEquation()  # This model will simulate a biological system with The Reaction Rate Equations (RRE) algorithm
    
    model5 = CLE.ChemLong()  # This model will simulate a biological system with The Chemical Langevin Equation (CLE) algorithm

    


    plot = Visualization.Plot()  # The Visualization module will contain different methods to visualize the simulation outputs using matplotlib & seaborn
```

[GitHub Homepage](https://github.com/LoqmanSamani/biostoch)


