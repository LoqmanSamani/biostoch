import numpy as np
import matplotlib.pyplot as plt
import seaborn

import model
import ode
import ssa
import tau_leaping
import cle

class Visualization(object):
    def __init__(self, model=None, model_name="Simulation Result", *args):

        self.model = model
        self.model_name = model_name

    def extract_species(self, model):

        simulation_result = None

        if model.species and isinstance(model.species, dict):
            simulation_result = model.species
            simulation_result["Model Name"] = model.model_name
        elif isinstance(model, dict):
            simulation_result = model
            simulation_result["Model Name"] = self.model_name

        else:
            raise TypeError("Please provide the simulation result either as a model object (simulated with biostoch) or in dictionary format.")

        return simulation_result

    def check_length_species(self, simulation_result):

        check = False
        result_length = 0
        step = 0
        for specie, result in simulation_result.items():
            if specie != "Model Name":
                if step == 0:
                    result_length = len(result)
                if len(result) != result_length:
                    check = True

        if check:
            raise ValueError("All species concentrations must have the same length.")

    def plot(self, model=None, x=None, y=None, plot_size=None, num_species=1):

        biostoch_model = False
        another_model = False
        if self.model:
            simulation_result = self.extract_species(model=self.model)
            self.check_length_species(simulation_result=simulation_result)
            biostoch_model = True

        elif model:
            simulation_result = self.extract_species(model=self.model)
            self.check_length_species(simulation_result=simulation_result)
            another_model = True

        else:
            raise ValueError("PLease provide a model!")

        if biostoch_model:
            if simulation_result["Model Name"] == "Runge-Kutta Algorithm" or simulation_result["Model Name"] == "Euler Method":
                if plot_size:
                    plt.figure(figsize=plot_size)
                else:
                    plt.figure(figsize=(10, 8))
                for specie, concentration in simulation_result.items():
                    if specie != "Time" and specie == "Model Name":
                        plt.plot(simulation_result["Time"], concentration, label=specie)
                        plt.xlabel("Time")
                        plt.ylabel("Concentration")
                        plt.title(f"Simulation with {simulation_result['Model Name']}")

                plt.legend()
                plt.show()

            elif simulation_result["Model Name"] == "Stochastic Simulation Algorithm":
                if plot_size:
                    plt.figure(figsize=plot_size)
                else:
                    plt.figure(figsize=(10, 8))
                for specie, concentration in simulation_result.items():
                    if specie != "Time" and specie == "Model Name":
                        plt.plot(simulation_result["Time"], concentration, label=specie)
                        plt.xlabel("Time")
                        plt.ylabel("Concentration")
                        plt.title(f"Simulation with {simulation_result['Model Name']}")

                plt.legend()
                plt.show()

            elif simulation_result["Model Name"] == "Tau-Leaping Algorithm":

                if plot_size:
                    plt.figure(figsize=plot_size)
                else:
                    plt.figure(figsize=(10, 8))

                for specie, concentration in simulation_result.items():
                    if specie != "Time" and specie == "Model Name":

                        plt.plot(
                            simulation_result["Time"],
                            concentration,
                            label=specie,
                            marker='o',
                            linestyle='dashed',
                            markersize=6
                        )

                        plt.xlabel("Time")
                        plt.ylabel("Concentration")
                        plt.title(f"Simulation with {simulation_result['Model Name']}")

                plt.legend()
                plt.show()

            elif simulation_result["Model Name"] == "Chemical Langevin Equation":

                if plot_size:
                    plt.figure(figsize=plot_size)
                else:
                    plt.figure(figsize=(10, 8))

                for specie, concentration in simulation_result.items():

                    if specie != "Time" and specie == "Model Name":

                        plt.plot(
                            simulation_result["Time"],
                            concentration,
                            label=specie,
                            marker='o',
                            linestyle='dashed',
                            markersize=6
                        )

                        plt.xlabel("Time")
                        plt.ylabel("Concentration")
                        plt.title(f"Simulation with {simulation_result['Model Name']}")

                plt.legend()
                plt.show()

            if another_model:

                if simulation_result["Model Name"] == self.model_name:

                    if plot_size:
                        plt.figure(figsize=plot_size)
                    else:
                        plt.figure(figsize=(10, 8))

                    for specie, concentration in simulation_result.items():

                        if "Time" in simulation_result.keys():

                            if specie != "Time" and specie == "Model Name":

                                plt.plot(simulation_result["Time"], concentration, label=specie)
                                plt.xlabel("Time")
                                plt.ylabel("Concentration")
                                plt.title(f"Simulation with {simulation_result['Model Name']}")

                            plt.legend()
                            plt.show()

                        else:
                            if num_species == 1:

                                if plot_size:
                                    plt.figure(figsize=plot_size)
                                else:
                                    plt.figure(figsize=(10, 8))

                                plt.plot(x, y)
                                plt.show()

                            else:
                                for i in range(num_species):

                                    plt.plot(x, y[i])

                                plt.show()






m = model.Model()
m.parameters({"K1": 0.1, "K2": 0.05})

m.species({"A": 100.0, "B": 0.0}, {"A": "K2 * B - K1 * A", "B": "K1 * A - K2 * B"})

m.reactions({"reaction1": "A -> B", "reaction2": "B -> A"},
           {"reaction1": "K1 * A", "reaction2": "K2 * B"})


model1 = ode.EulerSimulator(model=m, start=0, stop=100, epochs=1000)
model2 = ode.RungeKuttaSimulator(model=m, start=0, stop=100, epochs=1000)
model3 = ssa.GillespieSimulator(model=m, start=0, stop=100, max_epochs=1000)
model4 = tau_leaping.TauLeaping(model=m, start=0, stop=100, max_epochs=100)
model5 = cle.ChemicalLangevin(model=m, max_epochs=100)

model1.simulate()
model2.simulate()
model3.simulate()
model4.simulate()

model11 = Visualization(model1, model_name="Euler_Method")
model12 = Visualization(model2, model_name="Runge_Kutta_Algorithm")
model13 = Visualization(model3, model_name="Stochastic_Simulation_Algorithm")
model14 = Visualization(model4, model_name="Tau_Leaping_Algorithm")
model15 = Visualization(model5, model_name="Chemical_Langevin_Equation")

model11.plot()
model12.plot()
model13.plot()
model14.plot()
model15.plot()

















