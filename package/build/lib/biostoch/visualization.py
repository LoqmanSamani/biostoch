import matplotlib.pyplot as plt
import numpy as np


class Visualization(object):
    def __init__(
            self,
            model=None,
            model_name="Simulation Result",
            **kwargs
    ):

        self.model = model
        self.model_name = model_name

    def extract_species(self, model):

        simulation_result = None

        if model.species and isinstance(model.species, dict):
            simulation_result = model.species
            simulation_result["Model Name"] = str(model.model_name)
        elif isinstance(model, dict):
            simulation_result = model
            simulation_result["Model Name"] = str(self.model_name)

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
                title = simulation_result["Model Name"]
                if plot_size:
                    plt.figure(figsize=plot_size)
                else:
                    plt.figure(figsize=(10, 8))
                for specie, concentration in simulation_result.items():

                    if specie != "Time" and specie != "Model Name":
                        plt.plot(simulation_result["Time"], concentration, label=specie)
                        plt.xlabel("Time")
                        plt.ylabel("Concentration")
                        plt.title(f"Simulation with {title}")

                plt.legend()
                plt.show()

            elif simulation_result["Model Name"] == "Stochastic Simulation Algorithm":
                if plot_size:
                    plt.figure(figsize=plot_size)
                else:
                    plt.figure(figsize=(10, 8))
                for specie, concentration in simulation_result.items():
                    if specie != "Time" and specie != "Model Name":
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
                    if specie != "Time" and specie != "Model Name":

                        plt.plot(
                            simulation_result["Time"],
                            concentration,
                            label=specie,
                            marker='o',
                            linestyle='dashed',
                            markersize=4
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

                    if specie != "Time" and specie != "Model Name":

                        plt.plot(
                            simulation_result["Time"],
                            concentration,
                            label=specie,
                            marker='o',
                            linestyle='dashed',
                            markersize=2
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



