import numpy as np
import time


class EulerSimulator(object):

    """ Simulation using Euler method """

    def __init__(
        self,
        model=None,
        start=0,
        stop=10,
        epochs=1000,
        seed=42,
        model_name="Euler Method",
        **kwargs
    ):

        self.model = model
        self.start = start
        self.stop = stop
        self.epochs = epochs
        self.seed = seed
        self.model_name = model_name

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)
        else:
            raise ValueError("Before simulating a model, please ensure that you have instantiated the biostoch.model.Model() object.")

        self.species = {}
        self.parameters = {}
        self.time = {}

        """
        Args:
            model: a class created by "biostoch.model.Model",
                   contains all necessary information used in simulation with SSA.
            start: an integer or a float that defines the start time of the simulation.
            stop: an integer or a float that defines the stop time of the simulation.
            epochs: an integer defines the number of iterations.
            seed: an integer parameter used to initialize the random number generator.
            model_name: the name of the simulation method "Stochastic Simulation Algorithm".
            **kwargs: a special parameter that allows passing additional keyword arguments to the function.

            species: an empty dictionary in which the calculated concentrations of the species are stored.
            parameters: an empty dictionary in which the rate constants of the model are stored.
            time: an empty dictionary in which the simulation duration is stored. 
        """

    def reset(self):

        """ Resets the model species and parameters dictionaries"""

        self.species = {}
        self.parameters = {}
        self.time = {}

    def initialize_parameters(self, model, start, stop, epochs):

        """
        Initializes the model species dictionary

        Args:
            model: a class created by "biostoch.model.Model"
                   contains all necessary information used in simulation with SSA.
            start: an integer or a float that defines the start time of the simulation.
            stop: an integer or a float that defines the stop time of the simulation.
            epochs: an integer defines the number of iterations.
        Returns:
            species: a dictionary contains initialized concentration of each
                     species and also initialized simulation time in the system.
                     each dictionary's key is the name of one species
                     and each value the corresponding initialized concentration.
            parameters: a dictionary contains the rate constant of each reaction each key
                        correspond to the name of the rate constant and each value is its value.
        """

        species = {}
        parameters = {}
        species["Time"] = np.linspace(start, stop, epochs)

        for specie in model.components:
            species[specie] = np.zeros(epochs)
            species[specie][0] = getattr(model, specie)
        for parameter in model.params:
            parameters[parameter] = getattr(model, parameter)

        return species, parameters

    def compute_rates(self, species, model, step):

        """

        Args:
            species: a dictionary in which the calculated concentrations of the species are stored.
            model: a class created by "biostoch.model.Model".
            step: integer that represents the current iteration number or time step in the simulation process.

        Returns:
            rates: a dictionary containing the calculated rates of each reaction at the current time step.

        """

        rates = {}
        for specie in species.keys():
            if specie != "Time":
                rate = ""
                split_rate = model.ROC_[specie].split()
                for component in split_rate:
                    if component in self.model.params.keys():
                        rate += " " + str(self.model.params[component])
                    elif component in model.signs:
                        rate += " " + component
                    elif component in self.model.components:
                        rate += " " + str(species[component][step - 1])
                    else:
                        ValueError(f"This component: {component} is not a valid component!")
                rates[specie] = rate

        return rates

    def simulate(self):

        """Runs the simulation"""

        start_simulation = time.time()

        species, parameters = self.initialize_parameters(
            model=self.model,
            start=self.start,
            stop=self.stop,
            epochs=self.epochs
        )

        tau = species["Time"][3] - species["Time"][2]

        for i in range(1, self.epochs):

            rates = self.compute_rates(
                species=species,
                model=self.model,
                step=i
            )

            for specie, concentration in species.items():
                if specie != "Time":
                    if specie in rates.keys():
                        species[specie][i] = species[specie][i - 1] + (eval(rates[specie]) * tau)
                    else:
                        raise ValueError(f"The rate equation for '{specie}' is not defined!")

        self.species = species
        self.parameters = parameters
        stop_simulation = time.time()
        self.time["Simulation Duration"] = stop_simulation - start_simulation


class RungeKuttaSimulator(object):

    """ Simulation using Runge Kutta method """

    def __init__(
        self,
        model=None,
        start=0,
        stop=10,
        epochs=1000,
        seed=42,
        model_name="Runge-Kutta Algorithm",
        **kwargs
    ):

        self.model = model
        self.start = start
        self.stop = stop
        self.epochs = epochs
        self.seed = seed
        self.model_name = model_name

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)
        else:
            raise ValueError("Before simulating a model, please ensure that you have instantiated the biostoch.model.Model() object.")

        self.species = {}
        self.parameters = {}
        self.time = {}

        """
        Args:
            model: a class created by "biostoch.model.Model",
                   contains all necessary information used in simulation with SSA.
            start: an integer or a float that defines the start time of the simulation.
            stop: an integer or a float that defines the stop time of the simulation.
            epochs: an integer defines the number of iterations.
            seed: an integer parameter used to initialize the random number generator.
            model_name: the name of the simulation method "Stochastic Simulation Algorithm".
            **kwargs: a special parameter that allows passing additional keyword arguments to the function.

            species: an empty dictionary in which the calculated concentrations of the species are stored.
            parameters: an empty dictionary in which the rate constants of the model are stored.
            time: an empty dictionary in which the simulation duration is stored. 
        """

    def reset(self):

        """ Resets the model species and parameters dictionaries"""

        self.species = {}
        self.parameters = {}
        self.time = {}

    def initialize_parameters(self, model, start, stop, epochs):

        """
        Initializes the model species dictionary

            Args:
                model: a class created by "biostoch.model.Model"
                       contains all necessary information used in simulation with SSA.
                start: an integer or a float that defines the start time of the simulation.
                stop: an integer or a float that defines the stop time of the simulation.
                epochs: an integer defines the number of iterations.
            Returns:
                species: a dictionary contains initialized concentration of each
                         species and also initialized simulation time in the system.
                         each dictionary's key is the name of one species
                         and each value the corresponding initialized concentration.
                parameters: a dictionary contains the rate constant of each reaction each key
                            correspond to the name of the rate constant and each value is its value.
        """

        species = {}
        parameters = {}
        species["Time"] = np.linspace(start, stop, epochs)

        for specie in model.components:
            species[specie] = np.zeros(epochs)
            species[specie][0] = getattr(model, specie)
        for parameter in model.params:
            parameters[parameter] = getattr(model, parameter)

        return species, parameters

    def compute_rates(self, species, model, step):

        """

        Args:
            species: a dictionary in which the calculated concentrations of the species are stored.
            model: a class created by "biostoch.model.Model".
            step: integer that represents the current iteration number or time step in the simulation process.

        Returns:
            rates: a dictionary containing the calculated rates of each reaction at the current time step.

        """

        rates = {}
        for specie in species.keys():
            if specie != "Time":
                rate = ""
                split_rate = model.ROC_[specie].split()
                for component in split_rate:
                    if component in self.model.params.keys():
                        rate += " " + str(self.model.params[component])
                    elif component in model.signs:
                        rate += " " + component
                    elif component in self.model.components:
                        rate += " " + str(species[component][step - 1])
                    else:
                        ValueError(f"This component: {component} is not a valid component!")
                rates[specie] = rate

        return rates

    def simulate(self):

        """Runs the simulation"""

        start_simulation = time.time()

        species, parameters = self.initialize_parameters(
            model=self.model,
            start=self.start,
            stop=self.stop,
            epochs=self.epochs
        )

        tau = species["Time"][3] - species["Time"][2]

        for i in range(1, self.epochs):
            rates = self.compute_rates(
                species=species,
                model=self.model,
                step=i
            )

            k1 = {}
            k2 = {}
            k3 = {}
            k4 = {}

            for specie, concentration in species.items():
                if specie != "Time":
                    k1[specie] = eval(rates[specie]) * tau
                    k2[specie] = eval(rates[specie]) * tau
                    k3[specie] = eval(rates[specie]) * tau
                    k4[specie] = eval(rates[specie]) * tau

            for specie, concentration in species.items():
                if specie != "Time":
                    species[specie][i] = species[specie][i - 1] + (1 / 6) * (k1[specie] + 2 * k2[specie] + 2 * k3[specie] + k4[specie])

        self.species = species
        self.parameters = parameters
        stop_simulation = time.time()
        self.time["Simulation Duration"] = stop_simulation - start_simulation


