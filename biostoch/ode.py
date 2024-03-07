import numpy as np
import model
import time


class EulerSimulator(object):
    def __init__(self, model=None, start=0, stop=10, epochs=1000, seed=42, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.epochs = epochs
        self.seed = seed

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)
        else:
            raise ValueError("Before simulating a model, please ensure that you have instantiated the biostoch.model.Model() object.")

        self.model_name = "Euler Method"
        self.species = None
        self.parameters = None
        self.time = {}

    def reset(self):
        self.species = None
        self.parameters = None
        self.time = {}

    def initialize_paramsters(self, model, start, stop, epochs):

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
        start_simulation = time.time()

        species, parameters = self.initialize_paramsters(
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
    def __init__(self, model=None, start=0, stop=10, epochs=1000, seed=42, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.epochs = epochs
        self.seed = seed

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)
        else:
            raise ValueError("Before simulating a model, please ensure that you have instantiated the biostoch.model.Model() object.")

        self.model_name = "Runge-Kutta Algorithm"
        self.species = None
        self.parameters = None
        self.time = {}

    def reset(self):
        self.species = None
        self.parameters = None
        self.time = {}

    def initialize_paramsters(self, model, start, stop, epochs):

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
        start_simulation = time.time()

        species, parameters = self.initialize_paramsters(
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



m = model.Model()
m.parameters({"K1": 0.1, "K2": 0.05})

m.species({"A": 100.0, "B": 0.0}, {"A": "K2 * B - K1 * A", "B": "K1 * A - K2 * B"})

m.reactions({"reaction1": "A -> B", "reaction2": "B -> A"},
           {"reaction1": "K1 * A", "reaction2": "K2 * B"})


model1 = EulerSimulator(model=m, start=0, stop=100, epochs=1000)

print(model1.model_name)