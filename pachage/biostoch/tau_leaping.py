import numpy as np
import time
import model



class TauLeaping(object):
    """ Simulation using Tau-Leaping method """

    def __init__(self, model=None, start=0.0, stop=10.0, max_epochs=100, seed=42,
                 steady_state=None, epsilon=0.03, call_tau=None, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.steady_state = steady_state
        self.epsilon = epsilon
        self.tau = (self.stop-self.start) / self.max_epochs
        self.call_tau = call_tau

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)
        else:
            raise ValueError("Before simulating a model, please ensure that you have instantiated the biostoch.model.Model() object.")

        self.model_name = "Tau-Leaping Algorithm"
        self.species = None
        self.parameters = None
        self.time = {}

    def reset(self):
        self.species = None
        self.parameters = None
        self.time = {}


    def initialize_parameters(self, model, start, max_epochs):

        species = {}
        parameters = {}

        species["Time"] = np.zeros(max_epochs)
        species["Time"][0] = start
        for specie in model.components:
            species[specie] = np.zeros(max_epochs)
            species[specie][0] = getattr(model, specie)
        for parameter in self.model.params:
            parameters[parameter] = getattr(model, parameter)

        return species, parameters


    def compute_propensity_sum(self, species, parameters, propensities, step):

        propensity_sum = 0.0
        propensities_ = {}
        last_step = {}
        for specie, concentration in species.items():
            if specie != "Time":
                last_step[specie] = concentration[step-1]
        for parameter, value in parameters.items():
            last_step[parameter] = value

        for reaction, propensity in propensities.items():
            propensity_ = eval(propensity, last_step)
            propensity_sum += propensity_
            propensities_[reaction] = propensity_

        return propensity_sum, propensities_


    def compute_tau(self, species, model, step, epsilon):

        X = np.array([species[con][step - 1] for con in species.keys() if con != "Time"])
        v = []

        for key, val in model.coeffs_.items():
            d = []
            for sp in model.react_sps[key]:
                d.append(val[sp])
            v.append(d)

        v = np.array(v)
        R = []

        comp = model.params
        X1 = {key: val[step - 1] for key, val in species.items() if key != "Time"}
        comp.update(X1)

        s = 0
        for react, rate in model.rates_.items():
            if react == model.react_names[s]:
                R.append(eval(rate, comp))
                s += 1
        R = np.array(R)

        mu_values = []
        sigma_squared_values = []

        for i in range(len(X)):
            mu_i = np.sum(v[i] * R)
            sigma_squared_i = np.sum((v[i] ** 2) * R)
            mu_values.append(mu_i)
            sigma_squared_values.append(sigma_squared_i)

        g_values = [np.argmax(v[i]) + 1 for i in range(len(X))]

        tau_values = []
        for i in range(len(X)):
            tau_i = min(max(epsilon * X[i] / g_values[i], 1) / abs(mu_values[i]),
                        max(epsilon * X[i] / g_values[i], 1) ** 2 / sigma_squared_values[i])
            tau_values.append(tau_i)

        return min(tau_values)


    def compute_lambdas(self, species, parameters, propensities, tau, step):

        last_step = {}

        for specie, concentration in species.items():
            if specie != "Time":
                last_step[specie] = concentration[step - 1]

        for parameter, value in parameters.items():
            last_step[parameter] = value

        lambdas = {}

        for reaction, propensity in propensities.items():
            lambda_value = eval(propensity, last_step) * tau

            if lambda_value < 0.0:
                lambdas[reaction] = 0
            else:
                lambdas[reaction] = lambda_value

        return lambdas


    def num_reaction(self, lambdas):
        num_reaction_ = {}
        for reaction, lambda_ in lambdas.items():
            num_reaction_[reaction] = np.random.poisson(lambda_)

        return num_reaction_


    def update(self, species, model, num_reaction, step, tau):

        species["Time"][step] = species["Time"][step - 1] + tau

        for reaction, formula in model.reacts_.items():
            split_formula = formula.split()
            index = [index for index, value in enumerate(split_formula) if value == '->']

            if len(index) != 1:
                print(f"Error: Each reaction should have exactly one '->', but there are {len(index)} in {reaction}.")

        component_reaction = {}
        for component in model.components:
            num_reaction_ = sum([num_reaction[reaction_] * model.coeffs_[reaction_][component] for reaction_ in model.react_names if component in model.react_sps[reaction_]])
            component_reaction[component] = num_reaction_

        for component, value in component_reaction.items():
            species[component][step] = species[component][step - 1] + value

        return species


    def resize_species(self, species, step):

        if step >= len(species["Time"]):

            new_max_steps = len(species["Time"]) * 2

            for specie, concentration in species.items():
                pad_width = (0, new_max_steps - len(concentration))
                species[specie] = np.pad(concentration, pad_width, mode='constant')

        return species


    def final_resize_species(self, species, final_step):
        for specie in species.keys():
            species[specie] = species[specie][:final_step]
        return species


    def simulate(self):

        start_simulation = time.time()

        species, parameters = self.initialize_parameters(model=self.model, start=self.start, max_epochs=self.max_epochs)

        step = 1
        while step < self.max_epochs:

            propensity_sum, propensities_ = self.compute_propensity_sum(
                species=species,
                parameters=parameters,
                propensities=self.model.rates_,
                step=step
            )

            if propensity_sum == 0 and self.steady_state:
                print(f"Simulation reached steady state (iteration: {step}). No further changes are occurring.")
                break

            if self.call_tau:
                tau = self.compute_tau(species=species, model=self.model, step=step, epsilon=self.epsilon)
            else:
                tau = self.tau

            lambdas = self.compute_lambdas(
                species=species,
                parameters=self.model.params,
                propensities=self.model.rates_,
                tau=tau,
                step=step
            )

            num_reaction = self.num_reaction(
                lambdas=lambdas
            )

            species = self.update(
                species=species,
                model=self.model,
                num_reaction=num_reaction,
                step=step,
                tau=tau
            )

            step += 1

            species = self.resize_species(
                species=species,
                step=step
            )

        species = self.final_resize_species(
            species=species,
            final_step=step
        )

        self.species = species
        self.parameters = parameters
        stop_simulation = time.time()
        self.time["Simulation Duration"] = stop_simulation - start_simulation





