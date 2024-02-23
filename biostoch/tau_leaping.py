import numpy as np
import time
import model
import matplotlib.pyplot as plt


class TauLeaping(object):
    def __init__(self, model=None, start=0.0, stop=100.0, max_epochs=100, seed=42, steady_state=None, epsilon=0.03, call_tau=None, **kwargs):

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

        self.species = None
        self.parameters = None
        self.lams = []
        self.num_r = []

    def param_init(self, model, start, max_epochs):

        species = {}
        parameters = {}
        epochs = max_epochs

        species["Time"] = np.zeros(epochs)
        species["Time"][0] = start
        for specie in model.components:
            species[specie] = np.zeros(epochs)
            species[specie][0] = getattr(model, specie)
        for parameter in self.model.params:
            parameters[parameter] = getattr(model, parameter)

        return species, parameters

    def propensity_sum(self, step, propensities, species, params):

        propensity_sum = 0.0
        props = {}
        last_step = {}
        for key, val in species.items():
            if key != "Time":
                last_step[key] = val[step-1]
        for key, val in params.items():
            last_step[key] = val

        for reaction, propensity in propensities.items():
            propensity = eval(propensity, last_step)
            propensity_sum += propensity
            props[reaction] = propensity

        return propensity_sum, props

    def calculate_tau(self, species, model, step, epsilon):
        """
        Calculate the time step tau based on the given parameters.

        Args:
        epsilon (float): A small positive constant.
        X (ndarray): Array containing state variables X_i.
        v (ndarray): Array containing auxiliary values v_ij.
        R (ndarray): Array containing reaction rates R_j.

        Returns:
        float: The calculated time step tau.
        """
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

        # Initialize lists to store μ_i and σ_i^2 values
        mu_values = []
        sigma_squared_values = []

        # Calculate μ_i and σ_i^2 for each state variable X_i
        for i in range(len(X)):
            mu_i = np.sum(v[i] * R)
            sigma_squared_i = np.sum((v[i] ** 2) * R)
            mu_values.append(mu_i)
            sigma_squared_values.append(sigma_squared_i)

        # Determine the highest order event g_i for each state variable
        g_values = [np.argmax(v[i]) + 1 for i in range(len(X))]

        # Calculate the time step τ based on the formula
        tau_values = []
        for i in range(len(X)):
            tau_i = min(max(epsilon * X[i] / g_values[i], 1) / abs(mu_values[i]),
                        max(epsilon * X[i] / g_values[i], 1) ** 2 / sigma_squared_values[i])
            tau_values.append(tau_i)

        # Return the minimum value of τ
        return min(tau_values)

    def calculate_lambda(self, propensities, species, params, tau, step):
        last_step = {}
        for key, val in species.items():
            if key != "Time":
                last_step[key] = val[step - 1]
        for key, val in params.items():
            last_step[key] = val
        lams = []
        lambdas = {}
        for react, prop in propensities.items():
            lambda_val = eval(prop, last_step) * tau
            lams.append(lambda_val)
            if lambda_val <= 0.0:
                lambdas[react] = 1
            else:
                lambdas[react] = lambda_val

        return lambdas, lams

    def num_reacts(self, lambdas):
        n_reacts = {}
        for react, lam in lambdas.items():
            n_reacts[react] = np.random.poisson(lam)

        return n_reacts

    def update(self, species, model, num_reacts, step, tau):

        species["Time"][step] = species["Time"][step - 1] + tau

        for react, prop in model.reacts_.items():
            split_prop = prop.split()
            split_index = [index for index, value in enumerate(split_prop) if value == '->']

            if len(split_index) != 1:
                print(
                    f"Error: Each reaction should have exactly one '->', but there are {len(split_index)} in {react}.")
                continue

        comp_react = {}
        for comp in model.components:
            r = sum([num_reacts[e] * model.coeffs_[e][comp] for e in model.react_names if comp in model.react_sps[e]])
            comp_react[comp] = r

        for comp, val in comp_react.items():
            species[comp][step] = species[comp][step-1] + val

        return species

    def change_size(self, species, step):

        if step >= len(species["Time"]):

            new_max_steps = len(species["Time"]) * 2

            for key, val in species.items():
                pad_width = (0, new_max_steps - len(val))
                species[key] = np.pad(val, pad_width, mode='constant')

        return species

    def resize_species(self, species, final_step):
        for key in species.keys():
            species[key] = species[key][:final_step]
        return species

    def simulate(self):

        species, parameters = self.param_init(
            model=self.model,
            start=self.start,
            max_epochs=self.max_epochs,
        )

        step = 1
        while step < self.max_epochs:

            a_sum, props = self.propensity_sum(
                step=step,
                propensities=self.model.rates_,
                species=species,
                params=parameters
            )

            if a_sum == 0 and self.steady_state:
                print(f"Simulation reached steady state (iteration: {step}). No further changes are occurring.")
                break

            if self.call_tau:
                tau = self.calculate_tau(
                    species=species,
                    model=self.model,
                    step=step,
                    epsilon=self.epsilon
                )
            else:
                tau = self.tau

            lambdas, lams = self.calculate_lambda(
                propensities=self.model.rates_,
                species=species,
                params=self.model.params,
                tau=tau,
                step=step
            )
            self.lams.append(lams)
            num_reacts = self.num_reacts(
                lambdas=lambdas
            )
            self.num_r.append(num_reacts)

            species = self.update(
                species=species,
                model=self.model,
                num_reacts=num_reacts,
                step=step,
                tau=tau
            )

            step += 1

            species = self.change_size(
                species=species,
                step=step
            )
            species = species

        species = self.resize_species(
            species=species,
            final_step=step
        )

        self.species = species
        self.parameters = parameters



