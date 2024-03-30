import numpy as np
import time


class TauLeaping(object):

    """ Simulation using Tau-Leaping method """

    def __init__(
        self,
        model=None,
        start=0.0,
        stop=10.0,
        max_epochs=100,
        seed=42,
        steady_state=False,
        epsilon=0.03,
        call_tau=False,
        model_name="Tau-Leaping Algorithm",
        **kwargs
    ):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.steady_state = steady_state
        self.epsilon = epsilon
        self.tau = (self.stop-self.start) / self.max_epochs
        self.call_tau = call_tau
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
            max_epochs: an integer defines the maximum number of iterations.
            seed: an integer parameter used to initialize the random number generator.
            steady_state: Boolean value (True or False); if true, 
                          the simulation is stopped as soon as the model has reached the steady state.
            epsilon: a float value that is less than one and is used as a fixed tolerance for the calculation of tau.
            tau: a float or an integer that represents the time step size.
            cal_tau: Boolean value (True or False); if true, tau is calculated in each step.
            model_name: the name of the simulation method "Tau-Leaping Algorithm".
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

    def initialize_parameters(self, model, start):

        """
        Initializes the model species dictionary

        Args:
            model: a class created by "biostoch.model.Model"
                   contains all necessary information used in simulation with Tau Leaping Algorithm.
            start: an integer or a float that defines the start time of the simulation.
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

        species["Time"] = [start]

        for specie in model.components:
            species[specie] = [getattr(model, specie)]

        for parameter in self.model.params:
            parameters[parameter] = getattr(model, parameter)

        return species, parameters

    def compute_propensity_sum(self, propensities, species, parameters):

        """
        Computes sum of the propensities
            Args:
                propensities: a dictionary contains propensity functions of the reactions.
                species: a dictionary in which the calculated concentrations of the species are stored.
                parameters: a dictionary in which the rate constants of the model are stored.
            Returns:
                propensity_sum: a float value, sum of the propensities.
                propensities_: a dictionary contains the propensity values of the reactions.

        """

        propensity_sum = 0.0
        propensities_ = {}
        last_step = {}

        for specie, concentration in species.items():
            if specie != "Time":
                last_step[specie] = concentration[-1]

        for parameter, value in parameters.items():
            last_step[parameter] = value

        for reaction, propensity in propensities.items():
            propensity_ = eval(propensity, last_step)
            propensity_sum += propensity_
            propensities_[reaction] = propensity_

        return propensity_sum, propensities_

    def compute_tau(self, species, model, epsilon):

        """

        Args:
            species: a dictionary in which the calculated concentrations of the species are stored.
            model: a class created by "biostoch.model.Model"
            epsilon: a float value that is less than one and is used as a fixed tolerance for the calculation of tau.

        Returns:
            tau: a float or an integer, calculated tau.
        """

        X = np.array([species[con][-1] for con in species.keys() if con != "Time"])
        v = []

        for key, val in model.coeffs_.items():
            d = []
            for sp in model.react_sps[key]:
                d.append(val[sp])
            v.append(d)

        v = np.array(v)
        R = []

        comp = model.params
        X1 = {key: val[-1] for key, val in species.items() if key != "Time"}
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

    def compute_lambdas(self, species, parameters, propensities, tau):

        """

        Args:
            species: a dictionary in which the calculated concentrations of the species are stored.
            parameters: a dictionary in which the rate constants of the model are stored.
            propensities: a dictionary contains the propensity values of the reactions.
            tau: a float or an integer that represents the time step size.

        Returns:
            lambda: calculated poisson distribution parameter (lambda),
                    (the mean number of events within a given interval of time or space)

        """

        last_step = {}

        for specie, concentration in species.items():
            if specie != "Time":
                last_step[specie] = concentration[-1]

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

        """

        Args: poisson distribution parameter (lambda)
            lambdas:

        Returns:
            num_reaction_: a dictionary contains number of times ach reaction occurred in the time interval (tau).

        """

        num_reaction_ = {}
        for reaction, lambda_ in lambdas.items():
            num_reaction_[reaction] = np.random.poisson(lambda_)

        return num_reaction_

    def update(self, species, model, num_reaction, tau):

        """
        Args:
            species: a dictionary in which the calculated concentrations of the species are stored.
            model: a class created by "biostoch.model.Model".
            num_reaction: an integer value that indicates the number of reaction in the system.
            tau: a float or an integer, calculated tau.
        Returns:
            species: a dictionary in which the calculated concentrations of the species are stored.
        """

        species["Time"].append(species["Time"][-1] + tau)

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
            species[component].append(species[component][-1] + value)

        return species

    def simulate(self):

        """Runs the simulation"""

        start_simulation = time.time()

        species, parameters = self.initialize_parameters(
            model=self.model,
            start=self.start
        )

        step = 2
        while species["Time"][-1] <= self.stop:

            propensity_sum, propensities_ = self.compute_propensity_sum(
                species=species,
                parameters=parameters,
                propensities=self.model.rates_
            )

            if propensity_sum == 0 and self.steady_state:
                print(f"Simulation reached steady state (iteration: {step}). No further changes are occurring.")
                break

            if self.call_tau:
                tau = self.compute_tau(
                    species=species,
                    model=self.model,
                    epsilon=self.epsilon
                )
            else:
                tau = self.tau

            lambdas = self.compute_lambdas(
                species=species,
                parameters=self.model.params,
                propensities=self.model.rates_,
                tau=tau
            )

            num_reaction = self.num_reaction(
                lambdas=lambdas
            )

            species = self.update(
                species=species,
                model=self.model,
                num_reaction=num_reaction,
                tau=tau
            )

            step += 1
            if step == self.max_epochs:
                print(f"Simulation reached the maximum iteration (max_epochs={self.max_epochs})!")
                break

        self.species = species
        self.parameters = parameters
        stop_simulation = time.time()
        self.time["Simulation Duration"] = stop_simulation - start_simulation





