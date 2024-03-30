import numpy as np
import time


class GillespieSimulator(object):

    """ Simulation using Gillespie's Stochastic Simulation Algorithm (SSA) """

    def __init__(
        self,
        model=None,
        start=0,
        stop=10,
        max_epochs=100,
        seed=42,
        steady_state=False,
        gamma=1e-30,
        model_name="Stochastic Simulation Algorithm",
        **kwargs
    ):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.steady_state = steady_state
        self.gamma = gamma
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
            gamma: a small float value used to prevent zero division when calculating tau.
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

    def initialize_parameters(self, model, start):

        """
        Initializes the model species dictionary

        Args:
           model: a class created by "biostoch.model.Model"
                  contains all necessary information used in simulation with SSA.
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


    def compute_tau(self, propensity_sum, gamma):

        """
        Computes  the time until the next reaction event occurs (tau)
        based on the sum of reaction propensities in the system.

        Args:
            propensity_sum: a float value, sum of the propensities.
            gamma: a small float value used to prevent zero division when calculating tau.

        Returns:
            tau: a float or an integer, calculated tau.

        """

        tau = np.random.exponential(scale=1 / (propensity_sum + gamma))

        return tau

    def update(self, species, model, reaction, num_reaction, propensities, tau):

        """

        Args:
            species: a dictionary in which the calculated concentrations of the species are stored.
            model: a class created by "biostoch.model.Model".
            reaction: a float value that indicates which reaction is taking place.
            num_reaction: an integer value that indicates the number of reaction in the system.
            propensities: a dictionary contains the propensity values of the reactions.
            tau: a float or an integer, calculated tau.

        Returns:
            species: a dictionary in which the calculated concentrations of the species are stored.
        """

        species["Time"].append(species["Time"][-1] + tau)

        for i in range(num_reaction):
            reaction_name = model.react_names[i]
            if i == 0:
                if reaction <= propensities[reaction_name]:
                    split_reaction = model.reacts_[reaction_name].split()
                    index = [index for index, value in enumerate(split_reaction) if value == '->']
                    if len(index) > 1:
                        print(f"Each reaction should have exactly one '->', but there are more than one in the {reaction_name}.")

                    components_ = []
                    for j in range(index[0]):
                        if split_reaction[j] in model.components:
                            components_.append(split_reaction[j])
                            species[split_reaction[j]].append(species[split_reaction[j]][-1] - 1)
                    for k in range(index[0] + 1, len(split_reaction)):
                        if split_reaction[k] in model.components:
                            components_.append(split_reaction[k])
                            species[split_reaction[k]].append(species[split_reaction[k]][-1] + 1)
                    for specie_ in species.keys():
                        if specie_ not in components_ and specie_ != "Time":
                            species[specie_].append(species[specie_][-1])

            else:

                reaction_name_ = model.react_names[i-1]
                keys_to_sum = model.react_names[:i+1]
                sum_propensities_ = sum(propensities[react_name_] for react_name_ in keys_to_sum)

                if reaction > propensities[reaction_name_] and reaction <= sum_propensities_:

                    split_reaction = model.reacts_[reaction_name].split()
                    index = [index for index, value in enumerate(split_reaction) if value == '->']
                    if len(index) > 1:
                        print(f"Each reaction should have exactly one '->', but there are more than one in the {reaction_name}.")

                    components_ = []
                    for j in range(index[0]):
                        if split_reaction[j] in model.components:
                            components_.append(split_reaction[j])
                            species[split_reaction[j]].append(species[split_reaction[j]][-1] - 1)
                    for k in range(index[0] + 1, len(split_reaction)):
                        if split_reaction[k] in model.components:
                            components_.append(split_reaction[k])
                            species[split_reaction[k]].append(species[split_reaction[k]][-1] + 1)
                    for specie_ in species.keys():
                        if specie_ not in components_ and specie_ != "Time":
                            species[specie_].append(species[specie_][-1])

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
                propensities=self.model.rates_,
                species=species,
                parameters=parameters
            )

            if propensity_sum == 0 and self.steady_state:
                print(f"Simulation reached steady state (iteration: {step}). No further changes are occurring.")
                break

            tau = self.compute_tau(
                propensity_sum=propensity_sum,
                gamma=self.gamma
            )

            random_number = np.random.uniform(low=0, high=1)
            num_reactions = len(self.model.reacts_)
            reaction = propensity_sum * random_number

            species = self.update(
                species=species,
                model=self.model,
                reaction=reaction,
                num_reaction=num_reactions,
                propensities=propensities_,
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



