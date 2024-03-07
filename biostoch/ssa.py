import numpy as np
import model
import time
import matplotlib.pyplot as plt



class GillespieSimulator(object):
    def __init__(self, model=None, start=0, stop=10, max_epochs=100,
                 seed=42, steady_state=None, gamma=1e-30, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.steady_state = steady_state
        self.gamma = gamma

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)
        else:
            raise ValueError("Before simulating a model, please ensure that you have instantiated the biostoch.model.Model() object.")

        self.model_name = "Stochastic Simulation Algorithm"
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


    def compute_propensity_sum(self, step, propensities, species, parameters):

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


    def compute_tau(self, propensity_sum, gamma):

        tau = np.random.exponential(scale=1 / (propensity_sum + gamma))

        return tau


    def update(self, species, model, reaction, num_reaction, propensities, step, tau):

        species["Time"][step] = species["Time"][step - 1] + tau

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
                            species[split_reaction[j]][step] = species[split_reaction[j]][step - 1] - 1
                    for k in range(index[0] + 1, len(split_reaction)):
                        if split_reaction[k] in model.components:
                            components_.append(split_reaction[k])
                            species[split_reaction[k]][step] = species[split_reaction[k]][step - 1] + 1
                    for specie_ in species.keys():
                        if specie_ not in components_ and specie_ != "Time":
                            species[specie_][step] = species[specie_][step - 1]

            else:

                reaction_name_ = model.react_names[i - 1]
                keys_to_sum = model.react_names[:i + 1]
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
                            species[split_reaction[j]][step] = species[split_reaction[j]][step - 1] - 1
                    for k in range(index[0] + 1, len(split_reaction)):
                        if split_reaction[k] in model.components:
                            components_.append(split_reaction[k])
                            species[split_reaction[k]][step] = species[split_reaction[k]][step - 1] + 1
                    for specie_ in species.keys():
                        if specie_ not in components_ and specie_ != "Time":
                            species[specie_][step] = species[specie_][step - 1]

        return species


    def resize_species(self, species, step):

        if step >= len(species["Time"]):

            new_max_steps = len(species["Time"]) * 2

            for specie, concentration in species.items():
                pad_width = (0, new_max_steps - len(specie))
                species[specie] = np.pad(specie, pad_width, mode='constant')

        return species


    def final_resize_species(self, species, final_step):

        for specie in species.keys():
            species[specie] = species[specie][:final_step]

        return species


    def simulate(self):

        start_simulation = time.time()

        species, parameters = self.initialize_parameters(
            model=self.model,
            start=self.start,
            max_epochs=self.max_epochs
        )

        step = 1
        while species["Time"][step-1] < self.stop:

            propensity_sum, propensities_ = self.compute_propensity_sum(
                step=step,
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





