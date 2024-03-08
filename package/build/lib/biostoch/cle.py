import numpy as np
import time



class ChemicalLangevin(object):

    """ Simulation using Chemical Langevin Equation """

    def __init__(
        self,
        model=None,
        start=0.0,
        stop=10.0,
        max_epochs=100,
        seed=42,
        steady_state=None,
        **kwargs
    ):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.steady_state = steady_state

        self.tau = (self.stop - self.start) / self.max_epochs

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)
        else:
            raise ValueError("Before simulating a model, please ensure that you have instantiated the biostoch.model.Model() object.")

        self.model_name = "Chemical Langevin Equation"
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
            if getattr(model, specie) != 0:
                species[specie][0] = getattr(model, specie)
            else:
                species[specie][0] = getattr(model, specie)

        for parameter in self.model.params:
            parameters[parameter] = getattr(model, parameter)

        return species, parameters


    def compute_change(self, model, species, tau, step):

        changes = {}
        terms = {}

        for parameters, value in model.params.items():
            terms[parameters] = value

        for specie, concentration in species.items():
            terms[specie] = concentration[step - 1]

        for reaction, rate in model.rates_.items():
            changes[reaction] = eval(rate, terms) * tau

        return changes, terms


    def compute_noise(self, model, terms, tau):

        noises = {}

        for reaction, rate in model.rates_.items():
            random_number = np.random.normal()
            noises[reaction] = (tau**.5) * ((eval(rate, terms))**.5) * random_number

        return noises


    def compute_changes(self, noises, changes):

        changes_ = {}

        for reaction in changes.keys():
            changes_[reaction] = noises[reaction] + changes[reaction]

        return changes_


    def update(self, species, model, changes_, step, tau):

        species["Time"][step] = species["Time"][step - 1] + tau

        for reaction, formula in model.reacts_.items():
            split_formula = formula.split()
            index = [index for index, value in enumerate(split_formula) if value == '->']

            if len(index) != 1:
                print(f"Error: Each reaction should have exactly one '->', but there are {len(index)} in {reaction}.")

        component_reaction = {}

        for component in model.components:

            num_reaction_ = sum([changes_[reaction_] * model.coeffs_[reaction_][component] for reaction_ in model.react_names if component in model.react_sps[reaction_]])
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

        species, parameters = self.initialize_parameters(
            model=self.model,
            start=self.start,
            max_epochs=self.max_epochs
        )

        step = 1
        while species["Time"][step] < self.stop and step < self.max_epochs:

            changes, terms = self.compute_change(
                model=self.model,
                species=species,
                tau=self.tau,
                step=step
            )

            noises = self.compute_noise(
                model=self.model,
                terms=terms,
                tau=self.tau
            )

            changes_ = self.compute_changes(
                noises=noises,
                changes=changes
            )

            species = self.update(
                species=species,
                model=self.model,
                changes_=changes_,
                step=step,
                tau=self.tau
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




