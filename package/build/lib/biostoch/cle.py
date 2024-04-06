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
        steady_state=False,
        model_name="Chemical Langevin Equation",
        **kwargs
    ):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.steady_state = steady_state
        self.model_name = model_name

        self.tau = (self.stop - self.start) / self.max_epochs

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
            model_name: the name of the simulation method "Tau-Leaping Algorithm".
            **kwargs: a special parameter that allows passing additional keyword arguments to the function.
            tau: a float or an integer that represents the time step size.
            
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

    def compute_change(self, model, species, tau):

        """

        Args:
            model: a class created by "biostoch.model.Model"
            species: a dictionary in which the calculated concentrations of the species are stored.
            tau: a float or an integer, calculated tau.

        Returns:
            changes: a dictionary stores the computed changes in the concentrations
                     of species due to each reaction during the time step (tau)(without noise).
            terms: a dictionary is used to collect the terms needed for evaluating the rates of reactions.
        """

        changes = {}
        terms = {}

        for parameters, value in model.params.items():
            terms[parameters] = value

        for specie, concentration in species.items():
            terms[specie] = concentration[-1]

        for reaction, rate in model.rates_.items():
            changes[reaction] = eval(rate, terms) * tau

        return changes, terms

    def compute_noise(self, model, terms, tau):

        """

        Args:
            model: a class created by "biostoch.model.Model"
            terms: a dictionary is used to collect the terms needed for evaluating the rates of reactions.
            tau: a float or an integer, calculated tau.

        Returns:
            noises: a dictionary contains computed noise terms for each reaction in the system.

        """

        noises = {}

        for reaction, rate in model.rates_.items():
            random_number = np.random.normal()
            noises[reaction] = (tau**.5) * ((eval(rate, terms))**.5) * random_number

        return noises

    def compute_changes(self, noises, changes):

        """

        Args:
            noises: a dictionary contains computed noise terms for each reaction in the system.
            changes: a dictionary stores the computed changes in the concentrations
                     of species due to each reaction during the time step (tau) (without noise).

        Returns:
            changes: a dictionary stores the computed changes in the concentrations
                     of species due to each reaction during the time step (tau) (witt noise).
        """

        changes_ = {}

        for reaction in changes.keys():
            changes_[reaction] = noises[reaction] + changes[reaction]

        return changes_

    def update(self, species, model, changes_, tau):

        """
        Args:
            species: a dictionary in which the calculated concentrations of the species are stored.
            model: a class created by "biostoch.model.Model".
            changes_: a dictionary stores the computed changes in the concentrations
                     of species due to each reaction during the time step (tau).
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

            num_reaction_ = sum([changes_[reaction_] * model.coeffs_[reaction_][component] for reaction_ in model.react_names if component in model.react_sps[reaction_]])
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

            changes, terms = self.compute_change(
                model=self.model,
                species=species,
                tau=self.tau
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
                tau=self.tau
            )

            step += 1
            if step == self.max_epochs:
                print(f"Simulation reached the maximum iteration (max_epochs={self.max_epochs})!")
                break

        self.species = species
        self.parameters = parameters
        stop_simulation = time.time()
        self.time["Simulation Duration"] = stop_simulation - start_simulation

