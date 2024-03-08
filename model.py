import numpy as np
import matplotlib.pyplot as plt
import time


class Model(object):

    """ Define biological or chemical system """

    def __init__(self, signs=None):

        if not signs:
            self.signs = ["+", "->", "*", "-", "**"]

    def __repr__(self):
        if hasattr(self, 'params'):
            params = self.params
        else:
            params = None

        if hasattr(self, 'components'):
            components = self.components
        else:
            components = None

        if hasattr(self, 'ROC_'):
            ROC_ = self.ROC_
        else:
            ROC_ = None

        if hasattr(self, 'react_names'):
            react_names = self.react_names
        else:
            react_names = None

        if hasattr(self, 'reacts_'):
            reacts_ = self.reacts_
        else:
            reacts_ = None

        if hasattr(self, 'coeffs_'):
            coeffs_ = self.coeffs_
        else:
            coeffs_ = None

        if hasattr(self, 'react_sps'):
            react_sps = self.react_sps
        else:
            react_sps = None

        if hasattr(self, 'rates_'):
            rates_ = self.rates_
        else:
            rates_ = None

        return (f"Model: {self.signs}, {params}, {components}, {ROC_}, "
                f"{react_names}, {reacts_}, {coeffs_}, {react_sps}, {rates_}")

    def reset(self):

        self.params = None
        self.components = None
        self.ROC_ = None
        self.react_names = None
        self.reacts_ = None
        self.coeffs_ = None
        self.reset_sps = None
        self.rates = None

    def parameters(self, params):
        """
        params: A dictionary containing all rate constants in the system.
                Each key represents a rate constant (e.g., K1, K2, ..., Kn),
                and each corresponding value represents the value of that rate constant.
        """
        if not isinstance(params, dict):
            raise TypeError("params must be a dictionary.")

        self.params = params
        setattr(self, "param_names", list(params.keys()))

        for key, val in params.items():
            setattr(self, key, val)


    def species(self, components, rate_change):

        """
        components: A dictionary containing all species in the system.
                   Each key represents a specie (e.g. A, B, ...), and
                   each corresponding value represents the initial concentration.
        rate_change: A dictionary containing rate of change for each component.
        """

        if not isinstance(components, dict):
            raise TypeError("components must be a dictionary.")

        if not isinstance(rate_change, dict):
            raise TypeError("rate_change must be a dictionary.")

        setattr(self, "components", list(components.keys()))

        for key, val in components.items():
            setattr(self, key, val)

        ROC_ = dict()
        for key, val in rate_change.items():

            roc = []
            val = val.replace('+', ' + ').replace('->', ' -> ').replace('*', ' * ').replace('-', ' - ').replace('**', ' ** ')
            val = val.split()

            for v in val:
                if v in self.components or v in self.params.keys():
                    roc.append(v)

                elif v in self.signs:
                    roc.append(v)
                else:
                    try:
                        roc.append(float(v))
                    except ValueError:
                        print(f"This part of the rate_change ({v}) is not valid!")

            roc = " ".join([str(r) for r in roc])

            ROC_[key] = roc

        setattr(self, "ROC_", ROC_)

    def reactions(self, reacts, rates):
        """
        reactions: A dictionary containing reaction equation for each reaction.
                   Each key represents the reaction name, and each value represents the reaction equation.
                   exp:
                       {"reaction1": "A + B -> C", ...}

        rates: A dictionary containing rate equations (or propensity functions) for each reaction.
               Each key represents the reaction name, and each value represents the rate equation of that reaction.
               exp:
                   {reaction1: "K1 * A * B", ...}
        """

        if not isinstance(reacts, dict):
            raise TypeError("reacts must be a dictionary.")

        if not isinstance(rates, dict):
            raise TypeError("rates must be a dictionary.")

        setattr(self, "react_names", list(reacts.keys()))

        if not self.params or not self.components:
            raise ValueError("Please define Species and Parameters before Reactions!")

        reacts_ = dict()
        coefficients = dict()
        react_sps = dict()
        for reaction, formula in reacts.items():

            react = []
            react_eq = formula.replace('+', ' + ').replace('->', ' -> ').replace('*', ' * ').replace('**', ' ** ')
            components = react_eq.split()

            if '->' not in formula:
                raise ValueError(f"Missing '->' in the reaction equation for {reaction}.")

            if formula.count('->') != 1:
                raise ValueError(
                    f"Each reaction should have exactly one '->', but there are {formula.count('->')} in {reaction}.")

            for component in components:
                if component in self.components or component in self.params:
                    react.append(component)

                elif component in self.signs:
                    react.append(component)
                else:
                    try:
                        react.append(float(component))
                    except ValueError:
                        print(f"Unexpected component '{component}' in the reaction equation for {reaction}.")
            comp = []
            for component in components:
                if component in self.components:
                    comp.append(component)

            react = " ".join([str(r) for r in react])
            reacts_[reaction] = react
            react_sps[reaction] = comp

            coefficient = dict()
            index = [index for index, value in enumerate(components) if value == '->']
            if len(index) > 1:
                print(f"Each reaction should have exactly one '->', but there are more than one in the {reaction}.")

            for i in range(index[0]):
                if components[i] in self.components:
                    if components[i-1] not in self.signs and components[i-1] not in self.components:
                        coefficient[components[i]] = - eval(components[i-1])
                    else:
                        coefficient[components[i]] = - 1
            for j in range(index[0]+1, len(components)):
                if components[j] in self.components:
                    if components[j-1] not in self.signs and components[j-1] not in self.components:
                        coefficient[components[j]] = eval(components[j - 1])
                    else:
                        coefficient[components[j]] = 1

            coefficients[reaction] = coefficient


        setattr(self, "reacts_", reacts_)
        setattr(self, "coeffs_", coefficients)
        setattr(self, "react_sps", react_sps)


        rates_ = dict()

        for reaction, rate in rates.items():
            rate_equation = []
            rate_eq = rate.replace('*', ' * ').replace('+', ' + ').replace('->', ' -> ').replace('**', ' ** ')
            components = rate_eq.split()

            for component in components:
                if component in self.components or component in self.params:
                    rate_equation.append(component)
                elif component in self.signs:
                    rate_equation.append(component)
                else:
                    try:
                        rate_equation.append(float(component))
                    except ValueError:
                        print(f"This component {component} is not valid!")

            rate_equation = " ".join([str(rate) for rate in rate_equation])

            try:
                if reaction in self.reacts_.keys():
                    rates_[reaction] = rate_equation
            except ValueError:
                print("Reactions and rates do not match!")

        setattr(self, "rates_", rates_)




class EulerSimulator(object):

    """ Simulation using Euler method """

    def __init__(
        self,
        model=None,
        start=0,
        stop=10,
        epochs=1000,
        seed=42,
        **kwargs
    ):

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

    def initialize_parameters(self, model, start, stop, epochs):

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

        species, parameters = self.initialize_parameters(model=self.model, start=self.start, stop=self.stop,
                                                         epochs=self.epochs)

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
        **kwargs
        ):

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

    def initialize_parameters(self, model, start, stop, epochs):

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

        species, parameters = self.initialize_parameters(model=self.model, start=self.start, stop=self.stop,
                                                         epochs=self.epochs)

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




class GillespieSimulator(object):
    """ Simulation using Stochastic Simulation Algorithm """

    def __init__(
        self,
        model=None,
        start=0,
        stop=10,
        max_epochs=100,
        seed=42,
        steady_state=None,
        gamma=1e-30,
        **kwargs
    ):

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





class TauLeaping(object):

    """ Simulation using Tau-Leaping method """

    def __init__(
        self,
        model=None,
        start=0.0,
        stop=10.0,
        max_epochs=100,
        seed=42,
        steady_state=None,
        epsilon=0.03,
        call_tau=None,
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



