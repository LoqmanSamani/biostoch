import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class Model(object):

    def __init__(self, signs=None):

        if not signs:
            self.signs = ["+", "->", "*", "-"]

    def parameters(self, params):
        """
        params: A dictionary containing all rate constants in the system.
                Each key represents a rate constant (e.g., K1, K2, ..., Kn),
                and each corresponding value represents the value of that rate constant.
        """
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

        setattr(self, "components", list(components.keys()))

        for key, val in components.items():
            setattr(self, key, val)

        ROC_ = dict()
        for key, val in rate_change.items():

            roc = []
            val = val.replace('+', ' + ').replace('->', ' -> ').replace('*', ' * ').replace('-', ' - ')
            val = val.split()

            for v in val:
                if v in self.components or v in self.params.keys():
                    roc.append("self." + v)

                elif v in self.signs:
                    roc.append(v)
                else:
                    try:
                        roc.append(float(v))
                    except ValueError:
                        print(f"This part of the ROC-equation ({v}) is not valid!")

            roc = " ".join([str(r) for r in roc])

            ROC_[key] = roc

        setattr(self, "ROC_", ROC_)

    def reactions(self, reacts, rates=None):
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
        setattr(self, "react_names", list(reacts.keys()))

        if not self.params or not self.components:
            raise ("Please, first define Species and Parameters, then Reactions!")

        reacts_ = dict()

        for reaction, formula in reacts.items():
            react = []
            react_eq = formula.replace('+', ' + ').replace('->', ' -> ').replace('*', ' * ')
            components = react_eq.split()

            for component in components:
                if component in self.components or component in self.params:
                    react.append(component)
                elif component in self.signs:
                    react.append(component)
                else:
                    try:
                        react.append(float(component))
                    except ValueError:
                        print(f"This component {component} is not valid!")

            react = " ".join([str(r) for r in react])
            reacts_[reaction] = react

        setattr(self, "reacts_", reacts_)

        if rates:

            rates_ = dict()

            for reaction, rate in rates.items():
                rate_equation = []
                rate_eq = rate.replace('*', ' * ').replace('+', ' + ').replace('->', ' -> ')
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


class ODE(object):
    def __init__(self, model=None, start=0, stop=10, epochs=1000, seed=42, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.epochs = epochs
        self.seed = seed

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)

        self.species = None
        self.parameters = None

    def param_init(self, model, start, stop, epochs):

        species = {}
        parameters = {}
        species["Time"] = np.linspace(start, stop, epochs)

        for specie in model.components:
            species[specie] = np.zeros(epochs)
            species[specie][0] = getattr(model, specie)
        for parameter in self.model.params:
            parameters[parameter] = getattr(model, parameter)

        return species, parameters

    def der_species(self, model):

        der_sp = {}
        for specie, rate in model.ROC_.items():
            der_sp[specie] = eval(rate)

        return der_sp


    def simulate(self):

        species, parameters = self.param_init(
            model=self.model,
            start=self.start,
            stop=self.stop,
            epochs=self.epochs
        )

        der_sp = self.der_species(model=self.model)

        dt = species["Time"][3] - species["Time"][2]

        for i in range(1, self.epochs):
            for specie in species.keys():
                if specie != "Time":
                    species[specie][i] = species[specie][i - 1] + (der_sp[specie] * dt)

        self.species = species
        self.parameters = parameters


class SSA(object):
    def __init__(self, model=None, start=0, stop=10, epochs=None, max_epochs=None, seed=42, alpha=100, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.epochs = epochs
        self.max_epochs = max_epochs
        self.seed = seed
        self.alpha = alpha

        if self.model:
            model_attributes = vars(self.model)
            self.__dict__.update(model_attributes)

        self.species = None
        self.parameters = None

    def param_init(self, model, start, stop, max_epochs, alpha):

        species = {}
        parameters = {}

        if max_epochs:
            epochs = max_epochs
        else:
            epochs = (stop - start) * alpha

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
            last_step[key] = val[step]
        for key, val in params.items():
            last_step[key] = val

        for reaction, propensity in propensities.items():
            propensity = eval(propensity, last_step)
            propensity_sum += propensity
            props[reaction] = propensity

        return propensity_sum, props

    def simulate(self):

        num_reacts = len(self.model.reacts_)

        species, parameters = self.param_init(
            model=self.model,
            start=self.start,
            stop=self.stop,
            max_epochs=self.max_epochs,
            alpha=self.alpha
        )

        step = 0
        ind = 1
        while species["Time"][-1] < self.stop:

            a_sum, props = self.propensity_sum(
                step=step,
                propensities=self.model.rates_,
                species=species,
                params=parameters
            )

            tau = np.random.exponential(scale=1/a_sum)
            rand = np.random.uniform(low=0, high=1)

            react = tau * rand
            for i in range(num_reacts):
                react_name = self.model.react_names[i]
                if i == 0:
                    if react <= props[react_name]:
                        sp = self.model.reacts_[react_name].split()
                        index = [index for index, value in enumerate(sp) if value == '->']
                        for j in range(index[0]):
                            if sp[j] in self.model.components:
                                species[sp[j]][ind] = species[sp[j]][ind-1] + 1
                        for k in range(index[0], len(sp)):
                            if sp[k] in self.model.components:
                                species[sp[k]][ind] = species[sp[k]][ind-1] - 1
                    else:
                        pass
                    











                   #elif react > props[i] and react <= np.sum(props[:i+1])




























        







m = Model()
m.parameters({"K1": 3, "K2": 5})

m.species({"A": 1, "B": 100, "C": 5}, {"A": "-2*A*K1", "B": "-3*B+4*A*K2", "C":"-C"})

m.reactions({"reaction1": "A + B -> C", "reaction2": "C -> A + B"},
           {"reaction1": "2.0 * A * B", "reaction2": "4.0*C"})


print("Parameters: ", m.params)
print("Species: ", m.components)
print("Rate of Reactions: ", m.rates_)
print("Reactions: ", m.reacts_)
print("Rate of Changes:", m.ROC_)
print("param names", m.param_names)
print("react names", m.react_names)


print("K1: ", m.K1)
print("K2: ", m.K2)

print("A: ", m.A)
print("B: ", m.B)
print("C: ", m.C)

model = SSA(m)
print(model.A)
print(model.B)

s, p = model.simulate()
print(s)
print(p)
print(len(s))
print(len(s["A"]))
#s= ODE(m)

#s.simulate()

#print(s.species)
#print(s.parameters)
#print(s.components)
#print(s.params)




