import numpy as np
import matplotlib.pyplot as plt



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
                    roc.append(v)

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
        coefficients = dict()
        react_sps = dict()
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
                        coefficient[components[i]] = -1
            for j in range(index[0]+1, len(components)):
                if components[j] in self.components:
                    if components[j-1] not in self.signs and components[j-1] not in self.components:
                        coefficient[components[j]] =  eval(components[j - 1])
                    else:
                        coefficient[components[j]] = 1

            coefficients[reaction] = coefficient


        setattr(self, "reacts_", reacts_)
        setattr(self, "coeffs_", coefficients)
        setattr(self, "react_sps", react_sps)


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
            der_sp[specie] = rate

        return der_sp

    def simulate(self):

        species, parameters = self.param_init(
            model=self.model,
            start=self.start,
            stop=self.stop,
            epochs=self.epochs
        )

        dt = species["Time"][3] - species["Time"][2]

        for i in range(1, self.epochs):

            der_sp = self.der_species(model=self.model)

            for specie in species.keys():
                if specie != "Time":
                    rate = ""
                    split_rate = der_sp[specie].split()
                    for s in split_rate:
                        if s in self.model.params.keys():
                            rate += " " + str(self.model.params[s])
                        elif s in ["+", "-", "*", "/"]:
                            rate += " " + s
                        elif s in self.model.components:
                            rate += " " + str(species[s][i-1])
                        else:
                            ValueError(f"There is a problem with your model ({s}).")

                    species[specie][i] = species[specie][i - 1] + (eval(rate) * dt)

        self.species = species
        self.parameters = parameters


class SSA(object):
    def __init__(self, model=None, start=0, stop=10, max_epochs=1000,
                 seed=42, alpha=100, steady_state=None, gamma=1e-30, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.alpha = alpha
        self.steady_state = steady_state
        self.gamma = gamma

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
            if key != "Time":
                last_step[key] = val[step-1]
        for key, val in params.items():
            last_step[key] = val

        for reaction, propensity in propensities.items():
            propensity = eval(propensity, last_step)
            propensity_sum += propensity
            props[reaction] = propensity

        return propensity_sum, props

    def calculate_dt(self, a_sum, gamma):

        dt = np.random.exponential(scale=1 / (a_sum + gamma))

        return dt

    def update(self, species, model, react, num_reacts, props, step, dt):

        species["Time"][step] = species["Time"][step - 1] + dt

        for i in range(num_reacts):
            react_name = model.react_names[i]
            if i == 0:
                if react <= props[react_name]:
                    sp = model.reacts_[react_name].split()
                    index = [index for index, value in enumerate(sp) if value == '->']
                    if len(index) > 1:
                        print(f"Each reaction should have exactly one '->', but there are more than one in the {react_name}.")

                    comp = []
                    for j in range(index[0]):
                        if sp[j] in model.components:
                            comp.append(sp[j])
                            species[sp[j]][step] = species[sp[j]][step - 1] - 1
                    for k in range(index[0] + 1, len(sp)):
                        if sp[k] in model.components:
                            comp.append(sp[k])
                            species[sp[k]][step] = species[sp[k]][step - 1] + 1
                    for sps in species.keys():
                        if sps not in comp and sps != "Time":
                            species[sps][step] = species[sps][step - 1]

            else:

                name = model.react_names[i-1]
                keys_to_sum = model.react_names[:i + 1]
                sum_props = sum(props[key] for key in keys_to_sum)
                if react > props[name] and react <= sum_props:
                    sp = model.reacts_[react_name].split()
                    index = [index for index, value in enumerate(sp) if value == '->']
                    if len(index) > 1:
                        print(f"Each reaction should have exactly one '->', but there are more than one in the {react_name}.")

                    comp = []
                    for j in range(index[0]):
                        if sp[j] in model.components:
                            comp.append(sp[j])
                            species[sp[j]][step] = species[sp[j]][step - 1] - 1
                    for k in range(index[0] + 1, len(sp)):
                        if sp[k] in model.components:
                            comp.append(sp[k])
                            species[sp[k]][step] = species[sp[k]][step - 1] + 1
                    for sps in species.keys():
                        if sps not in comp and sps != "Time":
                            species[sps][step] = species[sps][step - 1]

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
            stop=self.stop,
            max_epochs=self.max_epochs,
            alpha=self.alpha
        )

        step = 1
        while species["Time"][step-1] < self.stop and step-1 < self.max_epochs:

            a_sum, props = self.propensity_sum(
                step=step,
                propensities=self.model.rates_,
                species=species,
                params=parameters
            )
            if a_sum == 0 and self.steady_state:
                print(f"Simulation reached steady state (iteration: {step}). No further changes are occurring.")
                break

            dt = self.calculate_dt(
                a_sum=a_sum,
                gamma=self.gamma
            )
            rand = np.random.uniform(low=0, high=1)
            num_reacts = len(self.model.reacts_)
            react = a_sum * rand

            species = self.update(
                species=species,
                model=self.model,
                react=react,
                num_reacts=num_reacts,
                props=props,
                step=step,
                dt=dt
            )

            step += 1

            species = self.change_size(
                species=species,
                step=step
            )

        species = self.resize_species(
            species=species,
            final_step=step
        )

        self.species = species
        self.parameters = parameters


class TauLeaping(object):
    def __init__(self, model=None, start=0, stop=None, max_epochs=100,
                 seed=42, alpha=100, steady_state=None, tau=None, epsilon=0.03, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.max_epochs = max_epochs
        self.seed = seed
        self.alpha = alpha
        self.steady_state = steady_state
        self.tau = tau
        self.epsilon = epsilon

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
        elif stop:
            epochs = (stop - start) * alpha
        else:
            epochs = (10 * start) * alpha

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
        X = X = np.array([species[con][step-1] for con in species.keys() if con != "Time"])
        v = []

        for key, val in model.coeffs_.items():
            d = []
            for sp in model.react_sps[key]:
                d.append(val[sp])
            v.append(d)

        v = np.array(v)
        R = []
        for p in model.param_names:
            R.append(model.params[p])
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
                last_step[key] = val[step-1]
        for key, val in params.items():
            last_step[key] = val

        lambdas = {}
        for react, prop in propensities.items():
            lambdas[react] = eval(prop, last_step) * tau

        return lambdas

    def num_reacts(self, lambdas):
        n_reacts = {}
        for react, lam in lambdas.items():
            n_reacts[react] = np.random.poisson(lam)

        return n_reacts

    def update(self, species, model, num_reacts, step, tau):

        species["Time"][step] = species["Time"][step - 1] + tau

        for react, prop in model.reacts_.items():
            split_prop = prop.split()
            index = [index for index, value in enumerate(split_prop) if value == '->']
            if len(index) > 1:
                print(f"Each reaction should have exactly one '->', but there are more than one in the {react}.")

            for i in range(index[0]):
                if split_prop[i] in model.components:
                    species[split_prop[i]][step] = species[split_prop[i]][step-1] - num_reacts[react]
            for j in range(index[0], len(split_prop)):
                if split_prop[j] in model.components:
                    species[split_prop[j]][step] = species[split_prop[j]][step-1] + num_reacts[react]

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
            stop=self.stop,
            max_epochs=self.max_epochs,
            alpha=self.alpha
        )

        step = 1
        while step-1 < self.max_epochs:

            a_sum, props = self.propensity_sum(
                step=step,
                propensities=self.model.rates_,
                species=species,
                params=parameters
            )

            if a_sum == 0 and self.steady_state:
                print(f"Simulation reached steady state (iteration: {step}). No further changes are occurring.")
                break

            if self.tau:
                tau = self.tau
            else:

                tau = self.calculate_tau(
                    species=species,
                    model=self.model,
                    step=step,
                    epsilon=self.epsilon
                )

            lambdas = self.calculate_lambda(
                propensities=self.model.rates_,
                species=species,
                params=self.model.params,
                tau=tau,
                step=step
            )
            num_reacts = self.num_reacts(
                lambdas=lambdas
            )

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

        species = self.resize_species(
            species=species,
            final_step=step
        )

        self.species = species
        self.parameters = parameters


class CLE(object):
    def __init__(self, model=None, start=0, stop=10, epochs=None, max_epochs=None,
                     seed=42, alpha=100, steady_state=None, tau=None, epsilon=0.3, **kwargs):

        self.model = model
        self.start = start
        self.stop = stop
        self.epochs = epochs
        self.max_epochs = max_epochs
        self.seed = seed
        self.alpha = alpha
        self.steady_state = steady_state
        self.tau = tau
        self.epsilon = epsilon

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

    def der_species(self, model):

        der_sp = {}
        for specie, rate in model.ROC_.items():
            der_sp[specie] = eval(rate)

        return der_sp

    def compute_tau(self, X, params, epsilon):

        eXi = epsilon * X
        gi = max(params.values())
        mu_X = sum(abs(rate) for rate in params.values())
        sigma2_X = sum(rate ** 2 for rate in params.values())

        tau1 = max(np.maximum(eXi / gi, 1) / np.abs(mu_X))
        tau2 = max(np.maximum(eXi, gi) ** 2 / sigma2_X)
        tau = min(tau1, tau2)

        return tau

    def update(self, species, der_sp, tau, step):

        species["Time"][step] = species["Time"][step - 1] + tau

        for specie in species.keys():
            if specie != "Time":
                rate = der_sp[specie] * tau
                noise = np.sqrt(der_sp[specie] * np.random.normal(loc=0, scale=1) * tau)
                d_specie = rate + noise
                species[specie][step] = species[specie][step - 1] + d_specie

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
            stop=self.stop,
            max_epochs=self.epochs,
            alpha=self.alpha
        )

        der_sp = self.der_species(
            model=self.model
        )

        if self.max_epochs:
            max_epochs = self.max_epochs
        else:
            max_epochs = len(list(species["Time"]))

        step = 1
        while species["Time"][step - 1] < self.stop and step - 1 < max_epochs:

            if self.tau:
                tau = self.tau
            else:
                X = np.array([species[con][step - 1] for con in species.keys() if con != "Time"])
                tau = self.compute_tau(
                    X=X,
                    params=self.model.params,
                    epsilon=self.epsilon
                )

            species = self.update(
                species=species,
                der_sp=der_sp,
                tau=tau,
                step=step
            )

            step += 1

            species = self.change_size(
                species=species,
                step=step
            )

        species = self.resize_species(
            species=species,
            final_step=step
        )

        self.species = species
        self.parameters = parameters








m = Model()
m.parameters({"K1": 0.1, "K2": 0.05})

m.species({"A": 100, "B": 0}, {"A": "K2 * B - K1 * A", "B": "K1 * A - K2 * B"})

m.reactions({"reaction1": "A -> B", "reaction2": "B -> A"},
           {"reaction1": "K1 * A", "reaction2": "K2 * B"})


print("Parameters: ", m.params)
print("Species: ", m.components)
print("Rate of Reactions: ", m.rates_)
print("Reactions: ", m.reacts_)
print("Rate of Changes:", m.ROC_)
print("param names", m.param_names)
print("react names", m.react_names)
print("Coefficients: ", m.coeffs_)
print("R", m.react_sps)

print("K1: ", m.K1)
print("K2: ", m.K2)

print("A: ", m.A)
print("B: ", m.B)


#model = ODE(model=m, start=0, stop=100, epochs=1000)
#model = SSA(model=m, start=0, stop=100, epochs=1000)
model = TauLeaping(model=m, start=0, stop=100, max_epochs=100)

model.simulate()


a = list(model.species["A"])
b = list(model.species["B"])
time = list(model.species["Time"])
print(a)
print(b)
print(time)
print(len(a))
print(len(b))
print(len(time))

plt.plot(time, a, label="A")
plt.plot(time, b, label="B")
plt.legend()
plt.show()



"""
model1 = CLE(m)
print(model1.A)
print(model1.B)

model1.simulate()
print("this is it ", model1.species)
#print(model.species)
print(len(model1.species["Time"]))
#s= ODE(m)

#s.simulate()

#print(s.species)
#print(s.parameters)
#print(s.components)
#print(s.params)
"""



