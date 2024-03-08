class Model(object):

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


