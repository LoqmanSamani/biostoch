import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class Model(object):

    def parameters(self, params):
        """
        params: A dictionary containing all rate constants in the system.
                Each key represents a rate constant (e.g., K1, K2, ..., Kn),
                and each corresponding value represents the value of that rate constant.
        """
        setattr(self, "params", list(params.keys()))

        for key, val in params.items():
            setattr(self, key, val)

    def species(self, components):

        """
        components: A dictionary containing all species in the system.
                   Each key represents a specie (e.g. A, B, ...), and
                   each corresponding value represents the initial concentration.
        """
        setattr(self, "components", list(components.keys()))

        for key, val in components.items():
            setattr(self, key, val)

    def reactions(self, reacts, rates):
        """
        reactions: A dictionary containing reaction equations for each reaction.
                   Each key represents the reaction name, and each value represents the reaction equation.
                   exp:
                       {reaction1: "A + B -> C", ...}

        rates: A dictionary containing rate equations (or propensity functions) for each reaction.
               Each key represents the reaction name, and each value represents the rate equation of that reaction.
               exp:
                   {reaction1: "K1 * A * B", ...}
        """
        signs = ["+", "->", "*"]

        for reaction, formula in reacts.items():
            react = []
            react_eq = formula.replace('+', ' + ').replace('->', ' -> ').replace('*', ' * ')
            components = react_eq.split()

            for component in components:
                if component in self.components or component in self.params:
                    react.append('[self.' + component + ']')
                elif component in signs:
                    react.append(component)
                else:
                    try:
                        react.append(float(component))
                    except ValueError:
                        print(f"This component {component} is not valid!")

            reaction_equation = ' '.join(str(r) for r in react)
            setattr(self, reaction, reaction_equation)

        for reaction, rate in rates.items():
            rate_equation = []
            rate_eq = rate.replace('*', ' * ').replace('+', ' + ').replace('->', ' -> ')
            components = rate_eq.split()

            for component in components:
                if component in self.components or component in self.params:
                    rate_equation.append('[self.' + component + ']')
                elif component in signs:
                    rate_equation.append(component)
                else:
                    try:
                        rate_equation.append(float(component))
                    except ValueError:
                        print(f"This component {component} is not valid!")

            rate_equation1 = ' '.join(str(r) for r in rate_equation)
            setattr(self, reaction + "_rate", rate_equation1)


class SSA(Model):
    pass

