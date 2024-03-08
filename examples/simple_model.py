from biostoch.model import Model
from biostoch.ode import EulerSimulator, RungeKuttaSimulator
from biostoch.ssa import GillespieSimulator
from biostoch.tau_leaping import TauLeaping
from biostoch.cle import ChemicalLangevin
from biostoch.visualization import Visualization

obj = Model()
obj.parameters({"K1": 0.1, "K2": 0.05})

obj.species(components={"A": 100.0, "B": 0.0}, rate_change={"A": "K2 * B - K1 * A", "B": "K1 * A - K2 * B"})

obj.reactions(reacts={"reaction1": "A -> B", "reaction2": "B -> A"},
              rates={"reaction1": "K1 * A", "reaction2": "K2 * B"})


model1 = EulerSimulator(model=obj, start=0, stop=100, epochs=1000)
model2 = RungeKuttaSimulator(model=obj, start=0, stop=100, epochs=1000)
model3 = GillespieSimulator(model=obj, start=0, stop=100, max_epochs=1000)
model4 = TauLeaping(model=obj, start=0, stop=100, max_epochs=100)
model5 = ChemicalLangevin(model=obj, stop=100, max_epochs=1000)

model1.simulate()
model2.simulate()
model3.simulate()
model4.simulate()
model5.simulate()

model11 = Visualization(model1)
model12 = Visualization(model2)
model13 = Visualization(model3)
model14 = Visualization(model4)
model15 = Visualization(model5)

model11.plot()
model12.plot()
model13.plot()
model14.plot()
model15.plot()
