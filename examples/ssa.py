import matplotlib.pyplot as plt
from biostoch.model import Model
from biostoch.ssa import GillespieSimulator

obj = Model()

obj.parameters({"K1": 0.1, "K2": 0.05})

obj.species(components={"A": 100.0, "B": 0.0}, rate_change={"A": "K2 * B - K1 * A", "B": "K1 * A - K2 * B"})

obj.reactions(reacts={"reaction1": "A -> B", "reaction2": "B -> A"},
              rates={"reaction1": "K1 * A", "reaction2": "K2 * B"})

# Simulate the model 10 times
num_simulations = 10
simulations = []

for _ in range(num_simulations):
    model = GillespieSimulator(model=obj, start=0, stop=100, max_epochs=1000)
    model.reset()
    model.simulate()
    simulations.append(model.species)

# Plot all simulations on the same plot
plt.figure(figsize=(10, 6))
for i, sim in enumerate(simulations):
    plt.plot(sim["Time"], sim["A"])
    plt.plot(sim["Time"], sim["B"])

plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Stochastic Simulation Algorithm')

plt.show()
