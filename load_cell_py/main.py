from PyLTSpice.LTSpiceBatch import LTSpiceBatch
import matplotlib.pyplot as plt

# Path to your LTspice executable
ltspice_path = "/Applications/LTspice.app/Contents/MacOS/LTspice"  # adjust for your OS

# Path to your .asc file
asc_file = "spice/andres_conditioning_circuit.asc"

# Create a batch simulator object
sim = LTSpiceBatch(ltspice_path)

# Run the simulation
raw_file = sim.run(asc_file)

# Parse the raw output
sim_data = sim.get_data(raw_file)

# Plot voltage of node "Vout" (replace with your node name)
plt.plot(sim_data['time'], sim_data['Vout'])
plt.xlabel("Time [s]")
plt.ylabel("Voltage [V]")
plt.title("LTspice Simulation: Vout")
plt.grid(True)
plt.show()
