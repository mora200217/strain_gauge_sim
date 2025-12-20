import ltspice
import matplotlib.pyplot as plt
import numpy as np
import os

# Ruta a tu archivo .asc
asc_file = "spice/andres_conditioning_circuit.asc"

# Crear objeto LTspice
l = ltspice.Ltspice(asc_file)

# Simular (esto genera un archivo .raw automáticamente)
l.run()  

# Leer resultados
l.parse()  

# Lista de señales disponibles
print(l.get_trace_names())

# Ejemplo: graficar tensión 'V(out)'
t = l.get_time()
vout = l.get_data('V(out)')

plt.plot(t, vout)
plt.xlabel("Tiempo [s]")
plt.ylabel("Voltaje [V]")
plt.title("Simulación LTspice")
plt.grid(True)
plt.show()
