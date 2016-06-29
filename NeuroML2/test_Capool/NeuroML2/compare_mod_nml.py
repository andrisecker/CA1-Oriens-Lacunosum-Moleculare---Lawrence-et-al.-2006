import matplotlib.pyplot as plt


print("Comparing behaviour of Capool...")

time = []
mod_v = []
mod_cai = []
mod_ca1 = []
mod_ca2 = []
mod_ca3 = []
nml_v = []
nml_cai = []
nml_ca1 = []
nml_ca2 = []
nml_ca3 = []

orig_data = "../NEURON/nrn_ca.dat"
nml_data = "lems_ca.dat"

for line in open(orig_data):
    time.append(float(line.split()[0]))
    mod_v.append(float(line.split()[1]))
    mod_cai.append(float(line.split()[3]))
    mod_ca1.append(float(line.split()[4]))
    mod_ca2.append(float(line.split()[5]))
    mod_ca3.append(float(line.split()[6]))

for line in open(nml_data):
    nml_v.append(float(line.split()[1])*1000)
    nml_cai.append(float(line.split()[3]))
    nml_ca1.append(float(line.split()[4]))
    nml_ca2.append(float(line.split()[5]))
    nml_ca3.append(float(line.split()[6]))

fig = plt.figure(figsize=(15, 12))

ax = fig.add_subplot(2, 2, 1)

ax.plot(time, mod_cai, 'r-', label="NEURON")
ax.plot(time[:-1], nml_cai, 'b-', label="jnml")
ax.set_xlim([time[0], time[-1]])
ax.set_xlabel('time [ms]')
ax.set_ylabel('cai [mM]')

plt.legend()
plt.grid('on')

ax2 = fig.add_subplot(2, 2, 2)

ax2.plot(time, mod_ca1, 'r-', label="NEURON")
ax2.plot(time[:-1], nml_ca1, 'b-', label="jnml")
ax2.set_xlim([time[0], time[-1]])
ax2.set_xlabel('time [ms]')
ax2.set_ylabel('ca1 [mM]')

plt.legend()
plt.grid('on')

ax3 = fig.add_subplot(2, 2, 3)

ax3.plot(time, mod_ca2, 'r-', label="NEURON")
ax3.plot(time[:-1], nml_ca2, 'b-', label="jnml")
ax3.set_xlim([time[0], time[-1]])
ax3.set_xlabel('time [ms]')
ax3.set_ylabel('ca2 [mM]')

plt.legend()
plt.grid('on')

ax4 = fig.add_subplot(2, 2, 4)

ax4.plot(time, mod_ca3, 'r-', label="NEURON")
ax4.plot(time[:-1], nml_ca3, 'b-', label="jnml")
ax4.set_xlim([time[0], time[-1]])
ax4.set_xlabel('time [ms]')
ax4.set_ylabel('ca3 [mM]')

plt.legend()
plt.grid('on')


fig2 = plt.figure(figsize=(10,8))

ax = fig2.add_subplot(1, 1, 1)

ax.plot(time, mod_v, 'r-', label="NEURON")
ax.plot(time[:-1], nml_v, 'b-', label="jnml")
ax.set_xlim([time[0], time[-1]])
ax.set_xlabel('time [ms]')
ax.set_ylabel('V [mM]')

plt.legend()
plt.grid('on')

plt.show()
