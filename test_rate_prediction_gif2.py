import siegert
import nest
import nest.voltage_trace
import numpy as np
import matplotlib.pyplot as plt
nest.ResetKernel()
nest.set_verbosity('M_WARNING')
try:
	nest.Install('gif2_module')
except:
	pass

# parameters
T = 10000.0
t_start = 1000.0
N = 10

nu = 55000.0
J = 0.125 * 87.8

# TODO: get I_0 from rate

# scale connection strength by tau_syn such that results for different
# synaptic time constants are comparable

neuron_params2 = {
	"tau_1":      100.0,
	"C_m":        500.0,
	"tau_syn_ex":   0.5,
	"tau_syn_in":   0.5,
	"g_rr":       25.0,
	"g":          25.0,
#	"V_m":        -70.0,
	"V_reset":    -56.0,
	"E_L":        -56.0,
	"V_th":       -50.0,
	"t_ref":        2.0}

tau_syn = neuron_params2[ 'tau_syn_ex' ]
PSC = J / tau_syn

nest.SetDefaults('gif2_psc_exp', neuron_params2)
neuron = nest.Create("gif2_psc_exp", N)

# connect poisson input
syn_dict = {'weight': PSC, 'delay': 1.5}
noise = nest.Create("poisson_generator")
nest.SetStatus(noise, [ {"rate": nu} ])
nest.Connect(noise, neuron, 'all_to_all', syn_spec=syn_dict)

# connect spike detectors
nest.SetDefaults('spike_detector', {'start': t_start})
spikedetector = nest.Create("spike_detector")
nest.Connect(neuron, spikedetector, 'all_to_all')

nest.Simulate(T)

# get firing rate of neurons
rate = nest.GetStatus(spikedetector, "n_events")[ 0 ]
print 'frequency neuron', 1000 * rate / (T - t_start) * 1 / float(N)

PSC *= 1e-3
tau_syn *= 1e-3
taum = 1e-3 * neuron_params2[ 'C_m' ] / neuron_params2[ 'g' ]
taur = 1e-3 * neuron_params2[ 't_ref' ]
C = 1e-6 * neuron_params2[ 'C_m' ]
V_th = neuron_params2[ 'V_th' ]
V_reset = neuron_params2[ 'V_reset' ]
gamma = neuron_params2[ 'g_rr' ] / neuron_params2[ 'g' ]

# TODO: is this correct???
sigma = PSC * tau_syn / C * np.sqrt(taum * nu)

mu = neuron_params2[ 'V_reset' ] + PSC * tau_syn / C * taum * nu
mu2 = siegert.mu_(J*nu/tau_syn, neuron_params2[ 'g_rr' ], V_th, V_reset, C, np.sqrt(sigma))
# mu = V0 + PSC * tau_syn / C * taum * nu
# sigma = PSC * tau_syn / C * np.sqrt(taum * nu)

# sigma = I_N / g * np.sqrt(tau_N / tau_m)  # (A9)
# => I_0 * sqrt(tau_N) = J * sqrt(nu)

J_0 = PSC / neuron_params2[ 'g' ] * nu

# standard siegert
rate_gif2 = siegert.nu0_gif2(taum, tau_syn, taur, V_th, V_reset, mu2, sigma, gamma, 0)
print 'rate gif2 siegert', rate_gif2
rate_th = siegert.nu_0(taum, taur, V_th, V_reset, mu, sigma)
print 'rate siegert', rate_th
# siegert for filtered synapses
rate_th = siegert.nu0_fb(taum, tau_syn, taur, V_th, V_reset, mu, sigma)
print 'rate siegert filtered synapse', rate_th
rate_th = siegert.nu0_fb433(taum, tau_syn, taur, V_th, V_reset, mu, sigma)
print 'rate siegert filtered synapse', rate_th

'''
Results:

tau_syn = 0.001
frequency neuron 17.1885858586
rate siegert 17.9756406223
rate siegert filtered synapse 17.7219741779
rate siegert filtered synapse 17.7213200161

tau_syn = 0.1
frequency neuron 15.9861616162
rate siegert 17.9756406223
rate siegert filtered synapse 15.5011012078
rate siegert filtered synapse 15.4324345605

tau_syn = 1.0
frequency neuron 11.3628282828
rate siegert 17.9756406223
rate siegert filtered synapse 10.6977601049
rate siegert filtered synapse 9.93331690802

'''
