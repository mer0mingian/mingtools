from nested_dict import nested_dict
import csv
import time

# numofsim = read last simnum from file TODO
dd = nested_dict()
numofsim = 0
dd['simulation']['num'] = numofsim
# dd['simulation']['timestamp'] =
dd['input']['used'] = flag_stim_type
dd['input']['poisson']['rate'] = r_0
dd['input']['poisson']['amplitude'] = r_1
dd['input']['poisson']['frequency'] = f
#dd['input']['poisson']['phase'] = 0.0
#dd['input']['current']['mean'] = 0.0
dd['input']['current']['std'] = I_N
#dd['input']['current']['phase'] = 0.0
dd['input']['current']['frequency'] = f
dd['input']['current']['amplitude'] = I_1
dd['input']['current']['offset'] = I_0
dd['input']['stimlabel']= 'Richardson low' #options: 'Richardson high/low (equivalent)', 'manual', 'realistic'
dd['parameters']['neuron']['g'] = g
dd['parameters']['neuron']['g_1'] = g_1
dd['parameters']['neuron']['tau_1'] = tau_1
dd['parameters']['neuron']['C'] = C_m
dd['parameters']['model']['V_th'] = V_th
dd['parameters']['model']['V_reset'] = V_reset
dd['parameters']['model']['E_L'] = E_L
dd['parameters']['model']['t_ref'] = t_ref
dd['parameters']['synapse']['synweight'] = synweight
dd['parameters']['simulation']['N'] = N
dd['parameters']['simulation']['dt'] = dt
dd['parameters']['simulation']['nbins'] = nbins
dd['parameters']['simulation']['t_sim'] = t_sim
dd['parameters']['simulation']['t_run'] = t_run
dd['predictions']['alpha'] = alpha
dd['predictions']['beta'] = beta
dd['predictions']['STR resonance occurs?'] = False # beta >? ...
dd['predictions']['f_R'] = exp_f_r * 1000.
dd['output']['r_0'] = r_0_reconstructed
dd['output']['r_1'] = r_1_reconstructed
dd['output']['gain'] = gain


# >>> from nested_dict import nested_dict
# >>> nd= nested_dict()
# >>> nd["one"] = 1
# >>> nd[1]["two"] = "1 / 2"
# #
# #   convert nested_dict -> dict and pickle
# #
# >>> nd.to_dict()
# {1: {'two': '1 / 2'}, 'one': 1}
# >>> import pickle
# >>> binary_representation = pickle.dumps(nd.to_dict())
#
# #
# #   convert dict -> nested_dict
# #
# >>> normal_dict = pickle.loads(binary_representation)
# >>> new_nd = nested_dict(normal_dict)
# >>> nd == new_nd

