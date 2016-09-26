import matplotlib.pyplot as plt

def figcplots(gains_low, gains_high, frequencies, exp_f_r, exp_r_0, C_m, g, g_1, tau_1, save=True):
    fig = plt.figure("Noise-dependent gain of the GIF2")
    plt.subplot(211)
    plt.semilogx(frequencies[1:], gains_high[ 1:, 0 ], color='blue', marker='o', ls=' ')
    plt.semilogx(frequencies[1:], gains_low[ 1:, 0 ], color='red', marker='o', ls=' ')
    plt.grid(True)
    plt.text(0.25, 3.3, 'C = {0}pF, g = {1}nS, g1 = {2}nS, t1 = {3}ms'.format(C_m, g, g_1, tau_1))
    plt.axvline(x=exp_f_r*1000, color='green')
    plt.text(exp_f_r*1000 - 1., 0.25, r"$f_R$", color='green')
    plt.axvline(x=exp_r_0, color='green')
    plt.text(exp_r_0 - 2., 0.25, r"$r_0$", color='green')
    plt.xlabel('Frequency [Hz]', size=9.)
    plt.ylabel('Normalised amplitude of signal gain', size=9.)
    plt.xlim([ 0.09, 100. ])
    plt.ylim([0, 3.5])
    plt.subplot(212)
    plt.grid(True)
    plt.semilogx(frequencies[1:], gains_high[ 1:, 1 ], color='blue', marker='o', ls=' ')
    plt.semilogx(frequencies[1:], gains_low[ 1:, 1 ], color='red', marker='o', ls=' ')
    plt.xlabel('Frequency [Hz]', size=9.)
    plt.ylabel('Phase of signal gain', size=9.)
    plt.ylim([ -10, 10 ])
    plt.axvline(x=exp_f_r*1000)
    plt.axvline(x=exp_r_0)
    if save==True:
    	plt.savefig('Richardson_Fig6CD.png')
    else:
    	plt.ion()
    	plt.show()
    return fig