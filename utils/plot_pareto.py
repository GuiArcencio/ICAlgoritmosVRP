import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('ggplot')

L1 = pd.read_csv('L1_normalized.csv', header=None)
L1_nadir = pd.read_csv('L1_nadir_normalized.csv', header=None)
Linf = pd.read_csv('tchebycheff_normalized.csv', header=None)
Linf_nadir = pd.read_csv('tchebycheff_nadir_normalized.csv', header=None)

fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
fig.suptitle('c50 (3 layers)')
fig.tight_layout(pad=1.5)

axs[0,0].plot(L1[0], L1[1], 'ro')
axs[0,0].set_title('L1 normalized')

axs[0,1].plot(L1_nadir[0], L1_nadir[1], 'bo')
axs[0,1].set_title('L1 normalized with nadir')

axs[1,0].plot(Linf[0], Linf[1], 'go')
axs[1,0].set_title('L∞ normalized')

axs[1,1].plot(Linf_nadir[0], Linf_nadir[1], 'mo')
axs[1,1].set_title('L∞ normalized with nadir')

for ax in axs.flat:
    ax.set(xlabel='distance traveled', ylabel='carbon emission')

plt.show()