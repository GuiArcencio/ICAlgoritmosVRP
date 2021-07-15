import matplotlib.pyplot as plt

L1_nadir_norm = 0.2702545366493173
L1_norm = 0.3157423421527177
L1_not_norm = 0.3157423430794484

Linf_nadir_norm = 0.2750499210737993
Linf_norm = 0.2932826828420635
Linf_not_norm = 0.30214806010654094

cell_text = [
    [L1_not_norm, L1_norm, L1_nadir_norm],
    [Linf_not_norm, Linf_norm, Linf_nadir_norm]
]

plt.style.use('ggplot')
fig, axs = plt.subplots(1,1)
axs.axis('tight')
axs.axis('off')

table = axs.table(cellText=cell_text, rowLabels=["L1", "L∞"], colLabels=["Sem normalizar", "Normalização", "Normalização com nadir"], loc='center')

plt.show()