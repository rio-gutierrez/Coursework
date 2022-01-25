import pandas as pd
from matplotlib import pyplot as plt

# from jupyterthemes import jtplot
# jtplot.style()

simpson_data = pd.read_csv("./output", header = None)

dx    = simpson_data[0]
error = simpson_data[1]

fig, ax = plt.subplots(1, 1)

ax.loglog(dx, error, 'y-*', label='Error')
plt.legend()
ax.set_xlabel(r'$\Delta x$ (log)')
ax.set_ylabel('Error (log)')
plt.title('Log-Log error of Simpson\'s Rule')
plt.show()

# Save the figure
# plt.savefig('Error.pdf')
