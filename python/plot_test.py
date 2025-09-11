import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Create data
x = np.linspace(0, 2*np.pi, 100)
y = np.linspace(0, 2*np.pi, 100)
X, Y = np.meshgrid(x, y)

# Six different patterns
patterns = [
    np.sin(X) * np.cos(Y),
    np.cos(X) * np.sin(Y),
    np.sin(X + Y),
    np.cos(X + Y),
    np.sin(X) + np.cos(Y),
    np.sin(X) - np.cos(Y)
]
cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo']

# Create figure
fig = plt.figure(figsize=(15, 3))

# Create GridSpec with custom spacing
# We'll make 7 columns: 3 plots | gap | 3 plots
gs = GridSpec(1, 7, width_ratios=[1, 1, 1, 0.2, 1, 1, 1],  # 0.2 creates the gap
             wspace=0, hspace=0,
             left=0.05, right=0.95,
             bottom=0.1, top=0.9)

# Plot positions: [0,1,2,4,5,6] (skip position 3 for gap)
plot_positions = [0, 1, 2, 4, 5, 6]

# Create subplots
for i, pos in enumerate(plot_positions):
    ax = fig.add_subplot(gs[pos])
    pc = ax.pcolormesh(X, Y, patterns[i], cmap=cmaps[i], shading='auto')
    
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Add label below each plot
    ax.text(0.5, -0.1, f'Plot {i+1}', 
            transform=ax.transAxes,
            ha='center', va='top')

# Add colorbar below all plots
cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.03])
cbar = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Value Scale')

plt.suptitle('Six Subplots with Gap Between 3rd and 4th', y=1.05)
plt.show()