import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Generate some data
mu, sigma = 12, 4  # Mean and standard deviation
x = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
y = norm.pdf(x, mu, sigma)

# Create the plot
plt.figure(figsize=(5, 3), dpi=250)
plt.plot(x, y, 'k')  # Black line for the probability distribution

# Shade the area under the curve
VaR = mu - 4
plt.fill_between(x, y, where=(x <= VaR), color='cornflowerblue', alpha=0.2)

plt.fill_between(x, y, where=(x >= VaR), color='lightgreen', alpha=0.2)

# Add the vertical dashed line
plt.axvline(x=VaR, color='red', linestyle='--')
plt.text(VaR + 0.4, 0.102, 'Resilience without uncertainty', color='red', fontsize=11, rotation=0)

# Add axis labels
plt.xlabel('Resilience', fontsize=12)
plt.ylabel('Probability', fontsize=12)

# Add "Value at Risk" text with an arrow pointing to the shaded area
plt.text(VaR - 7.25, 0.062, 'Value at Risk', fontsize=11, color='royalblue')
plt.annotate(
    '', xy=(VaR - 1, 0.01), xytext=(VaR - 4, 0.06),
    arrowprops=dict(arrowstyle='->', color='cornflowerblue', lw=1)
)

# Customize axes
ax = plt.gca()
ax.spines['top'].set_visible(False)  # Remove top border
ax.spines['right'].set_visible(False)  # Remove right border
ax.spines['left'].set_position('zero')  # Position y-axis at x=0
ax.spines['bottom'].set_position('zero')  # Position x-axis at y=0

# Remove ticks and tick labels
ax.set_xticks([])
ax.set_yticks([])

# Show the plot
plt.tight_layout()
plt.savefig(../output/resilience_VaR.png, dpi=250)
plt.show()
