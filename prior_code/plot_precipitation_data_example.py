import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# Load the data from the CSV file
data = pd.read_csv('../input/all_precipitation_data_BCC-CSM2-MR.csv')

# Ensure 'date_projected' column is in datetime format
data['date_projected'] = pd.to_datetime(data['date_projected'])

# Filter data for the selected year
# Define the year as a parameter
year_to_plot = 2029  # You can change this to any year you want to plot
data_selected_year = data[(data['date_projected'].dt.year == year_to_plot)]

# Main plot
fig, ax = plt.subplots(figsize=(8.5, 6.2), dpi=250)
ax.plot(data['date_projected'], data['ssp585'], label='Daily Precipitation (SSP 585)', linewidth=0.5)

# Set x-axis ticks to show the start of each year
yearly_ticks = pd.date_range(start=data['date_projected'].min(), end=data['date_projected'].max(), freq='YS')
ax.set_xticks(yearly_ticks)
ax.set_xticklabels([d.strftime('%Y') for d in yearly_ticks], rotation=45, fontsize=13)  # Set fontsize for x-tick labels in large figure
ax.set_ylim(top=100)

# Set y-axis tick labels and label font size
ax.tick_params(axis='y', labelsize=13)  # Set fontsize for y-tick labels
ax.set_ylabel('Projected Daily Precipitation (mm)', fontsize=17)  # Set y-axis label font size

# Set x-axis label font size
ax.set_xlabel('Year', fontsize=17)  # Set x-axis label font size

# Title for the main plot
# ax.set_title(f'Projected Daily Precipitation Data (SSP 585) for {year_to_plot}', fontsize=14)

# Inset plot for the selected year
inset_ax = inset_axes(
    ax, 
    width=3,  # Width in inches
    height=1.5,  # Height in inches
    loc='upper right',  # Position relative to the main plot
    bbox_to_anchor=(0.8, 0.98),  # Fine-tune position (x, y)
    bbox_transform=ax.transAxes  # Transform relative to ax
)
inset_ax.plot(data_selected_year['date_projected'], data_selected_year['ssp585'], color='orange', linewidth=0.5)
inset_ax.set_ylim(top=62)

# Adjust x-axis labels in the inset (frequency every 2 months)
inset_ax.set_xticks(pd.date_range(start=f'{year_to_plot}-01-01', end=f'{year_to_plot}-12-31', freq='2MS'))
inset_ax.set_xticklabels([d.strftime('%b') for d in pd.date_range(start=f'{year_to_plot}-01-01', end=f'{year_to_plot}-12-31', freq='2MS')], fontsize=11, rotation=45)  # Set fontsize for inset figure

# Set y-axis tick labels for inset figure
inset_ax.tick_params(axis='y', labelsize=12)  # Set fontsize for y-tick labels in inset

# Add title inside the inset figure
inset_ax.text(0.5, 0.75, f'{year_to_plot} Daily Data', ha='center', va='bottom', fontsize=13, transform=inset_ax.transAxes)

plt.tight_layout()
plt.savefig('./precipitation_data_example.png', dpi=300)
plt.show()

