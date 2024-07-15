import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# Generate sine wave data
num_cycles = 10
num_points = 200
t = np.linspace(0, num_cycles * 2 * np.pi, num_points)

y1 = np.sin(t)
y2 = 1.5 * np.sin(t + np.pi / 4)
y3 = np.sin(2 * t + np.pi / 2)

# Apply Hilbert transform to make the signals analytic
analytic_y1 = hilbert(y1)
analytic_y2 = hilbert(y2)
analytic_y3 = hilbert(y3)

# Stack the time series data
data = np.stack((analytic_y1, analytic_y2, analytic_y3))

# Window size for rolling calculation
window_size = 20

# Initialize lists to store the results
eigenvalues_list = []

# Perform the rolling window calculation
for i in range(num_points - window_size + 1):
    # Extract the window
    window_data = data[:, i:i+window_size]

    # Compute the covariance matrix for this window
    cov_matrix = np.cov(np.real(window_data))

    # Compute eigenvalues
    eigenvalues, _ = np.linalg.eig(cov_matrix)

    # Store results
    eigenvalues_list.append(eigenvalues)

# Convert list of eigenvalues to an array for easier plotting
eigenvalues_array = np.array(eigenvalues_list)

# Creating subplots
fig, axs = plt.subplots(2, 1, figsize=(12, 12))

# Plot the time series data
axs[0].plot(t, np.real(analytic_y1), label='Sine Wave 1')
axs[0].plot(t, np.real(analytic_y2), label='Sine Wave 2')
axs[0].plot(t, np.real(analytic_y3), label='Sine Wave 3')
axs[0].set_title('Analytic Time Series Data')
axs[0].set_xlabel('Time')
axs[0].set_ylabel('Amplitude')
axs[0].legend()
axs[0].grid(True)

# Plot the eigenvalues as a function of time
for i in range(eigenvalues_array.shape[1]):
    axs[1].plot(t[window_size-1:num_points], eigenvalues_array[:, i], label=f'Eigenvalue {i+1}')
axs[1].set_title('Eigenvalues Evolution Over Time')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('Eigenvalue')
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.show()