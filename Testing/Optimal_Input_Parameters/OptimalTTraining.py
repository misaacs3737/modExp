import numpy as np
import tensorflow as tf
from tensorflow import keras
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

# Load Data
data = np.load('optimal_t_data.npy', allow_pickle = True)
X = data[:, :-1]  # Input features
y = data[:, -1]   # Output

print("First 5 rows of the generated data:")
print(data[:5])  # Displaying Sample Rows

# Splitting the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Creating a simple neural network model
model = keras.Sequential([
    keras.layers.Dense(64, activation='relu', input_shape=(X_train.shape[1],)),
    keras.layers.Dense(64, activation='relu'),
    keras.layers.Dense(1)  # Output layer with single neuron for regression
])

# Compiling the model
model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mae'])

# Training the model
history = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2)

# Evaluating the model on test data
loss, mae = model.evaluate(X_test, y_test)
print("Test Loss:", loss)
print("Test MAE:", mae)

predictions = model.predict(X_test)

# Plotting actual vs. predicted values
plt.figure(figsize=(10, 6))
plt.plot(predictions, label='Predictions', color='blue', marker='o', linestyle='None', alpha=0.5)
plt.plot(y_test, label='Actuals', color='red', marker='x', linestyle='None', alpha=0.5)
plt.xlabel('Data Point Index')
plt.ylabel('t Value')
plt.title('Predictions vs Actuals')
plt.legend()
plt.grid(True)
plt.show()
