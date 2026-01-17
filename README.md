# README: Generic Linear Kalman Filter Implementation

This project provides a robust C++ implementation of a **Discrete-Time Linear Kalman Filter**. While specifically tested using the state-space model for **Position (), Velocity (), and Acceleration ()** from the paper *"Filtering of IMU Data using Kalman Filter"* by Naveen Prabu Palanisamy, the class architecture is entirely generic and can be adapted to any linear dynamics problem.

---

### 🚀 Implementation Overview

The core of this project is a `KalmanFilter` class that manages state estimation through a recursive feedback loop of **Prediction** and **Correction**.

* **State Vector ():** Tracks distance, velocity, and acceleration ().
* **Measurement ():** Processes raw accelerometer readings () to refine the entire state.
* **Architecture:** Uses raw pointer arithmetic for maximum performance on embedded hardware (e.g., Arduino, STM32) and memory-constrained environments.

---

### 🛠 Configuration & Parameters

The filter's behavior is defined by its matrices, which are initialized to follow physical kinematic laws:

1. **State Transition Matrix (A):** Implements .
2. **Measurement Matrix (H):** Mapped as `[0, 0, 1]` to ensure the accelerometer reading only updates the acceleration state.
3. **Process Noise (Q):** Represents uncertainty in the model (e.g., wind, vibration).
4. **Measurement Noise (R):** Represents the sensor's variance.

---

### ⚖️ Multi-Axis Handling (X, Y, Z)

To track motion in 3D space, this implementation recommends the **Parallel Filter** approach:

* **Approach:** Instantiate three separate `KalmanFilter` objects (one for each axis).
* **Pros:** * **Efficiency:** Significantly lower computational cost than a single 9x9 matrix.
* **Simplicity:** Leverages the orthogonality of X, Y, and Z axes.


* **Cons:** Does not account for cross-axis noise correlations (rarely necessary for standard IMU filtering).

---

### ⚠️ Current Limitations

* **Matrix Inversion:** The `invertMatrix` method currently supports up to **3x3** matrices using an analytical adjugate method. For higher-order systems (), a Gauss-Jordan elimination algorithm would be required.
* **Linearity:** This is a **Linear** Kalman Filter. It cannot handle non-linear relationships (e.g., orientation or circular motion) without being upgraded to an Extended Kalman Filter (EKF).

---

### 📈 Optimization & Real-Time Accuracy

To achieve high accuracy with real IMU hardware:

1. **Static Calibration:** Measure the sensor variance while stationary to set an accurate **R** value.
2. **Bias Removal:** Subtract constant gravity and sensor offsets before passing data to the `update()` function to prevent quadratic position drift.
3. **Q/R Tuning:** * **High Q / Low R:** Filter trusts the sensor more (faster response, more noise).
* **Low Q / High R:** Filter trusts the model more (smoother output, potential lag).