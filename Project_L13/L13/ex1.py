import numpy as np
import matplotlib.pyplot as plt

def predict_steps(A, x0, steps=5):
    #x_{k+1} = A*x_k
    A = np.asarray(A, dtype=float)
    x = np.asarray(x0, dtype=float).reshape(-1, 1)

    n, m = A.shape
    if n != m:
        raise ValueError(f"Matrix A must be square. Got {A.shape}")
    if x.shape[0] != n:
        raise ValueError(f"Vector x0 size must match A. Got x0={x.shape[0]}, A={n}")

    states = np.zeros((steps + 1, n), dtype=float)
    states[0] = x.flatten()

    for k in range(steps):
        x = A @ x
        states[k + 1] = x.flatten()

    return states

def plot_states(states):
    T, n = states.shape
    steps = np.arange(T)

    plt.figure()
    for i in range(n):
        plt.plot(steps, states[:, i], marker="o", label=f"x[{i}]")

    plt.title("Discrete Prediction Over Steps")
    plt.xlabel("Step k")
    plt.ylabel("Value")
    plt.xticks(steps)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    A = [
        [0.9, 0.1],
        [0.2, 0.8]
    ]
    x0 = [100, 50]

    states = predict_steps(A, x0, steps=5)

    for k, xk in enumerate(states):
        print(f"x_{k} = {xk}")

    plot_states(states)
