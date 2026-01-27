import jax, optax, h5py
from jax import grad, jit
import jax.numpy as jnp
import jax.random as jrandom
from jax.scipy.linalg import cholesky

# function definitions: 
def model(L, gi): 
    L0 = L[:, :, 0]
    L1 = L[:, :, 1]
    L = L0 + gi * L1
    return jnp.dot(L.T, L)

def loss(L, gmts, truths):
    y = jnp.float32(0.0)  # Initialize y as a scalar with float32 type
    for i in range(gmts.shape[0]):
        gi = gmts[i]
        truth = truths[:,:, i]
        prediction = model(L, gi)
        y += jnp.mean((prediction - truth)**2)
    return y

@jit
def update(params, x, y, opt_state):
    grads = grad(loss)(params, x, y)
    updates, new_opt_state = optimizer.update(grads, opt_state)
    new_params = optax.apply_updates(params, updates)
    return new_params, new_opt_state

# Load data 
directory = 'PLEASE/SET/YOUR/SAVE/PATH/HERE/'
name = 'tas_covariances'
file_save = h5py.File(directory + name + '_model.hdf5', 'w')
for j in range(1, 13):
    file = h5py.File(directory + name + '.hdf5', 'r')
    xs = jnp.array(file['temperature'])
    ys = jnp.transpose(jnp.array(file['covariances ' + str(j)]),axes=(2, 1, 0)) # permute from julia
    file.close()
    # Initialize Arrays and Parameters
    N = ys.shape[0]
    Savg = jnp.mean(ys, axis = 2)
    L0guess = cholesky(Savg)
    Lguess = jnp.stack((L0guess, jnp.zeros((N,N))), axis=2)
    L = Lguess
    S = ys 
    gmts = xs
    # Set up Optimizer
    optimizer = optax.adam(0.001)
    opt_state = optimizer.init(L)
    # Optimize 
    Lstart = L.copy()
    for i in range(100000):  # Number of iterations
        L, opt_state = update(L, gmts, S, opt_state)
        if i % 100 == 0:
            current_loss = loss(L, gmts, S)
            print(f"Iteration {i}: Loss {current_loss}")
    # Save results
    Lpermuted = jnp.transpose(L, axes=(2, 1, 0)) # permute to julia
    file_save.create_dataset('L'+str(j), data=Lpermuted)

file_save.close()
