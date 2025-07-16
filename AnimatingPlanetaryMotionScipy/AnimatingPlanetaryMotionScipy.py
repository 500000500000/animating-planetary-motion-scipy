import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp # differential equation solver

# Model Parameters
a = 5 # semi-major axis
b = 3 # semi-minor axis
c = 4 # focal distance such that a^2 = b^2 + c^2
T = 3 # period
posTol = 1e-4 # for solve_ivp
velTol = 1e-2 # for solve_ivp
rtol=1e-6 # for solve_ivp
P = 4*T # period of precession
H = 2*math.pi*a*b/P # when doing precession
#H = 0 # set to zero when not doing precession

# Animation Parameters
totalTime = P # total time of animation in seconds
timeStep = 1/50; # time step in seconds

##########################################################

# The Differential Equation
def diff_eqn(t, state):
    # Unpack the state to position and velocity.
    X = state[0] 
    Y = state[1]
    Xprime = state[2]
    Yprime = state[3]
    
    # Calculate the acceleration from the position.
    r = math.sqrt(X*X + Y*Y);
    coeff_inv_square =  ((2*math.pi/T)**2) * (a/r)**3
    coeff_inv_cube = (   (1/r**4) * H * ( (4*math.pi*a*b/T)  + H )   )
    coeff = -(coeff_inv_square + coeff_inv_cube)
    Xprimeprime = coeff*X
    Yprimeprime = coeff*Y

    # The derivative of position is velocity.
    # The derivative of velocity is acceleration.
    dX = Xprime
    dY = Yprime
    dXprime = Xprimeprime
    dYprime = Yprimeprime

    # Pack the derivative and return.
    return [dX, dY, dXprime, dYprime]
    
# Set the initial conditions,
x0 = a+c
y0 = 0
xprime0 = 0
yprime0 = b * (2*math.pi/T)*(a /(a+c))
gprime0 = H/ ( (a+c)**2 )
X0 = x0
Y0 = y0
Xprime0 = 0
Yprime0 = yprime0 + x0 * gprime0

# Pack the initial conditions for the solver.
state0 = [X0, Y0, Xprime0, Yprime0]

# Find the frame count.
frameCount = int(totalTime/timeStep);

# MAIN CALCULATION 
tArray = np.linspace(0, totalTime, frameCount);
t_span = (0, totalTime)
atol = np.array([posTol, posTol, velTol, velTol])
solution = solve_ivp(
    diff_eqn,
    t_span,
    state0,
    t_eval=tArray,
    method="LSODA",
    rtol=rtol,
    atol=atol
)

# See how hard the solver is working.
numCallsToDiffEqn = solution.nfev

# Unpack the solution for the animation.
states = solution.y
arr_T = states.T # transpose
split_arrays = [arr_T[:, i:i+1] for i in range(4)]
xArray = split_arrays[0]
yArray = split_arrays[1]

# Create the figure.
fig, ax = plt.subplots()
ax.set_aspect('equal')
border = 0.2
ax.set_xlim(-c-a - border, c+a + border)
ax.set_ylim(-c-a - border , c+a + border)

# Declare the objects to animate.
line, = ax.plot([], [], 'k-', lw=2) # k- is black solid
dot, = ax.plot([], [], 'bo', markersize=6) # bo is blue dot
sun_dot, = ax.plot([0], [0], 'yo', markersize=12)  # yo is yellow dot

# Animatation functions.
def init():
    line.set_data([], [])
    dot.set_data([], [])
    return line, dot, sun_dot

def update(frame):
    line.set_data(xArray[:frame+1], yArray[:frame+1])
    dot.set_data(xArray[frame], yArray[frame])
    return line, dot, sun_dot

# Run the animation.
timeStepInMilliseconds = timeStep * 1000;
ani = FuncAnimation(fig, update, frameCount, init_func=init, blit=True, interval=timeStepInMilliseconds)

###########################################################
# Save the animation if desired.
#ani.save("animation.gif", writer='pillow', fps=1000 // timeStepInMilliseconds)

plt.show()
