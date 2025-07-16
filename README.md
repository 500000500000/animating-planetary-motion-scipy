# Animating Planetary Motion Scipy

This is a companion project to animating-planetary-motion.

The function diff_eqn and array state0 are a differential equation and initial conditions
for a planet moving in an elliptical orbit with precession.

## Input Parameters

The input parameters are identical to those in animating-planetary-motion except that alphaTolerance
is replaced by tolerances on position and velocity.

a = semi-major axis

b = semi-minor axis

c = focal distance such that a^2 = b^2 + c^2

T = period

posTol = absolute tolerance for position

velTol = absolute tolerance for velocity

rtol = relative tolerance

P = period of precession

H = 2 * pi * a * b / P    use this formula, or set H = 0 to turn off precession

## Running The Program

Running the program calls solve_ivp with diff_eqn and state0 to produce an animation using parameters:

totalTime = time of animation in seconds

timeStep = time step in seconds

## Output

The animation may be saved as an animated gif, if desired.
