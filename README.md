# JitterTime
A Matlab toolbox for calculating transient performance of control loops with timing jitter

JitterTime is a spin-off from the Matlab toolbox Jitterbug (Lincoln and Cervin, 2022) and can be used for calculating the performance of a controller under non-ideal timing conditions. Examples of such conditions include
delay and jitter due to CPU and network scheduling, lost samples or lost controls due to packet loss or execution overruns, and aperiodic behavior due to clock drift, asynchronous nodes, and random sampling. 
Both Jitterbug and JitterTime evaluate a quadratic cost function for a mixed continuous-time/discrete-time linear system driven by white noise and/or deterministic impulse disturbances. The main
difference is the timing model. In Jitterbug, the timing of the discrete systems are governed by random delays with specified probability density functions. This allows the total system to be treated as a jump-linear system, and covariance can be calculated by solving a set of linear equations. In JitterTime, however, the timing is arbitrary and completely driven by the user. This allows for more complex timing scenarios to be analyzed, including scheduling algorithms with long-term timing dependencies and asynchronous execution in distributed control systems. For deterministic timing scenarios over a finite horizon (or a repeating hyper\-period), the performance is evaluated exactly. For stochastic timing scenarios, however, lengthy Monte Carlo simulations can be needed to obtain results with high confidence. This is a major drawback of the tool compared to Jitterbug.

Using the Toolbox

JitterTime consists of a small number of functions and requires Matlab with the Control System Toolbox (any reasonably recent version). Simply add the directory containing the JitterTime functions to the Matlab path.
