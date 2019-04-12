# CIVL 598L Project

# Preliminaries

This paper is divided into Part 1 and Part 2. In Part 1, I examine the propagation of weakly nonlinear internal waves in a closed basin; in Part 2, I examine surface gravity waves in counter-propagating flows over an obstacle.

# Part 1: Internal waves in a two-layer system
# Introduction

Internal waves play an important role in the transport of heat, nutrients, and other scalars in lakes and reservoirs. Internal wave breaking leads to enhanced mixing between layers in a density stratified fluid. Predicting the mixing processes in a lake requires an understanding of the mechanisms controlling the generation, evolution, and degeneration of internal waves. Here, I focus on interfacial gravity waves in a two-layer system.

When the wind blows over a lake, a stress is applied to the water surface, which tilts the water surface in the downwind direction. To maintain a hydrostatic pressure gradient, the thermocline tilts in the opposite direction. Since the density difference between the two layers is typically small, the deflection of the thermocline is much larger than that of the water surface. 

After the wind stops, there is no longer an applied stress at the water surface to maintain the surface and interfacial tilts, and as a result, basin-scale seiches develop — one on the water surface and one on the interface. The amplitude of the water surface seiche is  typically small and will no longer be considered here. The amplitude of the interfacial seiche, however, can be large, reaching the height of the epilimnion.

Interfacial waves resulting from wind setup are often described using linear wave theory [@Mortimer:1952]; however, field observations have shown that wave amplitudes can be large enough so that nonlinear effects become important [@Farmer:1978]. These nonlinear effects can have important consequences on the physical processes in the water body. For example, nonlinear steepening can transfer wave energy from the basin scale to shorter length scales, generating packets of shorter waves. These shorter waves, or ‘solitons’, are prone to shoal at sloping boundaries, resulting in enhanced boundary mixing [@Horn:2002].

The goal in Part 1 is to describe the weakly nonlinear model from @Horn:2002 that was used to predict the evolution of interfacial waves in a closed basin. As a component of this project, I wrote a computer code to try to reproduce the results presented in @Horn:2002. The rest of the paper is organized as follows. In section 2, I present the theoretical formulation used here to study weakly nonlinear interfacial waves resulting from wind setup. In section 3, I describe the numerical method used to compute the results presented herein. In section 4, I show a set of numerical results starting with the simplest and gradually becoming more complex. Finally, in section 5, I discuss the results, summarize the key findings, and propose some future directions.

# Theoretical formulation





# Numerical method


# Results

Two sets of numerical results are presented. First, using the code developed in the present study, numerical results are compared to the classical numerical experiments from @Zabusky:1965. Then, after adapting the code to simulate interfacial waves,  I reproduce numerical and laboratory experiments from @Horn:2002.   


## Solitons

To verify that the numerical model developed in the present study reproduces the expected results, I compare them to those presented in @Zabusky:1965. They solve the KdV equation of the form:

$$u_t + u u_x + \delta^2 u_{xxx} = 0$$,{#eq:zabusky}

with initial conditions $u(x, 0) = \cos(\pi x)$, periodic boundary conditions on the interval $0 \leq x < 2$, and $\delta^2 = 0.022$. A side-by-side comparison of the computed results from @Zabusky:1965  with those computed in the present study show good agreement (Fig. 1). The good agreement provides some confidence that the numerical model in the present study working as intended. More rigorous validation would include comparing the computed results to an analytical solution and calculating the error (e.g. L2-norm error) for a set of spatial resolutions to determine if (a) the computed results converge to the known solution and (b) at what spatial resolution do the results become mesh independent.  


![Temporal development of the wave $u(x,t)$ a $t = 0, \; \pi^{-1}, ; \text{and} \; 3.6 \pi^{-1}$. (a) Numerical results from @Zabusky:1965; (b) Numerical results from the present study.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555045768139_zabusky.png)



## Nonlinear long internal waves in closed basins


![](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555046160421_horn.png)







# Summary


# Part 2: Surface gravity waves in inhomogeneous flows


# Introduction


# Theoretical formulation


![Problem definition sketch. A surface gravity wave propagates to the right with a phase speed $c_p(k)$ in a counter-propagating current with velocity $v(x)$. The presence of the obstacle gives rise to an inhomogeneous flow.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555046859980_schematic.png)



Consider 



# Numerical method


# Results




![Velocity potential plotted as a function of space and time for a packet of waves propagating against a current with Froude numbers $Fr$ as indicated on the respective panels. $Fr_{min}$ occurs far away from the obstacle and $Fr_{max}$ occurs at $x=0.75$, at the crest of the obstacle. $Fr = Fr_{min}$ for $x < 0.65$ and $x > 0.85$ and $Fr_{min} \leq Fr \leq Fr_{max}$ for $0.65 \leq x \leq 0.85$. The vertical dashed lines indicate the extent of the obstacle.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555046851332_spacetime_jet.png)

# Summary


# Appendix




# References










