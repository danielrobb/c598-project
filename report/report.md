---
title: CIVL 598L Project
author: Daniel Robb
geometry: margin=1.25in
---

\section*{Preliminaries}

This paper is divided into Part 1 and Part 2. In Part 1, I examine the propagation of weakly nonlinear internal waves in a closed basin; in Part 2, I examine surface gravity waves in counter-propagating flows over an obstacle.

\section*{Part 1: Internal waves in a two-layer system}

# Introduction

Internal waves play an important role in the transport of heat, nutrients, and other scalars in lakes and reservoirs. Internal wave breaking leads to enhanced mixing between layers in a density stratified fluid. Predicting the mixing processes in a lake requires an understanding of the mechanisms controlling the generation, evolution, and degeneration of internal waves. Here, I focus on interfacial gravity waves in a two-layer system.

When the wind blows over a lake, a stress is applied to the water surface, which tilts the water surface in the downwind direction. To maintain a hydrostatic pressure gradient, the thermocline tilts in the opposite direction. Since the density difference between the two layers is typically small, the deflection of the thermocline is much larger than that of the water surface. 

After the wind stops, there is no longer an applied stress at the water surface to maintain the surface and interfacial tilts, and as a result, basin-scale seiches develop — one on the water surface and one on the interface. The amplitude of the water surface seiche is  typically small and will no longer be considered here. The amplitude of the interfacial seiche, however, can be large, reaching the height of the epilimnion.

Interfacial waves resulting from wind setup are often described using linear wave theory [@Mortimer:1952]; however, field observations have shown that wave amplitudes can be large enough so that nonlinear effects become important [@Farmer:1978]. These nonlinear effects can have important consequences on the physical processes in the water body. For example, nonlinear steepening can transfer wave energy from the basin scale to shorter length scales, generating packets of shorter waves. These shorter waves, or *solitons*, are prone to shoal at sloping boundaries, resulting in enhanced boundary mixing [@Horn:2002].

The goal in Part 1 is to describe the weakly nonlinear model from @Horn:2002 that was used to predict the evolution of interfacial waves in a closed basin. As a component of this project, I wrote a computer code to try to reproduce the results presented in @Horn:2002. The rest of Part 1 is organized as follows. In section 2, I present the theoretical formulation used here to study weakly nonlinear interfacial waves resulting from wind setup. In section 3, I describe the numerical method used to compute the results presented herein. Finally, in section 4, I show a set of numerical results starting with the simplest and gradually becoming more complex. 

# Theoretical formulation

The equations of motion used to compute the results in Part 1 are for long weakly nonlinear and weakly dispersive waves. I will first address the nonlinearity and then the dispersive nature of these waves.


## Long weakly nonlinear waves

Starting with the Saint-Venant equations in one dimension, conservation of mass and momentum are:

$$\frac{\partial h}{\partial t} + \frac{\partial (hu)}{\partial x} = 0,$$ {#eq:mass}

$$\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + g \frac{\partial h}{\partial x} = 0,$${#eq:momentum}

where $h(x, t)$ is the depth of water from the free surface to the bottom, $u(x, t)$ is velocity, $g$ is gravitational acceleration and $x$ and $t$ are coordinates of space and time, respectively. The depth of water $h = H + \eta(x, t)$, where $H$ is a constant depth. These equations describe nonlinear shallow flows with long (non-dispersive) waves.

Equations (@eq:mass) and (@eq:momentum) can be written in vector form:

$$\frac{\partial}{\partial t} \begin{pmatrix} u \\ h \end{pmatrix} + \begin{pmatrix} u & g \\ h & u \end{pmatrix} \frac{\partial}{\partial x} \begin{pmatrix} u \\ h \end{pmatrix}=0.$$ {#eq:swe}

We can transform the coordinate system to be moving with a right-moving characteristic curve $X_+ = x - c_+t$. Calculating the eigenvalues of the matrix in (@eq:swe), gives $c_{\pm} = u \pm \sqrt{gh}$, and taking the right-moving one, it can be shown the $u_{X_+} = \sqrt{(g/h)} h_{X_+}$ and that $r_{\pm} \equiv u \pm 2 \sqrt{gh}$ is constant along the right (left)-moving  characteristic curves. Following @Sutherland:2010 Chp. 4, we use $r_-$ as an initial condition for the rightward-propagating wave; therefore:

$$u = 2\sqrt{gh} - 2\sqrt{gH}$${#eq:ini}

along characteristic curves $X_+$. Substituting (@eq:ini) into (eq:mass)  and noting that $h = H + \eta$:

$$\frac{\partial \eta}{\partial t} + \left(3 \sqrt{g(H + \eta)} - 2 \sqrt{gH} \right)\frac{\partial \eta}{\partial x} = 0.$${#eq:nonlinear}

If $\eta / H \ll 1$, then the linear wave equation is recovered, where $c_0 = \sqrt{gH}$. If we consider $\eta$ to be non-negligible, then we can see that waves with $\eta > 0$ will move faster than those with $\eta < 0$. This means the wave will steepen as the crests catch up to the troughs and the wave may begin to break. 

Let’s consider $\eta/H$ to be a small parameter $\epsilon$, then 

$$\sqrt{g(H + \eta)} = c_0 \sqrt{1 + \epsilon} \approx c_0 \left( 1 + \frac{1}{2}\epsilon + ...\right),$${#eq:expansion}

and neglecting terms of $\mathcal{O}(\epsilon^2)$ or higher, (@eq:nonlinear) becomes

$$\eta_t + c_0 \eta_x + \frac{3}{2}\frac{c_0}{H}\eta \eta_x = 0.$${#eq:nonlinear2}

The nonlinear term in (@eq:nonlinear2) appears in the model equations that are solved in the present study. They account for weakly nonlinear effects of long waves.


## Weakly dispersive waves

Under a different set of assumptions as above, let’s consider linear dispersive waves with the dispersion relation 

$$\omega = \sqrt{gk \tanh{kH}},$${#eq:disp}

where $\omega$ is the frequency and $k$ is the wave number. For small values of $kH$,  $\tanh(kH) \approx kH$ and (@eq:disp) reduces to $\omega / k = \sqrt{gH} = c_0$, i.e., the shallow water wave speed, which does not depend on $k$, so these are non-dispersive waves. If we let $kH$ be a small parameter $\epsilon$ then:

$$\omega = \sqrt{gk \tanh{kH}} \approx \sqrt{gk} \left( \epsilon^{1/2} - \frac{1}{6} \epsilon^{5/2} + ... \right) = \sqrt{gk} \left( \sqrt{kH} - \frac{(kH)^{5/2}}{6}\right) =  c_0k - \frac{c_0 H^2}{6} k^3.$${#eq:dispapprox}

For linear waves, we look for solutions in the form:

$$\eta(x, t) = A \exp \left[i(kx - \omega t )\right].$${#eq:linear}

 From @Whitham:1974 Chp. 11, a linear equation with constant coefficients can be written as:

$$P \left( \frac{\partial}{\partial t}, \frac{\partial}{\partial x} \right) \eta = 0,$${#eq:poly}

where $P$ is a polynomial. Substituting the lefthand side of (@eq:linear) into (@eq:poly) gives:

$$P(-i\omega, ik)=0,$$  {#eq:generaldispersion}

which is the dispersion relation for (@eq:linear). Therefore, given the dispersion relation in (@eq:dispapprox), the corresponding linear equation is:

$$\eta_t + c_0 \eta_x + \frac{1}{6}c_0 H^2 \eta_{xxx} = 0.$${#eq:kdvlinear}



## The Korteweg-de Vries equation

For a weakly nonlinear, weakly dispersive wave where the nonlinear and dispersion terms or of the same order, we can use (@eq:nonlinear2)  and (@eq:kdvlinear) to give

$$\eta_t + c_0 \eta_x + \frac{3}{2}\frac{c_0}{H} \eta \eta_x + \frac{1}{6} c_0 H^2 \eta_{xxx} = 0.$${#eq:kdv}


Equation @eq:kdv is the Korteweg-de Vries (KdV) equation, which I solve numerically in section 4. If we neglect the third and fourth terms on the righthand side of (@eq:kdv), the linear non-dispersive wave equation is recovered. If the third and fourth terms are of the same order, then $a/H \sim (H/L)^2$, where $a$ is a measure of wave amplitude, $H$ is the vertical scale and $L$ is a  measure of the wave length. 

There are many variants of the KdV equation; here we focus on a variant in the form:

$$\eta_t + (c_0 + \alpha \eta) \eta_x + \beta \eta_{xxx} = 0.$$ {#eq:kdvgeneral}

For a surface gravity wave, $c_0 = \sqrt{gH}$, $\alpha = \frac{3}{2} \frac{c_0}{H}$, and $\beta = \frac{1}{6}c_0 H^2$ and $H$ is the water depth. 


## The Korteweg-de Vries equation for an interfacial gravity wave

For an interfacial gravity wave in a two-layer system where the density differences between layers are small, the KdV equation can have the same form as in (@eq:kdvgeneral) but the coefficients will differ. If we define a reduced gravity as:

$$g^\prime = g\frac{\rho_2 - \rho_1}{\bar{\rho}},$${#eq:gprime}

then 

$$c_0 = \sqrt{g^\prime \frac{h_1 h_2}{h_1 + h_2}}, \quad \alpha = \frac{3}{2} c_0 \frac{h_1 - h_2}{h_1 h_2}, \quad \text{and} \quad \beta = \frac{1}{6}c_0 h_1 h_2,$$ {#eq:coeffs}

as described in @Helfrich:2006.  In (@eq:gprime) and (@eq:coeffs), the subscripts 1 and 2 stand for the top and bottom layer, respectively.

In the next section, I will describe the numerical method used to solve (@eq:kdvgeneral).


# Numerical method

The numerical method in the present study follows @Fornberg:1978 and is summarized herein. Consider the KdV equation

$$u_t + u u_x + u_{xxx} = 0.$${#eq:fornberg} 

The variable $u(x, t)$ is transformed into Fourier space with respect to $x$ so that derivatives with respect to $x$ become algebraic, which reduces the partial differential equation to an ordinary differential equation with respect to $t$. The transformation between real space and wave number space, and vice versa, are performed using the Fourier transform:   

$$\hat{u}(k) = \mathcal{F} \left[ u(x) \right] = \int_{-\infty}^{\infty} u(x) \exp(-i kx) dx \approx \sum_{m=0}^{n-1} u(x) \exp \left( -i \frac{m k}{n} \right),$${#eq:fft}

and the inverse Fourier transform:

$$u(x) = \mathcal{F}^{-1} \left[ \hat{u}(k) \right] = \int_{-\infty}^{\infty} \hat{u}(k) \exp(i kx)dk \approx \frac{1}{n}\sum_{m=0}^{n-1} \hat{u}(k) \exp \left(i \frac{m k}{n} \right),$${#eq:ifft}
 
written here in both continuous and discrete forms. Equation @eq:fornberg can then be solved numerically by representing the temporal derivatives as finite differences, and with a leap-frog time stepping scheme, the discrete solution is given by:

$$\frac{u(x, t + \Delta t) - u(x, t - \Delta t)}{2 \Delta t} + u \mathcal{F}^{-1} \left[ik \mathcal{F}(u) \right] + \mathcal{F}^{-1} \left[ -ik^3 \mathcal{F}(u) \right] = 0,$${#eq:discrete}

where the $u = u(x, t)$ in the second and third terms. This numerical method above was used to compute the results presented in the next section.

# Results

Two sets of numerical results are presented. First, using the code developed in the present study, numerical results are compared to the classical numerical experiments from @Zabusky:1965. Then, after adapting the code to simulate interfacial waves,  I reproduce numerical and laboratory experiments from @Horn:2002.   


## Solitons

To verify that the numerical model developed in the present study reproduces the expected results, I compare them to those presented in @Zabusky:1965. They solve the KdV equation of the form:

$$u_t + u u_x + \delta^2 u_{xxx} = 0,$${#eq:zab}

with initial conditions $u(x, 0) = \cos(\pi x)$, periodic boundary conditions on the interval $0 \leq x < 2$, and $\delta^2 = 0.022$. Note, (@eq:zab) is equivalent to (@eq:kdvgeneral) for $c_0=0$, $\alpha = 1$, and $\beta = 0.022$.

A side-by-side comparison of the computed results from @Zabusky:1965  with those computed in the present study show good agreement (Fig. @fig:zabusky). The good agreement provides some confidence that the numerical model in the present study working as intended. More rigorous validation would include comparing the computed results to an analytical solution and calculating the error (e.g. L2-norm error) for increasingly fine spatial resolutions to determine (a) if the computed results converge to the known solution and (b) at what spatial resolution do the results become mesh independent.


![Temporal development of the wave $u(x,t)$ a $t = 0, \; \pi^{-1}, ; \text{and} \; 3.6 \pi^{-1}$. (a) Numerical results from @Zabusky:1965; (b) Numerical results from the present study.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555045768139_zabusky.png){#fig:zabusky}



## Weakly nonlinear long internal waves in closed basins

The numerical method applied in the previous section is now adapted for interfacial wave propagation in a closed basin. The equations of motion are nearly the same as in (@eq:kdvgeneral) with constant coefficients determined by (@eq:coeffs). A dissipative term is added to (@eq:kdvgeneral) to account for laminar boundary layers along solid surfaces (see @Horn:2002 section 2.2 for details).

The objective here is to simulate weakly nonlinear interfacial waves in a closed basin and to compare the predictions in the present study to laboratory and numerical results in @Horn:2002. The motivation is to understand two-layer wave propagation arising from an initial tilt of the interface --- an idealization of interfacial waves generated by wind setup in lakes.

To study weakly nonlinear internal wave propagation in a closed basin, @Horn:2002 conducted laboratory experiments in a long tank that could rotate about a horizontal axis so that the interface of a two-layer fluid could be initially tilted (Fig. @fig:lab). The experiment started at $t=0$ when the tank was suddenly returned to a horizontal position. The interfacial tilt resulted in internal wave propagation which was recorded by wave gauges at several locations (Fig. @fig:lab).


![From Horn et al. (2002). (a) Schematic diagram of the experimental set-up. The ultrasonic wavegauges were located at the positions marked A, B and C. (b, c) The tank and the density structure immediately before and after an experiment commences: (b) initially tilted tank, (c) initial condition with the tank horizontal and interface inclined.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555474583513_lab.png){#fig:lab}


In @Horn:2002, a numerical model was developed to study the internal wave propagation under a variety of conditions. They use a Fourier (pseudospectral) method like the one presented in section 3. Since the KdV equation is uni-directional and this physical problem has waves propagating to the right and to the left, the authors applied a method of reflection where they extend the physical bi-directional domain to form a uni-directional domain of twice the former’s length. Given the initial conditions and the method of reflection described in  @Horn:2002, the authors were able to simulate bi-directional wave propagation from solutions to the KdV equation. The same method was used to compute the results in the present study. 

Comparisons of results from @Horn:2002 to results from the present study show good agreement (Fig. @fig:gauge; Fig. @fig:snapshot). The results demonstrate the utility of the numerical model in the present study to simulate weakly nonlinear internal waves in a two-layer system.


![Comparison of observed and simulated interface displacements. Time series of interfacial displacements at wave gauges B ($x=3.0 \, \mathrm{m}$) and C ($x = 4.5 \, \mathrm{m}$). (a, b) from @Horn:2002; (c, d) present study; (a, c) wave gauge B; (b, d) wave gauge C. For this experiment $h_1 = 23.2 \, \mathrm{cm}$, $h_2 = 5.8 \, \mathrm{cm}$, $\theta = 0.5 ^\circ$, and $\Delta \rho = 20 \, \mathrm{kg \, m^{-3}}$.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555046160421_horn.png){#fig:gauge}

 
![Comparison of simulated interface displacements plotted at intervals of $T_i/4$ where $T_i$ is the basin seiche period. (top) results from @Horn:2002; (bottom) results from the present study. For this experiment $h_1 = 23.2 \, \mathrm{cm}$, $h_2 = 5.8 \, \mathrm{cm}$, $\theta = 0.5 ^\circ$, and $\Delta \rho = 20 \, \mathrm{kg \, m^{-3}}$.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555473665865_horn2.png){#fig:snapshot width=5in}



\section*{Part 2: Surface gravity waves in a steady flow over an obstacle}

# Introduction

The motivation of Part 2 was to apply the theory presented in *Low-frequency analogue Hawking radiation: The Korteweg–de Vries model* by @Coutant:2018 to the laboratory experiments described in @Weinfurtner:2011 and @Lawrence:2012. While some initial effort was put into linking the theory presented in these studies, given the time constraints for the CIVL 598L project this idea was abandoned. Instead, I will present a mathematical formulation and numerical method for studying wave-current interaction problems like those described in @Weinfurtner:2011 and @Lawrence:2012. The rest of Part 2 is organized as follows. In section 6, I present the theoretical formulation used here to study linear surface gravity waves in an inhomogeneous flow. In section 7, I describe the numerical method used to compute the results presented herein. In section 8, I show numerical results similar to those presented in @Unruh:1995. Finally, in section 9, I propose some future directions.

# Theoretical formulation

Consider an incompressible, inviscid, irrotational fluid in two dimensions in the $xz$-plane (Fig. @fig:prob). There is a no penetration condition along the bottom boundary and a free-surface condition on the top boundary. A steady current flows to the left and a counter-propagating linear surface gravity wave travels to the right. An obstacle gives rise to a water depth $h(x)$ and velocity $v(x)$ that vary over the obstacle. Since the wave properties,( e.g. frequency $\omega$, wave number $k$, phase speed $c_p$ , group speed $c_g$) depend on $h$ and $v$, they too vary over the obstacle.


![Problem definition sketch. A surface gravity wave propagates to the right with a phase speed $c_p(k)$ in a counter-propagating current with velocity $v(x)$. The presence of the obstacle gives rise to an inhomogeneous flow. ](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555046859980_schematic.png){#fig:prob}



In a stationary wave field, the dispersion relation for a linear surface gravity wave is:

$$\omega^2 = gk\tanh{kh}$${#eq:dispstill}  

but in the presence of a current, the wave is doppler shifted, so (@eq:dispstill) becomes

$$\Omega^2 = (\omega - vk)^2 = gk\tanh{kh},$${#eq:dispmove}

where $\Omega(k)$ is the frequency as seen from an observer travelling with the mean flow (i.e., the intrinsic frequency in the co-moving frame) and $\omega(k)$ is the frequency as seen by a stationary observer (i.e.,  the local frequency in the laboratory frame). 

The generated waves have wave lengths that are long compared to the depth of water, so they are shallow water waves and they are non-dispersive, but as they approach the obstacle, $v(x)$ becomes more negative. Since $\omega$ remains constant, $k$ must increase, resulting in shorter wavelengths. With shorter wavelengths compared to the water depth, dispersive effects become important and the waves are no longer shallow water waves. An approach that can account for these effects is presented in @Lawrence:2012 and @Robertson:2016.

Starting with equation (15) from @Lawrence:2012, the perturbation velocity potential $\phi$ caused by the waves on a background flow is

$$\frac{\partial^2 \phi}{\partial t^2} + \frac{\partial }{\partial t} \left(v(x) \frac{\partial \phi}{\partial x} \right) + \frac{\partial}{\partial x} \left( (v(x)^2 - c(k)^2) \frac{\partial \phi}{\partial x} \right) + \frac{\partial}{\partial x}\left(v(x) \frac{\partial \phi}{\partial t} \right) = 0,$${#eq:eq15}

which can be re-written to give

$$(\partial_t + \partial_x v(x))(\partial_t + v(x) \partial_x )\phi + f^2(i\partial_x) \phi = 0,$${#eq:rob}

where

$$f^2(k) = gk \tanh(kh) \approx gk (kh - \frac{1}{3}(kh)^3) = c_0^2 k^2 (1 - \frac{1}{3}(kh)^2),$${#eq:disss}

and

$$f^2(i\partial_x) \phi(x, t) = \mathcal{F}^{-1} \left[f^2(k) \hat{\phi}(k, t) \right]$$

and substituting the approximation on the righthand side of (@eq:disss) into (@eq:rob) gives

$$(\partial_t + \partial_x v(x))(\partial_t + v(x) \partial_x )\phi -  c_0^2 \partial^2_x\phi - \frac{1}{3}c_0^2 h^2 \partial^4_x \phi= 0.$${#eq:rob2}

In the next section, I will describe the numerical method used to solve (@eq:rob).


# Numerical method

The numerical method used to compute the results presented herein followed @Unruh:1995 and @Robertson:2011. To include dispersive effects resulting in the higher order spatial derivatives, they use a similar Fourier method as described in Part 1 of this paper.

# Results

I simulated a wave packet propagating against a current for a variety of background flow conditions. Simulation results are shown for six different background flows (Fig. @fig:blackhole). The location of the obstacle is indicated by the vertical dashed lines. The colour gradient shows the variation in perturbation velocity potential. Outside the vertical dashed lines, the background flow has a Froude number equal to $Fr_{min}$  and inside the dashed lines, the Froude number increases smoothly to $Fr_{max}$ at $x=0.75$.

The background velocity profile was

$$v(x) = \frac{1}{2}(Fr_{max}-Fr_{min}) \tanh\left[24\cos \left(2\pi(x-0.25) \right) + 23\right]-\frac{1}{2}(Fr_{max}+Fr_{min}),$${#eq:vb} 

which is the same shape as in @Unruh:1995. The initial wave velocity potential was

$$\phi(x, 0) = \frac{\sigma}{2} \exp \left(-\frac{(x-x_0)^2}{2 \sigma^2}  + 8i \frac{x-x_0}{\sigma} \right),$${#eq:phio}

where $\sigma = 0.05$, $x_0 = 0.3$, and $i = \sqrt{-1}$. Simulations were run with a spatial discretization of $\Delta x = 1/2048$ on a domain $x \in [0, 1)$ with periodic boundary conditions. 

From Fig. @fig:blackhole, the wave passes over the obstacle for $Fr_{max} \leq 0.5$ and  the wave is blocked for $Fr_{max} \geq 0.7$, where $Fr_{min} = Fr_{max} - 0.4$ in all cases. For $Fr_{max}=0.6$, the  incident wave is partially arrested such that shorter waves are reflected downstream to the left and longer waves are transmitted to the right over the obstacle.

Videos showing the variation in velocity potential with space and time are provided as supplementary material. They can be viewed at https://youtu.be/xKXcPcEFGmQ . 


![Velocity potential plotted as a function of space and time for a packet of waves propagating against a current with Froude numbers $Fr$ as indicated on the respective panels. $Fr_{min}$ occurs far away from the obstacle and $Fr_{max}$ occurs at $x=0.75$, at the crest of the obstacle. $Fr = Fr_{min}$ for $x < 0.65$ and $x > 0.85$ and $Fr_{min} \leq Fr \leq Fr_{max}$ for $0.65 \leq x \leq 0.85$. The vertical dashed lines indicate the extent of the obstacle.](https://paper-attachments.dropbox.com/s_43F38CB5250439B4F4818F7140C5D7A0A2B6D38F80AC8423474177A81FD505C5_1555046851332_spacetime_jet.png){#fig:blackhole}


# Future directions

The results presented herein could be extended in a number of ways. A few future directions are listed below:


- Simulate numerically the Vancouver experiment presented in @Weinfurtner:2011 and @Lawrence:2012
- Review @Robertson:2016 who discuss scattering of subcritical waves of obstacles of different shapes and they discuss the Vancouver experiment presented in @Weinfurtner:2011 and @Lawrence:2012
- Review the paper *Minimal analytical model for undular tidal bore profile; quantum and Hawking effect analogies* by @Berry:2018 who compares his analytical model to undular bore data from Chanson
- Predict the properties of the outgoing waves (e.g. amplitude, phase speed) for various background flows and incoming waves.
- Link the theory @Lawrence:1987 to the theory in @Lawrence:2012?
- Connect the theory presented in @Coutant:2018 on the KdV model to previous work in  @Weinfurtner:2011


\section*{Supplementary material}

Videos showing the variation in velocity potential for various background flows are provided here: https://youtu.be/xKXcPcEFGmQ . 

# References










