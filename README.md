# CHAO: Complex Ginzburg-Landau & Spiral Defect Turbulence Engine

> **Live Interactive Application:** [https://hapybeing.github.io/foundry-ginzburg-landau/](https://hapybeing.github.io/foundry-ginzburg-landau/)

---

## 1. Abstract & Theoretical Overview

**CHAO** is an academic-grade, real-time interactive computational laboratory simulating the **Two-Dimensional Complex Ginzburg-Landau Equation (CGLE)**. Built using modern ES6+ and a high-performance **2D Strang Split-Step Fourier Spectral Solver**, the platform numerically integrates the universal amplitude equation governing non-equilibrium pattern formation, topological phase singularities, and spatiotemporal chaos.

The engine provides an interactive computational environment for exploring fundamental non-linear physical phenomena:
1. **The Benjamin-Feir-Newell (BFN) Instability:** Spontaneous breakdown of uniform limit-cycle oscillations when $1 + c_1 c_3 < 0$.
2. **Spiral Defect Turbulence:** Continuous spontaneous nucleation, chaotic migration, and topological annihilation of vortex-antivortex phase singularities where the order parameter amplitude vanishes ($|A| = 0$).
3. **Phase Turbulence vs. Defect Turbulence:** Transitions between smooth phase fluctuations with strictly positive amplitude ($|A| > 0$) and defect-mediated spatiotemporal chaos.
4. **Archimedean Spiral Wave Emission & Target Patterns:** Stable rotating spiral arms and autonomous pacemaker target rings.
5. **Exact Analytic Non-Linear Integration:** Utilizing the exact closed-form algebraic solution for the non-linear reaction step, eliminating truncation error and guaranteeing unconditional numerical boundedness.
6. **Topological Phase Winding Detection:** Real-time contour integration identifying singular defect cores with integer topological charges $n \in \{\pm 1\}$.

---

## 2. Mathematical & Algorithmic Foundations

### 2.1 The Universal Complex Ginzburg-Landau Equation

Near a continuous supercritical **Hopf bifurcation** in spatially extended continuous media, any oscillatory non-linear system can be asymptotically reduced via multiscale perturbation expansion to the **Complex Ginzburg-Landau Equation (CGLE)**:

$$\frac{\partial A(\vec{r}, t)}{\partial t} = \mu A + (1 + i c_1) \nabla^2 A - (1 - i c_3) |A|^2 A$$

where:
* $A(\vec{r}, t) \in \mathbb{C}$ is the macroscopic complex order parameter.
* $\mu > 0$ is the linear growth rate controlling the instability threshold.
* $c_1 \in \mathbb{R}$ represents linear dispersion (governing wave propagation and phase curvature).
* $c_3 \in \mathbb{R}$ represents non-linear frequency detuning (amplitude-dependent frequency shift).
* $\nabla^2 = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}$ is the 2D spatial Laplacian.

### 2.2 The Benjamin-Feir-Newell Instability Criterion

The CGLE admits a trivial spatially uniform, limit-cycle plane wave solution:

$$A_0(t) = \sqrt{\mu} \, e^{i \omega_0 t}, \quad \omega_0 = c_3 \mu$$

Applying linear stability analysis to perturbations $\delta A(\vec{r}, t) = u(\vec{r}, t) + i v(\vec{r}, t) \sim e^{\sigma t + i \vec{k}\cdot\vec{r}}$ yields the dispersion relation:

$$\sigma(k) = - (1 + c_1 c_3) k^2 - \frac{1}{4} c_3^2 (1 + c_1^2) \frac{k^4}{\mu} + \mathcal{O}(k^6)$$

For small wavenumbers $k \ll 1$, the growth rate $\sigma(k) > 0$ if and only if the coefficient of $k^2$ is negative. This yields the celebrated **Benjamin-Feir-Newell (BFN) Criterion**:

$$1 + c_1 c_3 < 0$$

* **BFN Stable Regime ($1 + c_1 c_3 > 0$):** Spatially uniform states and coherent Archimedean spiral waves are stable against long-wavelength phase perturbations.
* **BFN Unstable Regime ($1 + c_1 c_3 < 0$):** Phase fluctuations grow exponentially, resulting in the spontaneous generation of topological defects and spatiotemporal defect chaos.

### 2.3 Topological Phase Singularities & Spiral Defects

Representing the complex field in polar coordinates:

$$A(\vec{r}, t) = R(\vec{r}, t) \, e^{i \phi(\vec{r}, t)}$$

where $R = |A|$ is the real amplitude and $\phi = \text{arg}(A) \in [-\pi, \pi]$ is the phase.

At spatial locations where $R(\vec{r}_0) = 0$, the phase $\phi$ becomes mathematically undefined. These nodal points correspond to **topological defects** (spiral wave cores).

The topological charge (winding number) $n$ enclosed by any simple closed path $\Gamma$ is strictly quantized:

$$n = \frac{1}{2\pi} \oint_\Gamma \nabla \phi \cdot d\vec{r} \in \mathbb{Z}$$

In the defect turbulent regime:
* Singularities nucleate spontaneously in dipole pairs carrying opposite charges ($n = +1$ clockwise, $n = -1$ counter-clockwise), ensuring conservation of total topological charge:
  $$Q_{\text{net}} = \sum_i n_i = N_+ - N_- = 0$$
* The phase singularities propagate, undergo mutual scattering, and annihilate upon collision with opposite-charge anti-defects.

### 2.4 Exact Analytic Integration of the Non-Linear Reaction Step

In a decoupled operator-splitting framework, the local non-linear reaction equation is:

$$\frac{dA}{dt} = \mu A - (1 - i c_3) |A|^2 A$$

Expressing this in radial and angular coordinates:

$$\frac{dR}{dt} = \mu R - R^3 = R (\mu - R^2)$$

$$\frac{d\phi}{dt} = c_3 R^2$$

The radial equation is a Bernoulli differential equation whose **exact analytical solution** over a finite time step $\Delta t$ is:

$$\int_{R_0}^{R(t)} \frac{dR}{R(\mu - R^2)} = \int_0^t dt' \implies R(t) = \frac{R_0 e^{\mu t}}{\sqrt{1 + \frac{R_0^2}{\mu}(e^{2\mu t} - 1)}}$$

Integrating the phase differential $d\phi = c_3 R^2 dt$:

$$\Delta \phi(t) = c_3 \int_0^t R(t')^2 dt' = \frac{c_3}{2} \ln\left(1 + \frac{R_0^2}{\mu}(e^{2\mu t} - 1)\right)$$

**Computational Significance:** This exact algebraic formula completely eliminates non-linear numerical truncation error, preserves stability across arbitrary time steps, and strictly bounds the amplitude $R \le \sqrt{\mu}$.

### 2.5 Second-Order Strang Split-Step Spectral Algorithm

The linear spatial operator $\frac{\partial A}{\partial t} = (1 + i c_1) \nabla^2 A$ is solved in Fourier space:

$$\mathcal{F}\left[(1 + i c_1) \nabla^2 A\right] = - (1 + i c_1) (k_x^2 + k_y^2) \tilde{A}(\vec{k})$$

$$\tilde{A}(t + \Delta t, \vec{k}) = \exp\left( - (1 + i c_1) k^2 \Delta t \right) \tilde{A}(t, \vec{k})$$

Because $\text{Re}(-(1 + i c_1) k^2 \Delta t) = -k^2 \Delta t \le 0$, high-frequency modes are exponentially dissipated, suppressing grid ringing.

The complete second-order symmetric **Strang Splitting** step is:

$$A(t + \Delta t) = \hat{\mathcal{N}}_{\Delta t / 2} \circ \mathcal{F}^{-1} \left\{ e^{-(1 + i c_1) k^2 \Delta t} \mathcal{F}\left\{ \hat{\mathcal{N}}_{\Delta t / 2} \left[ A(t) \right] \right\} \right\} + \mathcal{O}(\Delta t^3)$$

where $\hat{\mathcal{N}}$ denotes the exact non-linear evaluation.

---

## 3. Interactive Parameter Controls & Guide

| Control | Parameter | Permissible Range | Physical Description |
| :--- | :--- | :--- | :--- |
| **Linear Dispersion** | $c_1$ | $-3.0 - +3.0$ | Phase dispersion coefficient. Modulates wave curvature and diffusion velocity. |
| **Non-Linear Detuning** | $c_3$ | $-3.0 - +3.0$ | Amplitude-dependent frequency detuning. Drives the Benjamin-Feir-Newell instability. |
| **Growth Rate** | $\mu$ | $0.1 - 2.0$ | Linear growth parameter determining asymptotic limit-cycle saturation $|A|_\infty = \sqrt{\mu}$. |
| **Time Step** | $\Delta t$ | $0.01 - 0.10$ | Integration delta. The exact non-linear formulation ensures stability across large $\Delta t$. |
| **Diagnostic Colormap** | Field | 4 Modes | Selects between Cyclic Phase Domain Coloring, Amplitude Magnitude, Phase Gradient Speed, and Real Component. |
| **Defect Core Glyphs** | $\pm 1$ | Boolean | Overlays real-time plaquette winding detectors highlighting vortices ($+1$, cyan) and anti-vortices ($-1$, magenta). |
| **Phase Razor** | Tool | User Drag | Slices the phase across a cursor drag line, nucleating pairs of phase singularities at the endpoints. |
| **Amplitude Pin** | Tool | User Drag | Forces local amplitude to zero ($R \to 0$), anchoring a stationary topological defect. |
| **Pacemaker** | Tool | User Drag | Locally accelerates frequency $c_3$, creating an autonomous target wave emission source. |

---

## 4. Curated Non-Linear Presets

1. **Spiral Defect Turbulence ($c_1 = -1.50, c_3 = 1.20$):** Deep within the BFN unstable regime ($1 + c_1 c_3 = -0.80 < 0$), exhibiting ongoing spontaneous vortex nucleation, chaotic motion, and annihilation.
2. **Rotating Archimedean Spiral ($c_1 = 0.20, c_3 = 0.80$):** BFN stable regime ($1 + c_1 c_3 = 1.16 > 0$), producing an ordered, stably rotating Archimedean spiral wave.
3. **Phase Turbulence ($c_1 = -2.00, c_3 = 0.70$):** Defect-free spatiotemporal chaos where the phase field fluctuates chaotically while amplitude remains strictly bounded away from zero ($|A| > 0$).
4. **Spatiotemporal Intermittency / Chimera ($c_1 = -1.20, c_3 = 1.50$):** Coexistence of coherent, laminar oscillatory domains alongside turbulent, defect-rich regions.
5. **Target Wave Pacemaker ($c_1 = 0.00, c_3 = 1.00$):** Expanding concentric circular wave trains emitted from a central pacemaker core.
6. **Topological Bound Dipole Orbit ($c_1 = 0.50, c_3 = 1.00$):** Counter-circulating vortex-antivortex pair orbiting each other in a stable bound state.

---

## 5. Architectural & Software Design

```
foundry-ginzburg-landau/
├── index.html              # HTML5 responsive canvas viewport, controls drawer, and KaTeX theory modal
├── style.css               # Glassmorphic dark theme, custom sliders, responsive HUD
├── app.js                  # FastFFT2D, CGLESolver, DefectTracker, CGLEAudio, and CGLERenderer
└── README.md               # Exhaustive academic documentation with LaTeX theoretical derivations
```

### Computational Pipeline Overview
1. **Field Allocation:** Real and imaginary components $A_{\text{re}}, A_{\text{im}}$ stored in contiguous `Float32Array(N * N)` buffers on a $128 \times 128$ spatial grid.
2. **Fast Fourier Transform:** In-place 2D Radix-2 Cooley-Tukey FFT with precomputed bit-reversal and trigonometric twiddle tables.
3. **Symmetric Strang Splitting:**
   * Exact analytic non-linear half-step $\Delta t / 2$.
   * Forward 2D FFT to momentum space.
   * Exact exponential dispersion $\tilde{A} \leftarrow \exp(-(1 + i c_1) k^2 \Delta t) \tilde{A}$.
   * Inverse 2D FFT to physical space.
   * Exact analytic non-linear half-step $\Delta t / 2$.
4. **Topological Defect Tracking:** Plaquette phase winding integration across all $2 \times 2$ grid cells locating nodes where $\oint \nabla \phi \cdot d\vec{r} = \pm 2\pi$.
5. **Acoustic Sonification:** Web Audio API synthesizer dynamically driven by the mean oscillation frequency $\bar{\omega} = -c_3 \langle |A|^2 \rangle$ and defect count.
6. **Rendering:** Direct Canvas 2D `ImageData` generation using cyclical HSV domain coloring.

---

## 6. References & Literature

1. **Ginzburg, V. L., & Landau, L. D. (1950).** On the theory of superconductivity. *Zh. Eksp. Teor. Fiz.*, 20, 1064.
2. **Aranson, I. S., & Kramer, L. (2002).** The world of the complex Ginzburg-Landau equation. *Reviews of Modern Physics*, 74(1), 99-143.
3. **Benjamin, T. B., & Feir, J. E. (1967).** The disintegration of wave trains on deep water Part 1. Theory. *Journal of Fluid Mechanics*, 27(3), 417-430.
4. **Newell, A. C. (1974).** Envelope equations. *Lectures in Applied Mathematics*, 15, 157.
5. **Cross, M. C., & Hohenberg, P. C. (1993).** Pattern formation outside of equilibrium. *Reviews of Modern Physics*, 65(3), 851-1111.
6. **Bohr, T., Jensen, M. H., Paladin, G., & Vulpiani, A. (1998).** *Dynamical Systems Approach to Turbulence*. Cambridge University Press.
7. **Strang, G. (1968).** On the construction and difference of difference schemes. *SIAM Journal on Numerical Analysis*, 5(3), 506-517.

---

## 7. License

MIT License. Designed, built, and shipped autonomously for non-linear dynamics research and exploratory algorithmic mathematics.
