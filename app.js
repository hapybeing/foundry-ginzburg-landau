/**
 * FOUNDRY // GINZBURG-LANDAU
 * 2D Complex Ginzburg-Landau Equation Spectral Solver & Spiral Defect Engine
 * 
 * Features:
 * - Second-order Strang Split-Step Fourier spectral algorithm
 * - Exact closed-form analytic solution for non-linear reaction step (zero truncation error)
 * - Exact Fourier-space linear dispersion exp(-(1 + i c_1) k^2 dt)
 * - Real-time topological phase singularity & defect tracking (±1 winding)
 * - Phase razor, amplitude pin, and pacemaker perturbation tools
 * - Web Audio API non-linear acoustic synthesizer driven by mean frequency and defect density
 */

// ============================================================================
// 1. FAST 2D FOURIER TRANSFORM (RADIX-2 COOLEY-TUKEY)
// ============================================================================

export class FastFFT2D {
    constructor(n) {
        this.n = n;
        this.bitRev = new Uint32Array(n);
        let j = 0;
        for (let i = 0; i < n - 1; i++) {
            this.bitRev[i] = j;
            let k = n >> 1;
            while (k <= j) { j -= k; k >>= 1; }
            j += k;
        }
        this.bitRev[n - 1] = n - 1;

        this.cosTable = [];
        this.sinTable = [];
        for (let len = 2; len <= n; len <<= 1) {
            const half = len >> 1;
            const ang = -2.0 * Math.PI / len;
            const c = new Float32Array(half);
            const s = new Float32Array(half);
            for (let k = 0; k < half; k++) {
                c[k] = Math.cos(ang * k);
                s[k] = Math.sin(ang * k);
            }
            this.cosTable.push(c);
            this.sinTable.push(s);
        }
    }

    fft1D(re, im, offset, stride, inv) {
        const n = this.n;
        const bitRev = this.bitRev;

        for (let i = 0; i < n; i++) {
            const rev = bitRev[i];
            if (i < rev) {
                const idxA = offset + i * stride;
                const idxB = offset + rev * stride;
                const tr = re[idxA]; re[idxA] = re[idxB]; re[idxB] = tr;
                const ti = im[idxA]; im[idxA] = im[idxB]; im[idxB] = ti;
            }
        }

        let stage = 0;
        for (let len = 2; len <= n; len <<= 1) {
            const half = len >> 1;
            const cTable = this.cosTable[stage];
            const sTable = this.sinTable[stage];
            stage++;

            for (let i = 0; i < n; i += len) {
                for (let k = 0; k < half; k++) {
                    const idxA = offset + (i + k) * stride;
                    const idxB = offset + (i + k + half) * stride;
                    const u_r = re[idxA], u_i = im[idxA];
                    const wr = cTable[k];
                    const wi = inv ? -sTable[k] : sTable[k];

                    const v_r = re[idxB] * wr - im[idxB] * wi;
                    const v_i = re[idxB] * wi + im[idxB] * wr;

                    re[idxA] = u_r + v_r;
                    im[idxA] = u_i + v_i;
                    re[idxB] = u_r - v_r;
                    im[idxB] = u_i - v_i;
                }
            }
        }

        if (inv) {
            const invN = 1.0 / n;
            for (let i = 0; i < n; i++) {
                const idx = offset + i * stride;
                re[idx] *= invN;
                im[idx] *= invN;
            }
        }
    }

    transform(re, im, inv = false) {
        const n = this.n;
        for (let y = 0; y < n; y++) {
            this.fft1D(re, im, y * n, 1, inv);
        }
        for (let x = 0; x < n; x++) {
            this.fft1D(re, im, x, n, inv);
        }
    }
}

// ============================================================================
// 2. COMPLEX GINZBURG-LANDAU EQUATION SPECTRAL SOLVER
// ============================================================================

export class CGLESolver {
    constructor(N = 128, L = 64.0) {
        this.N = N;
        this.L = L;
        this.dx = L / N;
        this.dt = 0.05;

        // CGLE Parameters
        this.c1 = -1.50; // Linear dispersion
        this.c3 = 1.20;  // Non-linear frequency detuning
        this.mu = 1.00;  // Linear growth rate

        this.time = 0.0;
        this.fft = new FastFFT2D(N);

        // Complex order parameter A(x, y) = A_re + i * A_im
        this.A_re = new Float32Array(N * N);
        this.A_im = new Float32Array(N * N);

        // Precomputed wavenumbers K^2 = kx^2 + ky^2
        this.k2 = new Float32Array(N * N);
        this.initWavenumbers();
        this.initPreset('defect_chaos');
    }

    initWavenumbers() {
        const N = this.N;
        const dK = (2.0 * Math.PI) / this.L;

        for (let y = 0; y < N; y++) {
            const kyVal = (y < N / 2 ? y : y - N) * dK;
            const ky2 = kyVal * kyVal;
            const rowOffset = y * N;
            for (let x = 0; x < N; x++) {
                const kxVal = (x < N / 2 ? x : x - N) * dK;
                this.k2[rowOffset + x] = kxVal * kxVal + ky2;
            }
        }
    }

    initPreset(preset = 'defect_chaos') {
        const N = this.N;
        const dx = this.dx;
        const halfL = this.L * 0.5;
        const total = N * N;

        if (preset === 'defect_chaos') {
            // Deep Benjamin-Feir-Newell unstable regime: 1 + c1*c3 = -0.8 < 0
            this.c1 = -1.50;
            this.c3 = 1.20;
            this.mu = 1.00;
            this.dt = 0.05;

            // Small initial random perturbations around unstable zero state
            for (let i = 0; i < total; i++) {
                this.A_re[i] = 0.15 * (Math.random() - 0.5);
                this.A_im[i] = 0.15 * (Math.random() - 0.5);
            }
        } else if (preset === 'archimedean') {
            // Stable single spiral wave: 1 + c1*c3 = 1.16 > 0
            this.c1 = 0.20;
            this.c3 = 0.80;
            this.mu = 1.00;
            this.dt = 0.05;

            for (let y = 0; y < N; y++) {
                const yCoord = y * dx - halfL;
                const row = y * N;
                for (let x = 0; x < N; x++) {
                    const xCoord = x * dx - halfL;
                    const r = Math.hypot(xCoord, yCoord);
                    const phi = Math.atan2(yCoord, xCoord);
                    const amp = Math.tanh(r * 0.8);
                    const phase = phi + 0.35 * r;
                    const idx = row + x;
                    this.A_re[idx] = amp * Math.cos(phase);
                    this.A_im[idx] = amp * Math.sin(phase);
                }
            }
        } else if (preset === 'phase_turb') {
            // Defect-free phase turbulence: c1 = -2.0, c3 = 0.7
            this.c1 = -2.00;
            this.c3 = 0.70;
            this.mu = 1.00;
            this.dt = 0.04;

            for (let i = 0; i < total; i++) {
                const phi = 0.4 * (Math.random() - 0.5);
                this.A_re[i] = 0.95 * Math.cos(phi);
                this.A_im[i] = 0.95 * Math.sin(phi);
            }
        } else if (preset === 'chimera') {
            // Spatiotemporal intermittency & chimera states: c1 = -1.2, c3 = 1.5
            this.c1 = -1.20;
            this.c3 = 1.50;
            this.mu = 1.00;
            this.dt = 0.04;

            for (let y = 0; y < N; y++) {
                const row = y * N;
                for (let x = 0; x < N; x++) {
                    const idx = row + x;
                    if (x < N * 0.5) {
                        // Ordered laminar domain
                        this.A_re[idx] = 0.95;
                        this.A_im[idx] = 0.0;
                    } else {
                        // Turbulent seed
                        this.A_re[idx] = 0.5 * (Math.random() - 0.5);
                        this.A_im[idx] = 0.5 * (Math.random() - 0.5);
                    }
                }
            }
        } else if (preset === 'target') {
            // Target Wave Pacemaker: c1 = 0.0, c3 = 1.0
            this.c1 = 0.0;
            this.c3 = 1.0;
            this.mu = 1.0;
            this.dt = 0.05;

            for (let y = 0; y < N; y++) {
                const yCoord = y * dx - halfL;
                const row = y * N;
                for (let x = 0; x < N; x++) {
                    const xCoord = x * dx - halfL;
                    const r = Math.hypot(xCoord, yCoord);
                    const idx = row + x;
                    const amp = 0.95;
                    const phase = Math.sin(r * 0.5);
                    this.A_re[idx] = amp * Math.cos(phase);
                    this.A_im[idx] = amp * Math.sin(phase);
                }
            }
        } else if (preset === 'dipole') {
            // Bound rotating dipole
            this.c1 = 0.50;
            this.c3 = 1.00;
            this.mu = 1.00;
            this.dt = 0.05;

            for (let y = 0; y < N; y++) {
                const yCoord = y * dx - halfL;
                const row = y * N;
                for (let x = 0; x < N; x++) {
                    const xCoord = x * dx - halfL;
                    const th1 = Math.atan2(yCoord - 3.5, xCoord);
                    const th2 = - Math.atan2(yCoord + 3.5, xCoord);
                    const phase = th1 + th2;
                    const r1 = Math.hypot(xCoord, yCoord - 3.5);
                    const r2 = Math.hypot(xCoord, yCoord + 3.5);
                    const amp = Math.tanh(r1 * 0.8) * Math.tanh(r2 * 0.8);
                    const idx = row + x;
                    this.A_re[idx] = amp * Math.cos(phase);
                    this.A_im[idx] = amp * Math.sin(phase);
                }
            }
        }

        this.time = 0.0;
    }

    /**
     * Exact Closed-Form Analytic Non-Linear Step
     * \frac{dA}{dt} = \mu A - (1 - i c_3) |A|^2 A
     * R(t) = \frac{R_0 e^{\mu t}}{\sqrt{1 + (R_0^2/\mu)(e^{2\mu t} - 1)}}
     * \Delta \phi = \frac{c_3}{2} \ln\left(1 + \frac{R_0^2}{\mu}(e^{2\mu t} - 1)\right)
     */
    applyNonlinearStep(stepDt) {
        const total = this.N * this.N;
        const mu = this.mu;
        const c3 = this.c3;

        const expMuT = Math.exp(mu * stepDt);
        const exp2MuT_minus_1 = Math.exp(2.0 * mu * stepDt) - 1.0;
        const invMu = 1.0 / Math.max(1e-6, mu);

        for (let i = 0; i < total; i++) {
            const re = this.A_re[i];
            const im = this.A_im[i];
            const R0_sq = re * re + im * im;

            if (R0_sq < 1e-18) {
                // If amplitude is tiny, simple linear growth
                this.A_re[i] *= expMuT;
                this.A_im[i] *= expMuT;
                continue;
            }

            const R0 = Math.sqrt(R0_sq);
            const phi0 = Math.atan2(im, re);

            const denom = 1.0 + R0_sq * invMu * exp2MuT_minus_1;
            const R1 = (R0 * expMuT) / Math.sqrt(Math.max(1e-12, denom));
            const dPhi = 0.5 * c3 * Math.log(Math.max(1e-12, denom));
            const phi1 = phi0 + dPhi;

            this.A_re[i] = R1 * Math.cos(phi1);
            this.A_im[i] = R1 * Math.sin(phi1);
        }
    }

    /**
     * Exact Fourier-Space Linear Dispersion Step
     * \tilde{A}(t + \Delta t) = \exp\left( - (1 + i c_1) k^2 \Delta t \right) \tilde{A}(t)
     */
    applyLinearDispersion(stepDt) {
        const total = this.N * this.N;
        const c1 = this.c1;

        // Forward 2D FFT
        this.fft.transform(this.A_re, this.A_im, false);

        for (let i = 0; i < total; i++) {
            const k2Val = this.k2[i];
            // Real damping factor exp(-k^2 * dt)
            const decay = Math.exp(- k2Val * stepDt);
            // Imaginary dispersion phase: exp(-i * c1 * k^2 * dt)
            const phase = - c1 * k2Val * stepDt;
            const cp = Math.cos(phase);
            const sp = Math.sin(phase);

            const u = this.A_re[i] * decay;
            const v = this.A_im[i] * decay;

            this.A_re[i] = u * cp - v * sp;
            this.A_im[i] = u * sp + v * cp;
        }

        // Inverse 2D FFT back to physical space
        this.fft.transform(this.A_re, this.A_im, true);
    }

    /**
     * Second-Order Strang Operator Splitting
     */
    step() {
        const halfDt = this.dt * 0.5;

        // 1. Non-linear half-step (exact analytic)
        this.applyNonlinearStep(halfDt);

        // 2. Linear dispersion full-step (Fourier space)
        this.applyLinearDispersion(this.dt);

        // 3. Non-linear second half-step (exact analytic)
        this.applyNonlinearStep(halfDt);

        this.time += this.dt;
    }

    /**
     * Interactive Tools
     */
    cutPhaseLine(startX, startY, endX, endY) {
        const N = this.N;
        const dx = this.dx;
        const halfL = this.L * 0.5;

        // Distance along segment
        const dist = Math.hypot(endX - startX, endY - startY);
        const steps = Math.max(10, Math.ceil(dist / dx * 2.0));

        for (let s = 0; s <= steps; s++) {
            const t = s / steps;
            const px = startX + t * (endX - startX);
            const py = startY + t * (endY - startY);

            const gx = Math.floor((px + halfL) / dx);
            const gy = Math.floor((py + halfL) / dx);

            if (gx >= 0 && gx < N && gy >= 0 && gy < N) {
                const idx = gy * N + gx;
                // Reverse phase by \pi
                this.A_re[idx] = -this.A_re[idx];
                this.A_im[idx] = -this.A_im[idx];
            }
        }
    }

    pinAmplitudeZero(cx, cy, radius) {
        const N = this.N;
        const dx = this.dx;
        const halfL = this.L * 0.5;
        const r2Max = radius * radius;

        for (let y = 0; y < N; y++) {
            const yCoord = y * dx - halfL;
            const row = y * N;
            for (let x = 0; x < N; x++) {
                const xCoord = x * dx - halfL;
                const d2 = (xCoord - cx) * (xCoord - cx) + (yCoord - cy) * (yCoord - cy);
                if (d2 < r2Max) {
                    const factor = d2 / r2Max;
                    const idx = row + x;
                    this.A_re[idx] *= factor;
                    this.A_im[idx] *= factor;
                }
            }
        }
    }

    injectNoise(magnitude = 0.25) {
        const total = this.N * this.N;
        for (let i = 0; i < total; i++) {
            this.A_re[i] += magnitude * (Math.random() - 0.5);
            this.A_im[i] += magnitude * (Math.random() - 0.5);
        }
    }

    /**
     * Compute Telemetry Diagnostics
     */
    computeDiagnostics() {
        const N = this.N;
        const total = N * N;
        let sumAmp = 0.0;
        let maxAmp = 0.0;

        for (let i = 0; i < total; i++) {
            const amp = Math.hypot(this.A_re[i], this.A_im[i]);
            sumAmp += amp;
            if (amp > maxAmp) maxAmp = amp;
        }

        const meanAmp = sumAmp / total;
        const bfnParameter = 1.0 + this.c1 * this.c3;
        const meanOmega = - this.c3 * (meanAmp * meanAmp);

        return {
            meanAmplitude: meanAmp,
            maxAmplitude: maxAmp,
            bfn: bfnParameter,
            meanOmega: meanOmega,
            time: this.time
        };
    }
}

// ============================================================================
// 3. TOPOLOGICAL PHASE SINGULARITY & DEFECT TRACKER
// ============================================================================

export class DefectTracker {
    static wrap(dth) {
        return (dth + Math.PI) % (2.0 * Math.PI) - Math.PI;
    }

    static findDefects(solver, maxDefects = 64) {
        const N = solver.N;
        const dx = solver.dx;
        const halfL = solver.L * 0.5;
        const re = solver.A_re;
        const im = solver.A_im;

        const defects = [];
        let netCharge = 0;

        for (let y = 1; y < N - 2; y++) {
            const row0 = y * N;
            const row1 = (y + 1) * N;

            for (let x = 1; x < N - 2; x++) {
                const i00 = row0 + x;
                const i10 = row0 + (x + 1);
                const i11 = row1 + (x + 1);
                const i01 = row1 + x;

                // Phase singularity occurs where amplitude is low
                const rhoMin = Math.min(
                    re[i00] * re[i00] + im[i00] * im[i00],
                    re[i10] * re[i10] + im[i10] * im[i10],
                    re[i11] * re[i11] + im[i11] * im[i11],
                    re[i01] * re[i01] + im[i01] * im[i01]
                );

                if (rhoMin > 0.45) continue;

                const th0 = Math.atan2(im[i00], re[i00]);
                const th1 = Math.atan2(im[i10], re[i10]);
                const th2 = Math.atan2(im[i11], re[i11]);
                const th3 = Math.atan2(im[i01], re[i01]);

                const d0 = this.wrap(th1 - th0);
                const d1 = this.wrap(th2 - th1);
                const d2 = this.wrap(th3 - th2);
                const d3 = this.wrap(th0 - th3);

                const winding = Math.round((d0 + d1 + d2 + d3) / (2.0 * Math.PI));

                if (winding !== 0) {
                    defects.push({
                        x: (x + 0.5) * dx - halfL,
                        y: (y + 0.5) * dx - halfL,
                        charge: winding
                    });
                    netCharge += winding;
                    if (defects.length >= maxDefects) break;
                }
            }
            if (defects.length >= maxDefects) break;
        }

        return {
            defects: defects,
            count: defects.length,
            netCharge: netCharge
        };
    }
}

// ============================================================================
// 4. NON-LINEAR ACOUSTIC SYNTHESIZER
// ============================================================================

export class CGLEAudio {
    constructor() {
        this.ctx = null;
        this.isPlaying = false;
        this.gainNode = null;
        this.oscBase = null;
        this.oscDetuned = null;
        this.filter = null;
    }

    init() {
        if (this.ctx) return;
        const AudioContext = window.AudioContext || window.webkitAudioContext;
        if (!AudioContext) return;

        this.ctx = new AudioContext();

        this.gainNode = this.ctx.createGain();
        this.gainNode.gain.setValueAtTime(0.0001, this.ctx.currentTime);
        this.gainNode.connect(this.ctx.destination);

        this.filter = this.ctx.createBiquadFilter();
        this.filter.type = 'lowpass';
        this.filter.frequency.setValueAtTime(350, this.ctx.currentTime);
        this.filter.Q.setValueAtTime(3.5, this.ctx.currentTime);
        this.filter.connect(this.gainNode);

        this.oscBase = this.ctx.createOscillator();
        this.oscBase.type = 'sine';
        this.oscBase.frequency.setValueAtTime(110.0, this.ctx.currentTime);

        this.oscDetuned = this.ctx.createOscillator();
        this.oscDetuned.type = 'triangle';
        this.oscDetuned.frequency.setValueAtTime(111.5, this.ctx.currentTime);

        this.oscBase.connect(this.filter);
        this.oscDetuned.connect(this.filter);

        this.oscBase.start();
        this.oscDetuned.start();
    }

    toggle() {
        if (!this.ctx) this.init();
        if (this.ctx.state === 'suspended') this.ctx.resume();

        this.isPlaying = !this.isPlaying;
        const targetGain = this.isPlaying ? 0.20 : 0.0001;
        this.gainNode.gain.cancelScheduledValues(this.ctx.currentTime);
        this.gainNode.gain.exponentialRampToValueAtTime(targetGain, this.ctx.currentTime + 0.5);
        return this.isPlaying;
    }

    update(meanOmega, defectCount) {
        if (!this.isPlaying || !this.ctx) return;

        const now = this.ctx.currentTime;
        const baseFreq = 80.0 + Math.abs(meanOmega) * 45.0;
        const beatFreq = 0.5 + defectCount * 0.4;

        this.oscBase.frequency.setTargetAtTime(baseFreq, now, 0.1);
        this.oscDetuned.frequency.setTargetAtTime(baseFreq + beatFreq, now, 0.1);

        const cutoff = Math.max(160, Math.min(1400, 240 + defectCount * 80));
        this.filter.frequency.setTargetAtTime(cutoff, now, 0.15);
    }
}

// ============================================================================
// 5. HIGH-PERFORMANCE 2D CANVAS RENDERER
// ============================================================================

export class CGLERenderer {
    constructor(canvas, solver, defectTracker) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.solver = solver;
        this.defectTracker = defectTracker;

        this.colormap = 'phase_domain'; // 'phase_domain', 'amplitude', 'phase_gradient', 'real_part'
        this.showDefects = true;
        this.showStreamlines = false;

        this.imgBuffer = this.ctx.createImageData(solver.N, solver.N);
        this.offscreenCanvas = document.createElement('canvas');
        this.offscreenCanvas.width = solver.N;
        this.offscreenCanvas.height = solver.N;
        this.offscreenCtx = this.offscreenCanvas.getContext('2d');
    }

    hsvToRgb(h, s, v) {
        let r, g, b;
        const i = Math.floor(h * 6);
        const f = h * 6 - i;
        const p = v * (1 - s);
        const q = v * (1 - f * s);
        const t = v * (1 - (1 - f) * s);

        switch (i % 6) {
            case 0: r = v; g = t; b = p; break;
            case 1: r = q; g = v; b = p; break;
            case 2: r = p; g = v; b = t; break;
            case 3: r = p; g = q; b = v; break;
            case 4: r = t; g = p; b = v; break;
            case 5: r = v; g = p; b = q; break;
        }
        return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
    }

    renderFields() {
        const N = this.solver.N;
        const total = N * N;
        const re = this.solver.A_re;
        const im = this.solver.A_im;
        const data = this.imgBuffer.data;

        if (this.colormap === 'phase_domain') {
            for (let i = 0; i < total; i++) {
                const u = re[i];
                const v = im[i];
                const amp = Math.min(1.2, Math.hypot(u, v));
                const phase = Math.atan2(v, u);
                const hue = (phase + Math.PI) / (2.0 * Math.PI);
                const [r, g, b] = this.hsvToRgb(hue, 0.90, amp * 0.95);

                const pIdx = i * 4;
                data[pIdx]     = r;
                data[pIdx + 1] = g;
                data[pIdx + 2] = b;
                data[pIdx + 3] = 255;
            }
        } else if (this.colormap === 'amplitude') {
            for (let i = 0; i < total; i++) {
                const amp = Math.min(1.2, Math.hypot(re[i], im[i]));
                const t = Math.min(1.0, amp / 1.0);
                const pIdx = i * 4;
                data[pIdx]     = Math.round(Math.pow(t, 0.6) * 255);
                data[pIdx + 1] = Math.round(Math.pow(t, 1.2) * 220);
                data[pIdx + 2] = Math.round(Math.pow(t, 2.2) * 255);
                data[pIdx + 3] = 255;
            }
        } else if (this.colormap === 'phase_gradient') {
            const inv2dx = 0.5 / this.solver.dx;
            for (let y = 1; y < N - 1; y++) {
                const row = y * N;
                for (let x = 1; x < N - 1; x++) {
                    const idx = row + x;
                    const u = re[idx], v = im[idx];
                    const rho = u * u + v * v + 1e-8;
                    const du_dx = (re[idx + 1] - re[idx - 1]) * inv2dx;
                    const dv_dx = (im[idx + 1] - im[idx - 1]) * inv2dx;
                    const du_dy = (re[idx + N] - re[idx - N]) * inv2dx;
                    const dv_dy = (im[idx + N] - im[idx - N]) * inv2dx;

                    const gradX = (u * dv_dx - v * du_dx) / rho;
                    const gradY = (u * dv_dy - v * du_dy) / rho;
                    const speed = Math.min(1.0, Math.hypot(gradX, gradY) * 0.4);

                    const pIdx = idx * 4;
                    data[pIdx]     = Math.round(speed * 20);
                    data[pIdx + 1] = Math.round(speed * 240);
                    data[pIdx + 2] = Math.round(speed * 255);
                    data[pIdx + 3] = 255;
                }
            }
        } else if (this.colormap === 'real_part') {
            for (let i = 0; i < total; i++) {
                const val = Math.max(-1.0, Math.min(1.0, re[i]));
                const pIdx = i * 4;
                if (val > 0) {
                    data[pIdx]     = Math.round(val * 246);
                    data[pIdx + 1] = Math.round(val * 211);
                    data[pIdx + 2] = Math.round(val * 101);
                } else {
                    data[pIdx]     = Math.round(-val * 255);
                    data[pIdx + 1] = 0;
                    data[pIdx + 2] = Math.round(-val * 127);
                }
                data[pIdx + 3] = 255;
            }
        }

        this.offscreenCtx.putImageData(this.imgBuffer, 0, 0);
    }

    render() {
        const ctx = this.ctx;
        const width = this.canvas.width;
        const height = this.canvas.height;
        const solver = this.solver;
        const halfL = solver.L * 0.5;

        ctx.fillStyle = '#04060a';
        ctx.fillRect(0, 0, width, height);

        this.renderFields();

        const size = Math.min(width, height) * 0.92;
        const startX = (width - size) * 0.5;
        const startY = (height - size) * 0.5;

        ctx.save();
        ctx.imageSmoothingEnabled = true;
        ctx.imageSmoothingQuality = 'high';
        ctx.drawImage(this.offscreenCanvas, startX, startY, size, size);

        const toPixX = (x) => startX + ((x + halfL) / solver.L) * size;
        const toPixY = (y) => startY + ((y + halfL) / solver.L) * size;

        // Draw Topological Phase Singularity Glyphs
        if (this.showDefects) {
            const detected = this.defectTracker.findDefects(solver);

            for (const d of detected.defects) {
                const px = toPixX(d.x);
                const py = toPixY(d.y);
                const isPos = d.charge > 0;

                ctx.shadowColor = isPos ? '#00f2fe' : '#ff007f';
                ctx.shadowBlur = 10;
                ctx.strokeStyle = isPos ? '#00f2fe' : '#ff007f';
                ctx.lineWidth = 2.0;

                ctx.beginPath();
                ctx.arc(px, py, 6.5, 0, Math.PI * 2);
                ctx.stroke();

                ctx.fillStyle = isPos ? '#00f2fe' : '#ff007f';
                ctx.font = 'bold 10px monospace';
                ctx.fillText(isPos ? '+1' : '-1', px + 8, py - 4);
                ctx.shadowBlur = 0;
            }
        }

        // Domain border
        ctx.strokeStyle = 'rgba(45, 65, 95, 0.6)';
        ctx.lineWidth = 1.0;
        ctx.strokeRect(startX, startY, size, size);

        ctx.restore();
    }
}

// ============================================================================
// 6. APPLICATION ENTRYPOINT & CONTROLLER
// ============================================================================

export class CGLEApp {
    constructor() {
        this.canvas = document.getElementById('cgle-canvas');
        this.solver = new CGLESolver(128, 64.0);
        this.defectTracker = DefectTracker;
        this.audio = new CGLEAudio();
        this.renderer = new CGLERenderer(this.canvas, this.solver, this.defectTracker);

        this.paused = false;
        this.toolMode = 'razor';
        this.toolRadius = 2.5;

        this.lastFpsUpdate = performance.now();
        this.fpsHistory = [];

        this.cacheDom();
        this.bindEvents();
        this.resize();
        this.loadPreset('defect_chaos');
    }

    cacheDom() {
        this.dom = {
            presetSelect: document.getElementById('preset-select'),
            btnAudio: document.getElementById('btn-audio'),
            btnTheory: document.getElementById('btn-theory'),
            theoryModal: document.getElementById('theory-modal'),
            btnCloseModal: document.getElementById('btn-close-modal'),
            btnToggleControls: document.getElementById('btn-toggle-controls'),
            btnToggleHud: document.getElementById('btn-toggle-hud'),
            controlsSidebar: document.getElementById('controls-sidebar'),
            telemetryHud: document.getElementById('telemetry-hud'),

            sliderC1: document.getElementById('slider-c1'),
            valC1: document.getElementById('val-c1'),
            sliderC3: document.getElementById('slider-c3'),
            valC3: document.getElementById('val-c3'),
            sliderMu: document.getElementById('slider-mu'),
            valMu: document.getElementById('val-mu'),
            sliderDt: document.getElementById('slider-dt'),
            valDt: document.getElementById('val-dt'),
            btnPlayPause: document.getElementById('btn-play-pause'),
            btnStep: document.getElementById('btn-step'),
            btnReset: document.getElementById('btn-reset'),
            btnInjectNoise: document.getElementById('btn-inject-noise'),

            selectColormap: document.getElementById('select-colormap'),
            checkDefects: document.getElementById('check-defects'),
            checkStreamlines: document.getElementById('check-streamlines'),

            selectTool: document.getElementById('select-tool'),
            sliderToolRadius: document.getElementById('slider-tool-radius'),
            valToolRadius: document.getElementById('val-tool-radius'),

            telBfn: document.getElementById('tel-bfn'),
            telDefects: document.getElementById('tel-defects'),
            telCharge: document.getElementById('tel-charge'),
            telMeanAmp: document.getElementById('tel-mean-amp'),
            telMaxAmp: document.getElementById('tel-max-amp'),
            telCoherence: document.getElementById('tel-coherence'),
            telOmega: document.getElementById('tel-omega'),
            telFps: document.getElementById('tel-fps'),
            telTime: document.getElementById('tel-time')
        };
    }

    bindEvents() {
        window.addEventListener('resize', () => this.resize());

        this.dom.presetSelect.addEventListener('change', (e) => this.loadPreset(e.target.value));

        this.dom.btnAudio.addEventListener('click', () => {
            const isPlaying = this.audio.toggle();
            this.dom.btnAudio.classList.toggle('active', isPlaying);
            this.dom.btnAudio.innerHTML = isPlaying 
                ? '<svg width="18" height="18" viewBox="0 0 24 24" fill="currentColor"><path d="M3 9v6h4l5 5V4L7 9H3zm13.5 3c0-1.77-1.02-3.29-2.5-4.03v8.05c1.48-.73 2.5-2.25 2.5-4.02zM14 3.23v2.06c2.89.86 5 3.54 5 6.71s-2.11 5.85-5 6.71v2.06c4.01-.91 7-4.49 7-8.77s-2.99-7.86-7-8.77z"/></svg> <span>AUDIO ON</span>'
                : '<svg width="18" height="18" viewBox="0 0 24 24" fill="currentColor"><path d="M16.5 12c0-1.77-1.02-3.29-2.5-4.03v2.21l2.45 2.45c.03-.2.05-.41.05-.63zm2.5 0c0 .94-.2 1.82-.54 2.64l1.51 1.51C20.63 14.91 21 13.5 21 12c0-4.28-2.99-7.86-7-8.77v2.06c2.89.86 5 3.54 5 6.71zM4.27 3L3 4.27 7.73 9H3v6h4l5 5v-6.73l4.25 4.25c-.67.52-1.42.93-2.25 1.18v2.06c1.38-.31 2.63-.95 3.69-1.81L19.73 21 21 19.73l-9-9L4.27 3zM12 4L9.91 6.09 12 8.18V4z"/></svg> <span>AUDIO OFF</span>';
        });

        this.dom.btnTheory.addEventListener('click', () => {
            this.dom.theoryModal.classList.add('visible');
            if (window.renderMathInElement) {
                try {
                    window.renderMathInElement(this.dom.theoryModal, {
                        delimiters: [
                            {left: '$$', right: '$$', display: true},
                            {left: '$', right: '$', display: false}
                        ]
                    });
                } catch (e) {
                    console.warn('KaTeX render warning:', e);
                }
            }
        });
        this.dom.btnCloseModal.addEventListener('click', () => this.dom.theoryModal.classList.remove('visible'));
        this.dom.theoryModal.addEventListener('click', (e) => {
            if (e.target === this.dom.theoryModal) this.dom.theoryModal.classList.remove('visible');
        });

        this.dom.btnToggleControls.addEventListener('click', () => {
            this.dom.controlsSidebar.classList.toggle('collapsed');
        });
        this.dom.btnToggleHud.addEventListener('click', () => {
            this.dom.telemetryHud.classList.toggle('collapsed');
        });

        this.dom.sliderC1.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            this.solver.c1 = val;
            this.dom.valC1.textContent = (val >= 0 ? '+' : '') + val.toFixed(2);
        });

        this.dom.sliderC3.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            this.solver.c3 = val;
            this.dom.valC3.textContent = (val >= 0 ? '+' : '') + val.toFixed(2);
        });

        this.dom.sliderMu.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            this.solver.mu = val;
            this.dom.valMu.textContent = val.toFixed(2);
        });

        this.dom.sliderDt.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            this.solver.dt = val;
            this.dom.valDt.textContent = val.toFixed(3);
        });

        this.dom.btnPlayPause.addEventListener('click', () => {
            this.paused = !this.paused;
            this.dom.btnPlayPause.textContent = this.paused ? 'RESUME' : 'PAUSE';
        });

        this.dom.btnStep.addEventListener('click', () => {
            this.solver.step();
            this.renderer.render();
            this.updateTelemetry();
        });

        this.dom.btnReset.addEventListener('click', () => {
            this.loadPreset(this.dom.presetSelect.value);
        });

        this.dom.btnInjectNoise.addEventListener('click', () => {
            this.solver.injectNoise(0.35);
        });

        this.dom.selectColormap.addEventListener('change', (e) => {
            this.renderer.colormap = e.target.value;
        });

        this.dom.checkDefects.addEventListener('change', (e) => {
            this.renderer.showDefects = e.target.checked;
        });

        this.dom.selectTool.addEventListener('change', (e) => {
            this.toolMode = e.target.value;
        });

        this.dom.sliderToolRadius.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            this.toolRadius = val;
            this.dom.valToolRadius.textContent = val.toFixed(1);
        });

        // Mouse Drag Perturbations
        let isInteracting = false;
        let lastPhysX = 0;
        let lastPhysY = 0;

        const getPhysCoords = (clientX, clientY) => {
            const rect = this.canvas.getBoundingClientRect();
            const width = this.canvas.width;
            const height = this.canvas.height;
            const size = Math.min(width, height) * 0.92;
            const startX = (width - size) * 0.5;
            const startY = (height - size) * 0.5;

            const px = (clientX - rect.left) * (width / rect.width);
            const py = (clientY - rect.top) * (height / rect.height);

            const xNorm = (px - startX) / size;
            const yNorm = (py - startY) / size;

            return {
                x: xNorm * this.solver.L - this.solver.L * 0.5,
                y: yNorm * this.solver.L - this.solver.L * 0.5,
                valid: xNorm >= 0.0 && xNorm <= 1.0 && yNorm >= 0.0 && yNorm <= 1.0
            };
        };

        this.canvas.addEventListener('mousedown', (e) => {
            const pos = getPhysCoords(e.clientX, e.clientY);
            if (pos.valid) {
                isInteracting = true;
                lastPhysX = pos.x;
                lastPhysY = pos.y;
                if (this.toolMode === 'amplitude_pin') {
                    this.solver.pinAmplitudeZero(pos.x, pos.y, this.toolRadius);
                }
            }
        });

        window.addEventListener('mousemove', (e) => {
            if (!isInteracting) return;
            const pos = getPhysCoords(e.clientX, e.clientY);
            if (pos.valid) {
                if (this.toolMode === 'razor') {
                    this.solver.cutPhaseLine(lastPhysX, lastPhysY, pos.x, pos.y);
                } else if (this.toolMode === 'amplitude_pin') {
                    this.solver.pinAmplitudeZero(pos.x, pos.y, this.toolRadius);
                } else if (this.toolMode === 'noise') {
                    this.solver.injectNoise(0.15);
                }
                lastPhysX = pos.x;
                lastPhysY = pos.y;
            }
        });

        window.addEventListener('mouseup', () => {
            isInteracting = false;
        });

        this.canvas.addEventListener('touchstart', (e) => {
            if (e.touches.length === 1) {
                const pos = getPhysCoords(e.touches[0].clientX, e.touches[0].clientY);
                if (pos.valid) {
                    isInteracting = true;
                    lastPhysX = pos.x;
                    lastPhysY = pos.y;
                }
            }
        }, { passive: true });

        window.addEventListener('touchmove', (e) => {
            if (isInteracting && e.touches.length === 1) {
                const pos = getPhysCoords(e.touches[0].clientX, e.touches[0].clientY);
                if (pos.valid) {
                    if (this.toolMode === 'razor') {
                        this.solver.cutPhaseLine(lastPhysX, lastPhysY, pos.x, pos.y);
                    } else if (this.toolMode === 'amplitude_pin') {
                        this.solver.pinAmplitudeZero(pos.x, pos.y, this.toolRadius);
                    }
                    lastPhysX = pos.x;
                    lastPhysY = pos.y;
                }
            }
        }, { passive: true });

        window.addEventListener('touchend', () => {
            isInteracting = false;
        });

        window.addEventListener('keydown', (e) => {
            if (e.code === 'Space') {
                e.preventDefault();
                this.paused = !this.paused;
                this.dom.btnPlayPause.textContent = this.paused ? 'RESUME' : 'PAUSE';
            } else if (e.key === 'r' || e.key === 'R') {
                this.loadPreset(this.dom.presetSelect.value);
            } else if (e.key === 'h' || e.key === 'H') {
                this.dom.telemetryHud.classList.toggle('collapsed');
            } else if (e.key === 't' || e.key === 'T') {
                this.dom.theoryModal.classList.toggle('visible');
            } else if (e.key === 'Escape') {
                this.dom.theoryModal.classList.remove('visible');
            }
        });
    }

    resize() {
        this.canvas.width = window.innerWidth;
        this.canvas.height = window.innerHeight;
    }

    loadPreset(name) {
        this.solver.initPreset(name);
        this.syncUi();
    }

    syncUi() {
        this.dom.sliderC1.value = this.solver.c1.toFixed(2);
        this.dom.valC1.textContent = (this.solver.c1 >= 0 ? '+' : '') + this.solver.c1.toFixed(2);

        this.dom.sliderC3.value = this.solver.c3.toFixed(2);
        this.dom.valC3.textContent = (this.solver.c3 >= 0 ? '+' : '') + this.solver.c3.toFixed(2);

        this.dom.sliderMu.value = this.solver.mu.toFixed(2);
        this.dom.valMu.textContent = this.solver.mu.toFixed(2);

        this.dom.sliderDt.value = this.solver.dt.toFixed(3);
        this.dom.valDt.textContent = this.solver.dt.toFixed(3);
    }

    updateTelemetry() {
        const diag = this.solver.computeDiagnostics();
        const detected = this.defectTracker.findDefects(this.solver);

        const isUnstable = diag.bfn < 0;
        this.dom.telBfn.textContent = `${diag.bfn.toFixed(3)} (${isUnstable ? 'UNSTABLE' : 'STABLE'})`;
        this.dom.telBfn.className = `metric-value ${isUnstable ? 'magenta' : 'green'}`;

        this.dom.telDefects.textContent = detected.count;
        this.dom.telCharge.textContent = (detected.netCharge > 0 ? '+' : '') + detected.netCharge;
        this.dom.telMeanAmp.textContent = diag.meanAmplitude.toFixed(3);
        this.dom.telMaxAmp.textContent = diag.maxAmplitude.toFixed(3);
        this.dom.telOmega.textContent = `${diag.meanOmega.toFixed(2)} rad/s`;
        this.dom.telTime.textContent = `t = ${diag.time.toFixed(2)}`;

        if (this.audio.isPlaying) {
            this.audio.update(diag.meanOmega, detected.count);
        }
    }

    updateFps(now) {
        const delta = now - this.lastFpsUpdate;
        if (delta > 0) {
            const fps = 1000.0 / delta;
            this.fpsHistory.push(fps);
            if (this.fpsHistory.length > 20) this.fpsHistory.shift();
            const avg = Math.round(this.fpsHistory.reduce((a, b) => a + b, 0) / this.fpsHistory.length);
            this.dom.telFps.textContent = `${avg} FPS`;
        }
        this.lastFpsUpdate = now;
    }

    start() {
        const loop = (now) => {
            requestAnimationFrame(loop);

            if (!this.paused) {
                this.solver.step();
            }

            this.renderer.render();
            this.updateTelemetry();
            this.updateFps(now);
        };

        requestAnimationFrame(loop);
    }
}

if (typeof window !== 'undefined') {
    window.addEventListener('DOMContentLoaded', () => {
        const app = new CGLEApp();
        app.start();
    });
}
