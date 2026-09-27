# Simulation redesign: showing collective predator response honestly

Status: **implemented** (designed 2026-09-26, built and calibrated 2026-09-27). Section 10 records what
was built, how it was calibrated, the results, and where the implementation departs from this design. It builds on
[CRITICALITY_AND_DYNAMICS.md](CRITICALITY_AND_DYNAMICS.md), which explains why the current simulation
can't show the intended behavior. The literature behind it is summarized here. More detailed
AI-generated reading notes are kept locally in `docs/research/`; they are deliberately not committed
(gitignored) because no human has reviewed them.

## 1. The point we want to get across

Agreed statement:

> With near-critical speed control and inertial turning dynamics, one bird's reaction to a predator
> spreads across the whole flock: the turn travels as an undamped wave and the speed change becomes
> flock-wide. Stiffen individual speed control, or make the turning overdamped, and the same reaction
> stays local and dies out.

This gives two independent contrasts (turning dynamics; speed control), plus a trade-off we show
honestly (responsiveness also spreads false alarms).

## 2. What the research supports

Labels used throughout: **[data]** measured in real flocks; **[model]** a published model result;
**[new]** published recently, or a preprint, and not yet independently tested; **[ours]** our own
modeling choice.

| Claim | Status |
|---|---|
| Birds interact with ~6–8 topological neighbors | [data] Ballerini 2008; review 2018 |
| Direction and speed correlations reach across the whole flock (ξ ∝ L) | [data] Cavagna 2010, 2022 |
| Scale-free *direction* correlations need no tuning (Goldstone modes) | [model] standard physics; review 2018 |
| Scale-free *speed* correlations need special speed control: near-critical g (Bialek 2014) or flat-bottomed "marginal" confinement (Cavagna 2022) | [data] + [model]; mechanism debated |
| Turns start with a few birds and propagate linearly in time with negligible attenuation, at 20–40 m/s (flock speed 7–12 m/s) | [data] Attanasi 2014 (published values; the preprint table lists 9.4–21.3 m/s) |
| Overdamped dynamics (Vicsek, and the paper's Appendix G) can't reproduce this; inertial turning (ISM) can | [model] Attanasi 2014; Cavagna 2015 |
| *Spontaneous* fluctuations in real flocks are overdamped (no spin-wave peaks) | [data][new] Cavagna et al., PRL 2026 (arXiv:2505.19665) |
| A quartic alignment term reconciles overdamped fluctuations with propagating turns | [model][new] same paper; authors call their nonlinear analysis "very crude"; untested in 3D with noise |
| Only birds near the predator detect it; reaction delay ~50–90 ms | [data] Hemelrijk 2015; Papadopoulou 2026 |
| Critical speed control improves predator evasion | **Untested.** Bialek et al.: "seems advantageous … it remains to be seen" |
| Being at *other* critical points is optimal against predators | **Contradicted/qualified**: Klamser 2021 (escape improves with stronger order); Poel 2022 (critical helps only at intermediate risk; false alarms spread too) |

Consequence: the simulation can show what *the model* predicts about predator response under each
condition. It must not present "critical flocks evade better" as a known fact about starlings.

## 3. Model

### 3.1 State and units

Physical units: meters and seconds. Each bird i has a position x_i, a unit heading u_i, a speed s_i
(velocity v_i = s_i u_i) and a spin σ_i ⊥ u_i (its rotational momentum; its turn rate is |σ_i|/χ).

### 3.2 Energy

For a fixed neighbor network, all conservative interactions come from one energy function:

```
H = H_dir + H_sp − Σ_i h_i·u_i − Σ_i b_i s_i

H_dir = (J/4)  Σ_ij n_ij |u_i − u_j|²  +  (J4/4) Σ_ij n_ij |u_i − u_j|⁴        (alignment; J4 ≥ 0)
H_sp  = (J/4)  Σ_ij n_ij (s_i − s_j)² / v0²  +  Σ_i V(s_i)                      (speed alignment + control)
V(s)  = (g/2)(s − v0)²/v0²              linear control (Bialek 2014)
      = (λ/v0⁶)(s² − v0²)⁴              marginal control (Cavagna 2022; their exact form)
```

- n_ij = (n̂_ij + n̂_ji)/2 is the symmetrized topological neighbor matrix, as in Bialek 2014.
- Splitting direction and speed follows the paper's own small-fluctuation decomposition (its Eq. 5),
  with one alignment strength J for both, as in its "unified" model.
- h_i and b_i are external fields. They are nonzero only for birds reacting to a predator (§4).
  This is the same mathematical role that boundary birds play in the paper's Appendix D.2.
- The quartic term J4 is the PRL 2026 addition. With J4 = 0 we recover the paper's alignment exactly
  (at fixed speed).

### 3.3 Turning dynamics: one equation family, three presets

Inertial Spin Model (Cavagna 2015), written with a unit heading:

```
du_i/dt = (1/χ) σ_i × u_i
dσ_i/dt = u_i × F_i  −  (η/χ) σ_i  +  u_i × ξ_i            ⟨ξ ξ⟩ = 2 η T
F_i = −∂H/∂u_i + f_i^coh    (pairwise:  Σ_j n_ij u_j [J + 4 J4 (1 − u_i·u_j)] + h_i + f_i^coh)
```

Its overdamped limit (χ → 0) is the paper's Appendix G / Vicsek-type dynamics:

```
η du_i/dt = P⊥(u_i) F_i + noise
```

| Preset | Setting | Matches | Contradicts | Label |
|---|---|---|---|---|
| **A. Overdamped** (current code / paper App. G) | χ → 0, J4 = 0 | static correlations | turns die out and flocks split instead of turning (Attanasi 2014; Cavagna 2015) | [model] reference/contrast |
| **B. Inertial** (ISM) | low friction η, J4 = 0 | linear, undamped turn propagation | spontaneous fluctuations oscillate, whereas real ones are overdamped (PRL 2026) | [model] |
| **C. Nonlinear inertial** (ISM + quartic) | high η, J4 > 0 | overdamped spontaneous fluctuations; large turns propagate | not yet demonstrated for realistic polarization or in 3D with noise (see §6) | [new] experimental |

The GUI switches presets; the underlying parameters are also exposed.

### 3.4 Speed dynamics

Speed dynamics have never been measured in real flocks, so this part is **[ours]**. We use the simplest
choice consistent with the energy: overdamped Langevin dynamics on the speed.

```
γ_s ds_i/dt = −∂H/∂s_i + 2T/(N s̄) + f_i^coh·u_i + ξ_i^s      ⟨ξ^s ξ^s⟩ = 2 γ_s T
```

The 2T/(N s̄) term is the 3D volume factor (the "entropic push" described in Cavagna 2022). In the
full-velocity models it acts once, on the flock's mean velocity: each bird's own factor cancels against
its sideways velocity fluctuations. (The first implementation applied it per bird, which blew up the
group speed under weak speed control; see section 10.) Speed control is switchable:

| Speed preset | Setting | Behavior |
|---|---|---|
| **S1. Stiff** | linear, g/(J·n_c) ≈ 1 | speed changes stay local (non-critical) |
| **S2. Near-critical linear** | linear, g/(J·n_c) ≈ 10⁻³ (Bialek 2014) | flock-wide speed response; known side effect: group speed drifts or inflates in small flocks (Cavagna 2022) |
| **S3. Marginal** | quartic V, λ from Cavagna 2022 | flock-wide speed correlations with bounded speeds **[data-backed statics]** |

### 3.5 Cohesion

We keep the paper's Appendix G attraction/repulsion forces (Eq. G5), rescaled so that the nearest-
neighbor distance is about 1 m and the hard core about 0.4 m, matching real starlings. Each force f_i
acts on heading as a torque (u_i × f_i) and on speed through its forward component (f_i·u_i).
These forces are not derived from H (as in the paper), so they are a non-equilibrium ingredient; we
label them as such. **[model]**

### 3.6 Neighbors

We keep topological neighbors with n_c ≈ 7. There are two options: the existing "balanced" rule
(Camperi 2012; its culling loop has an iterator bug to fix), or plain n_c-nearest. Either way the
matrix is symmetrized. Neighbors are recomputed every frame. Real neighbor turnover takes seconds, so
this is effectively a slowly varying network.

### 3.7 Formal consistency with the paper

For a fixed neighbor network, with no cohesion forces and no predator, **all three presets have the
same stationary distribution P ∝ exp(−H/T)**. With J4 = 0 that is the paper's maximum-entropy
distribution: the ISM's spin part factors out (Cavagna 2015, Eq. 42), and the speed part includes the
3D volume factor.

So the presets differ only in *dynamics*: how disturbances travel, not what snapshots look like. This
makes the comparison fair, and it is exactly the distinction in CRITICALITY_AND_DYNAMICS.md. Moving
birds, cohesion forces and the predator make the real system non-equilibrium, so this holds only
approximately. We verify it numerically (§5).

## 4. Predator

- **Motion** **[ours, informed by data]**: the predator shadows the flock, attacks at 1.3–2× flock speed
  (Papadopoulou 2026), passes through, and retreats. Attacks are triggered by a key, or run on an
  automatic cycle. Manual steering is optional.
- **Detection** **[data]**: only birds within a detection radius R_d of the predator sense it
  directly. Everyone else learns only through the model's own couplings. Detection starts after a
  reaction delay of about 70 ms (50–90 ms reported).
- **Response** **[model]**: a detected bird gets a sustained external field for as long as it
  detects the predator:
  - a heading field h_i = h_p · (direction away from the predator), i.e. a torque turning it away;
  - an optional speed field b_i > 0, i.e. an escape speed-up.

  In the energy these are just the −h·u and −b·s terms, so the predator enters the model the same
  way the paper handles boundary birds. The field is sustained because the PRL 2026 analysis
  suggests large turns propagate only while they are driven.
- **False alarms** (toggle) **[data-motivated]**: occasionally one random bird receives the same
  field for a short time with no predator present. In fish schools, false alarms spread like real
  ones (Rosenthal 2015), and a more responsive group pays for it in false alarms (Poel 2022). Showing
  this is part of being honest about the trade-off.
- **Outcome measures** (labeled as model predictions):
  - fraction of the flock that has turned within 1 s;
  - lowest polarization during the escape;
  - number of fragments (connected components of the neighbor graph);
  - turn-front speed and attenuation;
  - captures (predator within a capture radius).

  Klamser 2021 warns that capture counts are confounded by where birds sit in the flock, so captures
  are secondary.
- **Not modeled** (documented as out of scope for now): StarDisplay-style copying of a discrete escape
  maneuver, which produces the dark "agitation" (banking) waves seen in starlings. It is a different
  mechanism from alignment and isn't part of the max-ent framework. It could become an optional,
  labeled channel later.

## 5. Calibration targets and live validation panel

Our flock is ~300 birds; real flocks are 10–4000. We match dimensionless ratios, not raw sizes.

| Observable | Target | Source | How we check |
|---|---|---|---|
| Polarization Φ | 0.96 ± 0.03 (0.84–0.995) | Bialek 2014 Table I; review 2018 | live |
| Interacting neighbors | 6–8 | Ballerini 2008 | live |
| Group speed | ~12 m/s, independent of N | review 2018; Cavagna 2022 | live + N sweep |
| Individual speed spread SD/V | ~0.07–0.2 (debated) | Bialek 2014 Fig. 1a; review figures | live |
| Correlation length ξ (zero crossing) of C_dir, C_sp | ∝ L across N (coefficient 0.1–0.4 depending on L's definition) | Cavagna 2010, 2022 | live + N = 100/300/1000 headless sweep |
| Turn front | linear in time, negligible attenuation; 20–40 m/s, i.e. ~2–5× flock speed; faster in more ordered flocks | Attanasi 2014 | measured automatically after each attack |
| Spontaneous fluctuations | overdamped (no oscillation peak) | PRL 2026 | autocorrelation of heading fluctuations |
| Alignment relaxation / neighbor turnover | ~0.1–0.5 s / ~3–20 s | Mora 2016 | headless |
| Stationary statics identical across presets | same Φ, C(r) within noise | §3.7 | headless |

Each row shows pass/fail or not-applicable in the GUI. Clicking a row shows its source. A preset is
described as "consistent with the research" only for the rows it passes, and the panel shows which rows
fail. Known mismatch we won't fix: real flocks are flat (axis ratios about 1 : 2.8 : 5.6), while model
flocks like ours are round.

Free parameters (χ, η, J, J4, T, γ_s, g or λ, R_d, h_p, b) are chosen per preset by headless calibration
runs against this table. The chosen values and each one's justification go in the docs.

## 6. Open risks

1. **Preset C may not work in our setting.** The PRL analysis relies on extremely polarized flocks (1−Φ
   ~ 10⁻⁷) and driven turns. At realistic polarization its own constraints limit the quartic term to
   J4 ≲ J/(1−Φ). A free turn pulse is also drained by friction faster than it can cross the flock.
   Whether a *predator-driven* turn propagates in this regime is an open question we would be testing,
   not assuming. If it fails, we report that, and presets A vs. B remain the main contrast.
2. **Measured turn speeds:** the published paper gives 20–40 m/s; the preprint table and the group's
   2026 paper cite 10–20 m/s. We target the published range and note the discrepancy.
3. **Speed contrast may be subtle:** turning stays coherent in the ordered phase whatever g is. The
   speed contrast will be shown with a speed-deviation color mode and a C_sp(r) plot, not by expecting
   dramatically different escapes.
4. **The ISM formula for wave speed overestimated its own simulation** (Cavagna 2015). We measure
   the wave speed rather than trust the formula.

## 7. Numerical scheme

- Heading and spin: a symplectic splitting integrator (BAOAB-style):
  - the kick step adds torque to spin;
  - the drift step rotates the heading exactly about σ, so |u| = 1 is preserved;
  - the friction+noise step solves the spin's Ornstein–Uhlenbeck process exactly: σ ← e^(−ηdt/χ) σ
    plus noise of variance χT(1 − e^(−2ηdt/χ)), which stays stable at high friction.

  The spin is re-projected perpendicular to u after each step.
- Preset A: Euler–Maruyama on the sphere, then renormalize.
- Speed: Euler–Maruyama, with substeps and clamping for the stiff quartic potential.
- Time step ~1 ms, with several substeps per 60 fps frame. At 300 birds × ~7 neighbors that is cheap.
- A **time-scale slider** (slow motion), since a real turn crosses the flock in under a second.

## 8. Implementation plan

1. **New simulation core** in `src/sim/`: a clean C++ module (state, neighbors, energy/forces,
   integrators, predator, measurements) replacing the legacy `include/swarming_spp` code. The legacy
   code stays until preset A reproduces its statistics (a parity check), then it is removed.
2. **Measurement and headless mode:** correlation functions, ξ, polarization, turn-front tracker,
   fluctuation autocorrelation, and a command-line batch mode for calibration and N sweeps.
3. **Calibration** of presets A, B, C and S1–S3 against §5; record the parameter values and
   justifications.
4. **Predator and false alarms.**
5. **UI:** presets, parameter panel with honesty labels, validation panel, color modes (turn rate /
   speed deviation / detection state), predator rendering, time scale.
6. **Docs:** update this file with what was built, the calibrated values, and the validation results.

## 9. Sources

Primary papers (per-paper details and verification status are in the local, uncommitted `docs/research/` notes):

- Bialek et al., arXiv:1307.5563; PNAS 111, 7212 (2014)
- Ballerini et al., PNAS 105, 1232 (2008), arXiv:0709.1916
- Cavagna et al., "Scale-free correlations in starling flocks", PNAS 107, 11865 (2010), arXiv:0911.4393
- Attanasi et al., Nature Physics 10, 691 (2014), arXiv:1303.7097
- Cavagna et al., "Flocking and turning", J. Stat. Phys. 158, 601 (2015), arXiv:1403.1202
- Cavagna et al., Phys. Biol. 13, 065001 (2016), arXiv:1605.09628
- Mora et al., "Local equilibrium in bird flocks", Nature Physics 12, 1153 (2016), arXiv:1511.01958
- Cavagna, Giardina, Grigera, "The physics of flocking", Physics Reports 728, 1 (2018)
- Cavagna et al., "Marginal speed confinement…", Nature Communications 13, 2315 (2022), arXiv:2101.09748
- Cavagna et al., fully conserved ISM, J. Phys. A 57 (2024), arXiv:2403.07644
- Cavagna et al., "Spin-Waves without Spin-Waves…", arXiv:2505.19665; Phys. Rev. Lett. 137, 098301 (2026)
- Hemelrijk et al., Behav. Ecol. Sociobiol. 69, 755 (2015) (StarDisplay agitation waves)
- Papadopoulou et al., Communications Biology (2026) (StarEscape)
- Rosenthal et al., PNAS 112, 4690 (2015); Poel et al., Science Advances 8, eabm6385 (2022)
- Klamser & Romanczuk, PLoS Comput. Biol. 17, e1008832 (2021)

## 10. Implementation and calibration record (2026-09-27)

### 10.1 What was built

- `src/sim/`: a new simulation core with no Cinder dependency.
  - `Flock`: state, topological neighbors, forces, and the BAOAB and overdamped integrators.
  - `Params`: parameters and the presets.
  - `Predator`: attack runs, detection with reaction delay, fields, false alarms, per-attack report.
  - `Measure`: statistics, correlation functions, turn-front tracker, fluctuation probe.
- `src/FlockingApp.cpp`: the app, rewritten on the new core. It has the presets, predator controls,
  color modes (including "turned since attack"), slow motion, and the live validation panel.
- `tools/flocksim_cli.cpp`: a headless driver with `stats`, `turn`, `speedpush`, `attack` and `probe`
  commands. All numbers below come from it.
- The legacy simulation code (`include/swarming_spp/`, `Particle`, `ParticleController`) was removed
  after the parity check (10.2).

### 10.2 Verification

- **Equilibrium test (section 3.7):** with a fixed neighbor network and no cohesion, the overdamped
  preset and the inertial preset at η = 0.5 and η = 20 all give polarization 0.9923–0.9927 at the same
  J and T. The friction differs by 280× across these runs. So the presets differ only in dynamics, as
  designed, and the integrators are correct.
- **Parity with the legacy code:** with the legacy parameters mapped into the new core (overdamped),
  the results match the old app.

  | | Legacy app | New core |
  |---|---|---|
  | Neighbors | ~7.0 | 7.03 |
  | Q_int | 0.012–0.014 | 0.0130 |
  | Mean speed | 12.2–12.5 m/s | 12.34 m/s |
  | Polarization | 0.989 | 0.995 |

  The small polarization gap is attributed to the legacy code's coarse time step (dt = 0.03) and to
  readings taken before it had fully settled.

### 10.3 Calibration and departures from the design

Final values are in `src/sim/Params.cpp`, with a comment on each. The points that needed real decisions:

1. **Inertial steering plus cohesion is unstable at low friction.** A bird steers toward its neighbors,
   moves, and overshoots. With inertia and little friction this relative motion grows: the dynamics are
   third order, and the Routh–Hurwitz condition needs roughly η·J·n_c > χ·K·v0, where K is the cohesion
   stiffness. At η = 0.5 the flock heated up to polarization 0.63 and fell apart.
   - Neither ingredient does this alone: moving birds without cohesion stay at 0.9915, and cohesion
     without movement stays at 0.9914.
   - Fix: η ≥ 1 for preset B (we use 2), and a moderate cohesion core.
2. **J was raised from 100 to 200** so that preset B's turn fronts land in the measured range (a
   brief pulse gives about 14 m/s, a sustained push about 22 m/s; data: 20–40 m/s). J/T = 10 gives
   polarization about 0.975.
3. **Cohesion** keeps the shape of eq. G5 but uses rEq = 1.3 m, rAttract = 1.9 m, far attraction 0.15
   (1 in G5) and core strength 10. G5's strong far attraction compressed the flock to r1 ≈ 0.3 m; real
   r1 is 0.68–1.51 m.
   - A fore-aft cohesion gain of 30 in the speed equation **[ours]** keeps the flock together under
     weak speed control. Birds can hold their place along the flight direction only by changing speed.
   - Cohesion strength is calibrated **per preset** (A 320, B 40, C 110) so every flock has similar size
     and spacing at rest. That's a documented difference beyond the turning dynamics.
4. **Preset A's friction (η = 37 = √(J·n_c·χ))** matches B's local response time, making inertia the
   only intended difference. A tension appears here:
   - An overdamped flock tuned to the *measured* spontaneous relaxation (η ≈ 840) cannot stay cohesive
     at all, even with 12× cohesion.
   - A much faster-steering one (η = 15) does carry turns across the flock, but only with relaxation
     10–20× faster than measured.
5. **The speed volume factor is collective** (section 3.4). The per-bird version pushed the group speed
   to 33 m/s under near-critical control.
6. **The predator strikes at the edge and breaks away.** It breaks away after closing to 1.5 m, or 0.6 s
   after the first detection. The first version plowed through the flock, so up to 40% of birds sensed
   it directly, which hid the social propagation the demo is about.
7. **Preset C:** at J = 200 we use η = 60, J4 = 30000. We found no parameters that give realistic
   spontaneous relaxation, realistic polarization *and* propagating turns at once:
   - lowering J to slow the relaxation stopped turns from spreading (1–18% turned);
   - the quartic term pushes polarization up to about 0.997.

### 10.4 Results (N = 300 unless noted; headless runs)

**Snapshot statistics** at rest, speed control S3 (marginal):

| Preset | Polarization | Speed (m/s) | Speed SD/mean | Neighbors | r1 (m) | L (m) | Q_int | ξ_dir/L | ξ_sp/L | Pieces |
|---|---|---|---|---|---|---|---|---|---|---|
| A overdamped | 0.978 | 11.0 | 0.16 | 7.0 | 0.83 | 21 | 0.051 | 0.25 | 0.20 | 1 |
| B inertial | 0.975 | 12.2 | 0.14 | 7.0 | 0.98 | 28 | 0.049 | 0.23 | 0.21 | 1 |
| C nonlinear | 0.997 | 12.4 | 0.14 | 6.8 | 0.94 | 20 | 0.019 | 0.31 | 0.25 | 1 |

All rows pass every static check except C's polarization, which is above the 0.995 upper end.
Stiff (S1) and near-critical (S2) speed control are similar, except for speed spread: S1 is about 0.09,
S2 about 0.14–0.17.

**Turn started by 3 edge birds** (heading push 1 × J·n_c; S3):

| Preset | Brief 0.3 s push: share turned in 3 s | Front speed | Sustained push: share turned | Front speed |
|---|---|---|---|---|
| A overdamped | 54% | not measurable (partial) | 61% | not measurable |
| B inertial | 100% | ~14 m/s | 100% | ~22 m/s |
| C nonlinear | 1% | — | 30–39% | ~3 m/s |

**Predator attacks**, 6 random attacks per preset (S3):

| Preset | Birds sensing the predator directly | Share turned > 30° within 2 s | Front speed | Most pieces |
|---|---|---|---|---|
| A overdamped | 3–39% | 69–100% (100% only when ≥ 20% sensed it directly) | 6–17 m/s | 2 |
| B inertial | 4–19% | 100% in every attack | 27–46 m/s | 1 |
| C nonlinear | 6–20% | 2–100% (all-or-nothing) | 12–28 m/s when it spreads | 1 |

**Speed push:** 3 edge birds pushed to speed up for 4 s (preset B):

| Speed control | Change near the pushed birds | Whole-flock mean change |
|---|---|---|
| S1 stiff | +0.2 m/s | −0.04 m/s (stays local) |
| S2 near-critical | +2.4 m/s, +0.6 m/s at 9 m | +0.98 m/s (flock-wide) |
| S3 marginal | +2.6 m/s, +0.5–0.85 m/s out to 12 m | +0.69 m/s (flock-wide) |

**Flock-size sweep** (preset B):

| Speed control | N = 100 | N = 300 | N = 1000 |
|---|---|---|---|
| S3 marginal | ξ_sp = 2.7 m, L = 17 m | 5.7 m, L = 28 m | 11.2 m, L = 53 m |
| S1 stiff | 1.6 m, L = 10 m | 2.1 m, L = 14 m | 2.8 m, L = 20 m |

S3's ξ_sp/L stays at 0.16–0.21. But S1's ratio is similar (0.14–0.16), because the stiff flocks stay
compact. This is the review's warning in action: at these sizes the zero crossing can't tell the two
apart (it's forced by a sum rule). The speed-push response test above is what separates them.

**Spontaneous fluctuations** at k = 1.5 m⁻¹ (target: overdamped, half-life ~0.2–0.4 s):

| Preset | Character | Half-life |
|---|---|---|
| A overdamped | overdamped | 0.04 s |
| B inertial | **oscillatory** (autocorrelation dips to −0.42) | 0.04 s |
| C nonlinear | overdamped | ≤ 0.02 s |

All three presets fail the half-life target.

### 10.5 What the simulation shows, and what it doesn't

- **Shown, within the model:**
  - With inertial turning, a threat sensed directly by a few birds turns the *whole* flock, fast
    (27–46 m/s, about the measured 20–40 m/s), and the flock stays in one piece.
  - With overdamped turning at the same statics and local response time, the same threat turns fewer
    birds, more slowly, and the flock sometimes splits. A brief stimulus mostly stays local.
  - With near-critical or marginal speed control (the paper's criticality), an escape speed-up spreads
    across the flock; with stiff control it stays with the birds that sensed the threat.
- **Not shown (and not claimed):**
  - That near-critical flocks evade real predators better. That is untested in the literature. The
    predator panel labels its outcomes as model predictions.
  - That any preset is "the" correct model. Each fails at least one published measurement, as listed
    above and shown live in the validation panel. The literature's tension (overdamped spontaneous
    fluctuations versus undamped turn propagation) is reproduced here, not resolved.
- **Known limitations:**
  - Model flocks are round; real ones are flat.
  - Speed dynamics are our assumption; they have never been measured.
  - Cohesion is non-equilibrium and calibrated per preset.
  - Turn-front speeds are noisy.
  - There is no copying of discrete escape maneuvers, so no agitation waves.
