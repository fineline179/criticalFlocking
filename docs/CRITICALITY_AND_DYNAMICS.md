# Criticality vs. dynamics: what the paper shows and what the simulation needs

Written 2026-09-26 while revisiting this project. It records an analysis of what
[Bialek et al. 2013 (arXiv:1307.5563)](https://arxiv.org/abs/1307.5563) actually establishes, compared
with the goal the simulation was originally written for.

The resulting redesign is in [SIMULATION_DESIGN.md](SIMULATION_DESIGN.md).

## The original goal

Build a real-time simulation where, with the flock's parameters tuned to the paper's critical point, a
predator near one part of the flock makes the whole flock respond coherently ("all as one"). Tuned away
from criticality, that coherent response should disappear and the flock should evade the predator less
effectively.

The concern that prompted this analysis: the paper's equations seem to describe the equilibrium
(stationary) distribution for given parameters, not a recipe for a dynamical simulation that shows the
"all as one" behavior.

## What the paper establishes

1. **The main result is static.** The maximum-entropy model (Sections II–IV) is a probability distribution
   over the birds' velocities in a *single snapshot*, with positions taken as given. Its parameters J, g
   and n_c are inferred from single snapshots. It predicts equal-time correlations, C_dir(r) and C_sp(r).
   It contains no time at all: no time scales, no propagation speed, no order of events.
2. **Direction and speed are scale-free for different reasons.**
   - *Direction:* long-range correlations come from Goldstone modes, meaning a broken rotational symmetry
     (the flock picked a heading). They occur anywhere in the ordered (flocking) phase, with no tuning.
     Real flocks are deep in that phase (polarization 0.84–0.99), not at the order–disorder transition.
   - *Speed:* choosing a speed breaks no symmetry, so scale-free speed correlations need tuning, with
     g/(J·n_c) → 0. Real flocks sit at about 10⁻³; this is the paper's "criticality." Formally it is a
     Gaussian field approaching its massless point: ξ ~ r_c·√(J·n_c/g) diverges, and so does the speed
     variance. It is not a boundary between two phases that you cross.
3. **In a finite flock, criticality is a regime, not a knife-edge.** Once ξ exceeds the flock size L,
   lowering g further changes nothing (the curves "pile up" in Fig. 3a).
4. **The dynamics (Section V, Appendix G) are an assumption, not a finding.** The authors chose overdamped
   Langevin motion on the model's energy as a "natural" dynamics. They checked it only against static
   quantities (Fig. 4: C_sp(r) and ξ versus g), never against how disturbances travel over time.

## Corrections to the original framing

1. **The code already is a real-time dynamical simulation.** `Bialek_consensus::sense_velocity`
   (`include/swarming_spp/behavior.cpp`) integrates G1–G2 forward in time: an Euler step with noise of
   variance 2T·dt. Nothing converges and then stops. In the stationary state every bird keeps moving and
   fluctuating, and the paper's correlations are statistics of those ongoing fluctuations. With positions
   frozen, these dynamics sample the max-ent distribution exactly. With birds moving, neighbors changing,
   and the cohesion forces f_ij (which aren't derived from the energy), they sample it only approximately.
   So the real question isn't "equilibrium vs. dynamics"; it's whether *these particular* dynamics show
   the behavior we want.
2. **The static result does say something about response, via fluctuation–dissipation.** For a stationary
   distribution P ∝ e^(−H/T), the steady response of bird j to a sustained push on bird i equals their
   equal-time correlation divided by T. "Speed fluctuations are correlated across the whole flock" is
   therefore the same statement as "a sustained nudge to one bird's speed ends up shifting every bird's
   speed." That is the paper's own claim in Section VI. But it covers only the *range and size of the
   eventual response*, not how fast the response arrives or what it looks like on the way.
3. **Criticality doesn't make birds interact "super quickly"; near a critical point the dynamics slow
   down.** Interactions stay local: each bird still responds only to its ~n_c neighbors, and what becomes
   long-range is the correlation. In the overdamped G1 dynamics, a speed disturbance of wavelength ℓ
   relaxes at a rate proportional to g + J·n_c·r_c²/ℓ². At g → 0 the time scales as ℓ², so disturbances
   spread diffusively (distance ~ √t) and are damped. This is critical slowing down. Direction
   disturbances behave the same way in G1.
4. **Predator evasion is mostly a turn, and in this model the critical knob g doesn't control turning.**
   Directional correlations are scale-free throughout the ordered phase, whatever g is. Tuning g mainly
   controls whether *speed changes* spread collectively. Raising T toward the order–disorder transition
   doesn't give "optimal coherence" either; the flock loses its alignment.
5. **The same group later showed that Appendix G-type dynamics are wrong for turning.** Attanasi et al.
   (2014), cited as ref. [1] in this paper, measured real starling flocks turning. They found the turn
   propagates *linearly in time with negligible attenuation*, contradicting overdamped Vicsek-type models
   (a standard family of flocking models, Appendix G's among them), which "predict a slower and
   dissipative transport of directional information." The explanation is a conserved "spin," a bird's
   rotational momentum. That makes the dynamics second-order in time, so turns travel as waves. Their
   Inertial Spin Model (Cavagna et al. 2015) reduces to the Vicsek model in the overdamped limit, which
   "dissipates rotational information and does not allow for polarized turns." In that theory,
   propagation gets faster with stronger order, not with criticality. The swerve envisioned above, where
   one bird sees the predator and a wave sweeps the whole flock, is exactly what G1-type dynamics fail to
   produce.
6. **The speed-criticality story was later refined.** With Gaussian (quadratic) speed control there is a
   trade-off: small g gives scale-free correlations but huge speed fluctuations. Cavagna et al. (2022)
   proposed "marginal speed confinement": a speed potential that ignores small deviations but sharply
   suppresses large ones. It gives scale-free correlations with bounded speeds, without fine-tuning g,
   and they validated it on flocks spanning over two orders of magnitude in size.

## The simulation's settings never leave the critical regime

The code's defaults (J = 19, g = 0.2, n_c = 8) give g/(J·n_c) ≈ 1.3×10⁻³, the same as the paper's real
flocks. For a 300-bird flock, L ≈ N^(1/3)·r_c ≈ 7·r_c, so the flock is effectively critical whenever
g/(J·n_c) ≲ (r_c/L)² ≈ 0.02. The g slider goes only to 0.5, so almost the entire range of the
parameter panel is in the critical regime. Seeing the non-critical contrast would need g/(J·n_c) around
0.1–1, i.e. g in the tens to hundreds. These are order-of-magnitude estimates; the code and the paper
differ by factors of about 2 in how J is normalized.

## Summary

- The paper pins down only the static statistics. It shows that real flocks' speed control is in a critical
  regime, and that their directional correlations are scale-free by symmetry.
- Statics plus fluctuation–dissipation give only the range of a slow, sustained response. They don't tell
  you whether the flock visibly evades "as one" in real time.
- What decides that is the dynamics, which the paper leaves as a free modeling choice (many dynamics share
  the same stationary distribution). The choice in Appendix G is dissipative, and the same group later
  showed that real flocks don't turn that way.
- Gaps between the current simulation and the goal:
  - its g range never leaves speed criticality;
  - g governs speed coordination, not turning;
  - its overdamped dynamics damp and slow the very signal we want to watch;
  - it has no predator term.
- A predator would enter the model as an external push on nearby birds, the same way Appendix D.2 treats
  boundary birds acting on the interior. Real evasion is a large disturbance, though, so
  fluctuation–dissipation gives intuition there, not guarantees.
- A candidate direction: replace the dynamics while keeping the paper's energy. The Inertial Spin Model
  does this for direction: its velocities settle into the same stationary distribution, but turns
  propagate as waves. Speed would need its own extension.

## Sources

- W. Bialek et al., "Social interactions dominate speed control in driving natural flocks toward
  criticality," [arXiv:1307.5563](https://arxiv.org/abs/1307.5563) (2013).
- A. Attanasi et al., "Superfluid transport of information in turning flocks of starlings,"
  Nature Physics 10, 691 (2014), [arXiv:1303.7097](https://arxiv.org/abs/1303.7097).
- A. Cavagna et al., "Flocking and turning: a new model for self-organized collective motion,"
  J. Stat. Phys. 158, 601 (2015), [arXiv:1403.1202](https://arxiv.org/abs/1403.1202).
- A. Cavagna et al., "Marginal speed confinement resolves the conflict between correlation and control in
  natural flocks of birds," Nature Communications 13, 2315 (2022),
  [arXiv:2101.09748](https://arxiv.org/abs/2101.09748).
