# Expanding bound complexes and many-body terms (step 6)

## What is expanded

Only a bound complex K is expanded. The added unit X is either

- one further defect, placed near A, near B, or near both members of K
  (codopex criteria `nearA_farB`, `nearB_farA`, `both_nn`), or
- a second bound complex (two copies of K, or K with another bound complex).

K is treated as a **composite defect**: it carries its own charge state q_c and
chemical identity, which need not be the sum of those of its constituents. Score
the contacts that X makes with K using the rules of `scoring.md`, with the
composite's q_c and oxidation states.

Examples of why the composite matters:

- Pb3(PO4)2: Cu_Pb2 (q = 0) + S_O1 (q = 0) forms a complex with q = −1, Cu²⁺ being
  reduced to Cu⁺. A new defect sees a charged complex and a Cu⁺, not two neutral
  defects.
- CuInSe2: In_Cu²⁺ + V_Cu⁻ has a net +1 charge that attracts a second V_Cu⁻,
  giving the neutral (2V_Cu + In_Cu) complex.
- Kesterites: Sn_Zn²⁺ + Cu_Zn⁻ (net +1) binds a second Cu_Zn⁻, giving the neutral
  2Cu_Zn + Sn_Zn cluster.

## Acceptance criterion

Growth is accepted only if the added unit binds:

    ΔE_b = E_b(K+X) − E_b(K) − E_b(X) < 0,   E_b(X) = 0 for a single defect.

For a pair (K and X both single defects) this is simply E_b < 0. A negative
E_b(K+X) alone is **not** sufficient: in PPO, Cu_Pb2 + 2S_O1 has a negative E_b
over the upper part of the gap but stays above the bound pair Cu_Pb2 + S_O1, so
the second S does not bind.

`python scripts/codopex_tools.py deltaeb --eb-kx ... --eb-k ... [--eb-x ...]`

Evaluate ΔE_b across the gap, not only at E_F^eq, and report where it changes sign.

## Many-body terms

The pair table predicts E_b(K+X) as the sum of the binding energies of all
nearest-neighbour pairs in contact inside the complex. The difference between the
computed E_b and that sum is the many-body departure. It combines interactions
beyond the contact shell with genuine three-body (or higher) terms.

`python scripts/codopex_tools.py manybody --eb 0.00 --pairs -0.17 -0.17`

In PPO the departures were as large as the pair energies and of either sign:

| complex | contact pairs | pair sum | computed E_b | departure |
|---|---|---|---|---|
| Cu_Pb2 + 2S_O1 (both S on the Cu) | 2 × Cu–S (−0.17) | −0.34 | ≈ 0.00 | ≈ +0.34 (saturation: Cu already reduced) |
| 2Cu_Pb2 + S_O1 (S on one Cu) | Cu–Cu (0.94) + Cu–S (−0.17) | +0.77 | 0.38 / 0.54 | −0.39 / −0.23 |

Consequences for the procedure:

- Pair scores **rank** the arrangements of a complex; they cannot decide whether
  it binds. Every arrangement considered must be computed.
- Qualitatively, the binding energy of a three-defect complex was bracketed by
  those of its constituent pairs in all nine PPO compositions. Use this to rank,
  not to skip.
- The arrangement can matter as much as the composition (spread of 0.62 eV at
  E_F^eq among arrangements of one composition in PPO), so always enumerate all
  placements codopex returns for a composition.

## Practical checks

- **Interaction cutoff**: discard arrangements in which a sub-pair has separated
  beyond the interaction range (they are not complexes of the intended order).
- **Supercell size**: a complex must stay well separated from its periodic images.
  In PPO's 78-atom cell the bound pair lies within 9.2 Å of its own image, so a
  four-defect (pair + pair) complex could not be computed explicitly; the
  pair–pair check was made through the three-defect motifs it contains. State such
  assumptions in the report.
- **Charge-state search**: around q_c = Σ q_i of the constituents, including the
  composite's own q_c.
- Stop when no added unit binds.
