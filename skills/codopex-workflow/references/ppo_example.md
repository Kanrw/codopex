# Worked example: Cu–S co-doping of Pb3(PO4)2 (the paper)

Use this to show the user what each step produces, or to check that a rerun
reproduces the published numbers. Values: PBE+U (U_eff = 5 eV on Cu 3d), 78-atom
supercell, eFNV corrections (centroid convention for complexes), condition A,
T = 900 K, E_F^eq = 1.22 eV.

## Steps 1–2

R-3m host. Sites: Pb1 (6c, C3v), Pb2 (3a, D3d), O1 (18h, Cs), O2 (6c, C3v), P (6c, C3v).
Representative defects (all inequivalent sites of the chemically allowed
sublattices; codopex labels in brackets):

| defect | codopex label | q at E_F^eq |
|---|---|---|
| Cu_Pb1 | Cu_Pb_C3v | −1 (acceptor) |
| Cu_Pb2 | Cu_Pb_D3d | 0 |
| S_O1 | S_O_Cs | 0 |
| S_O2 | S_O_C3v | 0 |
| S_P | S_P_C3v | +1 (donor) |

Cu_Pb1 is 0.48 eV more stable than Cu_Pb2, yet only Cu_Pb2 binds S: do not cut the
defect set by formation energy.

## Steps 3–5: the fifteen nearest-neighbour pairs

| pair | + | − | tier | E_b at E_F^eq (eV) | q |
|---|---|---|---|---|---|
| Cu_Pb1 + S_O1 | HSAB | | I | 0.27 | −1 |
| Cu_Pb1 + S_O2 | HSAB | | I | 0.18 | −1 |
| Cu_Pb2 + S_O1 | HSAB | | I | **−0.17** | −1 |
| Cu_Pb2 + S_O2 | HSAB | | I | 0.47 | 0 |
| Cu_Pb1 + Cu_Pb1 | hole coupling | insufficient proximity | II | 1.05 | −2 |
| Cu_Pb1 + Cu_Pb2 | hole coupling | insufficient proximity | II | 0.82 | −1 |
| Cu_Pb2 + Cu_Pb2 | hole coupling | insufficient proximity | II | 0.94 | 0 |
| Cu_Pb1 + S_P | charge compensation | polyanion pull | II | 0.35 | 0 |
| Cu_Pb2 + S_P | | | III | 0.46 | +1 |
| S_O1 + S_O1 | | | III | 0.57 | 0 |
| S_O1 + S_O2 | | | III | 0.81 | 0 |
| S_O2 + S_O2 | | | III | 0.77 | 0 |
| S_O1 + S_P | | | III | 0.53 | +1 |
| S_O2 + S_P | | | III | 0.61 | +1 |
| S_P + S_P | | | III | 0.73 | +2 |

Machine-readable marks: `templates/marks_ppo.json`.

What the failures taught:

- **Charge compensation** (Cu_Pb1 + S_P): positive over the whole gap; on relaxation
  the dopants move apart (3.37 → 4.37 Å), because the rigid PO4³⁻ units pull Cu
  more strongly than the SO4²⁻ unit.
- **Hole-mediated coupling** (Cu–Cu): positive at every Fermi level; Cu 3d bands
  unchanged and DOS simply doubled (no hybridization). The Pb sublattice keeps the
  Cu atoms 3.76–5.56 Å apart (≥ 3.67 Å after relaxation).
- **HSAB** (Cu_Pb2 + S_O1): bound, E_b = −0.17 eV at E_F^eq, down to −0.45 eV at
  higher E_F; Cu–S bond 2.21 Å; the complex is q = −1 while both constituents are
  neutral, i.e. Cu²⁺ → Cu⁺ on binding.

## Step 6: expansion

| complex | arrangement | E_b at E_F^eq (eV) | ΔE_b vs bound pair |
|---|---|---|---|
| Cu_Pb2 + 2S_O1 | 2nd S on the Cu, away from the 1st S | ≈ 0.00 | +0.17 (does not bind) |
| Cu_Pb2 + 2S_O1 | two S adjacent | ≈ 0.37 | +0.54 |
| 2Cu_Pb2 + S_O1 | S binds one Cu (Cu–Cu 5.56 Å) | 0.38, 0.54 | +0.55, +0.71 |

Pair + pair (two Cu_Pb2 + S_O1): every contact reproduces one of the unfavourable
motifs above and both pairs have q = −1, so they do not associate; an explicit
four-defect calculation would need a larger cell. The expansion stops at the pair.

Many-body departures from the pair sum: +0.34 eV (Cu_Pb2 + 2S_O1) and −0.23 to
−0.39 eV (2Cu_Pb2 + S_O1).

## Step 7: abundance and E_F^eq

- Open system: bound fraction of Cu_Pb2 ≈ 6 exp[−(1.01 − 0.17)/k_B T] ≈ 1e-4.
  With Cu_Pb2/Cu_Pb1 ≈ 0.5 exp(−0.48/k_B T) ≈ 1e-3, the complex is ≈ 1e-7 of the Cu.
- Closed system with S on a fraction x_S of O1 sites: bound/free Cu_Pb2 ≈ 54 x_S,
  i.e. about one third of Cu_Pb2 bound at x_S = 1 %, but still ≈ 5e-4 of Cu_Pb1⁻.
- Cu_Pb1⁻ holds almost all of the Cu, so r ≤ 1e-7 (open) or ≈ 5e-4 (closed, 1 %):
  E_F^eq is self-consistent (direct re-solve: shift < 1e-4 eV).
- Conclusion: Cu_Pb2 + S_O1 is the only configuration in which Cu and S associate,
  but a minority form of Cu in equilibrium.
