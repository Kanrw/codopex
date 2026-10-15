---
name: codopex-workflow
description: Guide a user through the physically guided co-doping workflow of codopex (Wang, Yao, He, Zhao) to find which defect complexes bind in a host crystal and how abundant they are. Use when someone wants to set up or run a co-doping / defect-complex search, enumerate symmetry-inequivalent dopant pairs or triples, decide which pairs to compute first, judge binding from DFT formation energies, expand bound pairs into larger complexes, or check whether complexes change the equilibrium Fermi level. Works with doped for point defects and the Fermi level, codopex for enumeration, and VASP or another DFT code for the energies.
---

# codopex workflow

This skill walks a user from "I want to co-dope host H with dopants X and Y" to a
list of bound defect complexes with their abundance. It does not replace any DFT
calculation: the rules below only decide **which** configurations are computed and
in what order. Every conclusion about binding must rest on computed energies.

The workflow has two parts (Fig. 5 of the paper):

- **(a) Pair prioritization**: point defects → representative defects →
  nearest-neighbour pairs → score each pair → rank into priority tiers.
- **(b) Evaluation and expansion**: charge-state search → binding energy
  E_b(E_F) → expand while the added unit binds → abundance and E_F^eq check.

Read `references/` only when the step needs it:

| File | Read when |
|---|---|
| `references/scoring.md` | assigning + / − marks and tiers (step 4) |
| `references/expansion.md` | growing complexes, ΔE_b, many-body terms (step 6) |
| `references/fermi_level.md` | choosing E_F^eq, abundance, self-consistency (steps 1, 7) |
| `references/ppo_example.md` | the user wants a worked example, or to check a result against the paper |

`scripts/codopex_tools.py` holds the small calculations (binding-energy curves,
ΔE_b, many-body departure, tiers, abundance, Fermi-level shift). It needs only
the Python standard library. Run `python scripts/codopex_tools.py -h`.
`templates/system.yaml` is the system description to fill in with the user;
`templates/marks_ppo.json` is the scored pair table of the paper's example.

## Working with the user

1. Start by filling `templates/system.yaml` with the user: host structure file,
   dopants, growth temperature, chemical-potential limits of interest, supercell,
   and the compute budget (how many complexes they can afford). Ask only for what
   is missing.
2. At each step, show the user what was decided and why before moving on. The
   scoring step in particular must be presented as a table with one line of
   justification per mark, so the user can overrule it.
3. Never invent DFT numbers. If an energy is not yet computed, say so and stop at
   the step that needs it.
4. Keep the user's terminology for their system; use the labels codopex produces
   (`{dopant}_{host}_{site symmetry}`, e.g. `Cu_Pb_D3d`) so files and tables match.

## Procedure

### Step 1. Point defects and E_F^eq (doped)

Run standard point-defect calculations for every inequivalent site and charge
state, including intrinsic defects, with doped:

```python
from doped.generation import DefectsGenerator
from doped.chemical_potentials import CompetingPhases, CompetingPhasesAnalyzer
from doped.analysis import DefectsParser

gen = DefectsGenerator("POSCAR_host", extrinsic={"Pb": "Cu", "O": "S", "P": "S"})
cp_gen = CompetingPhases("Pb3(PO4)2", extrinsic=["Cu", "S"])
# ... write inputs (doped.vasp.DefectsSet), run DFT, then:
cpa = CompetingPhasesAnalyzer("Pb3(PO4)2", entries="competing_phases/")
cpa.calculate_chempots(extrinsic_species=["Cu", "S"])
dp = DefectsParser("defect_calcs/", dielectric=eps)   # finite-size corrections
thermo = dp.get_defect_thermodynamics(chempots=cpa.chempots)
EF_eq = thermo.get_equilibrium_fermi_level(bulk_dos="vasprun_dos.xml",
                                           limit="A", temperature=900)
table = thermo.get_formation_energies(limit="A", fermi_level=EF_eq)
```

Record, for each defect: charge state q_i at E_F^eq and formation energy ΔH_i.
E_F^eq is solved with the isolated defects only; its validity is checked in step 7.
Use the growth (annealing) temperature, not room temperature, for E_F^eq, because
that is where the complexes form. See `references/fermi_level.md`.

API names above follow doped 3.2; check them against the installed version.

### Step 2. Representative defects

Choose sublattices on **chemical grounds** (electronegativity, ionic size): a
cation dopant on cation sites, an anion dopant on anion or polyanion-centre sites.
Keep **all** symmetry-inequivalent sites of those sublattices. Do not cut by
formation energy: in the paper the only bound pair involves the Cu site with the
*higher* isolated formation energy. Formation energies enter only in step 7.

Interstitials and vacancies are legitimate members of the set. codopex enumerates
substitutions only; generate interstitial sites separately (doped's
`DefectsGenerator` uses Voronoi tessellation) and treat them by hand in steps 3–6.
Tell the user if the defect set is restricted, and why.

### Step 3. Enumerate pairs (codopex)

```python
import codopex as cp
pairs = cp.generate_pairs("SPOSCAR", [("Pb", "Cu"), ("O", "S"), ("P", "S")])
pairs.type_labels                      # representative defect types
rows = cp.pair_rows(pairs)             # combo, shell, dist_AB, degeneracy
cp.io.export_pairs(pairs, "runs/pairs")  # POSCAR tree for DFT
cp.io.save(pairs, "runs/pairs.json.gz")  # reload later with cp.io.load
```

Command line: `codopex pairs SPOSCAR -r Cu@Pb -r S@O -r S@P --out runs/pairs`.
Use `shells="nnn"` when longer-ranged interactions are expected (charged dopants
in a weakly screening host). The input must be the supercell that will be used
for DFT; codopex does not build or check it.

### Step 4. Score and rank (read `references/scoring.md`)

For each pair:

- **+ (chemical environment)**: one per mechanism that the chemical identity and
  charge states of the two dopants favour: charge compensation, HSAB matching,
  hole-mediated coupling, or any other mechanism the user can name and justify.
- **− (local environment)**: one per precondition of such a mechanism that the
  surrounding lattice does not meet: competing electrostatic pull of neighbouring
  polyanions, restricted relaxation, insufficient proximity of the two dopants.
  A − only counts against an existing +.

Tiers: **I** = a + not offset by any −; **II** = every + offset by a −;
**III** = no +. Compute tier I first, then II, then III as the budget allows.
Local-environment marks only *lower priority*; they never delete a pair, because
they are judged on isolated defects and can be wrong once the complex forms.

```bash
python scripts/codopex_tools.py tiers marks.json      # marks.json: see -h
```

### Step 5. Compute in priority order (DFT)

For each pair taken in tier order:

1. Relax the pair supercell in each charge state around q_c = Σ q_i
   (at least q_c − 1, q_c, q_c + 1). The lowest state of the complex is often not
   the sum of its constituents.
2. Apply the same finite-size correction scheme as for the point defects. For
   eFNV on a complex, place the point charge and the sampling centre at the
   centroid of the constituent defects (pydefect).
3. Build E_b(E_F) from the formation-energy lines and report it across the gap,
   not only at E_F^eq:

```bash
python scripts/codopex_tools.py eb eb_input.json --ef-eq 1.22
```

A pair is bound where E_b < 0. Report the charge state of the complex and any
change of oxidation state (magnetic moment, Bader charges), since it often
explains the binding and matters for the next step.

### Step 6. Expand bound complexes (read `references/expansion.md`)

Only a bound complex K is expanded. Add a unit X, either one more defect (near A,
near B, or near both) or a second bound complex, and treat K as a composite defect
with its own charge state and chemical identity.

```python
triples = cp.generate_triples(pairs, pair="Cu_Pb_D3d+S_O_Cs", third=[("O", "S")])
cp.triple_rows(triples)   # criterion nearA_farB / nearB_farA / both_nn, d_AC, d_BC
cp.io.export_triples(triples, "runs/triples")
```

Score the **new contacts** as in step 4 to rank the arrangements, compute them,
and accept the growth only if the added unit binds:

    ΔE_b = E_b(K+X) − E_b(K) − E_b(X) < 0       (E_b(X) = 0 for a single defect)

A negative E_b(K+X) alone is not enough. Compare every computed complex with the
sum of its contact-pair energies (`codopex_tools.py manybody`) and report the
departure: it is the many-body part that pair scores cannot predict. Check that
the supercell still separates the complex from its periodic images; complex–
complex combinations usually need a larger cell. Stop when no added unit binds.

### Step 7. Abundance and E_F^eq check (read `references/fermi_level.md`)

For each converged complex: ΔH_c = Σ_i ΔH_i + E_b, then its concentration in the
dilute limit. Estimate both limits that bracket the experiment:

```bash
python scripts/codopex_tools.py abundance --eb -0.17 --z 6 --T 900 --dh-partner 1.01
python scripts/codopex_tools.py abundance --eb -0.17 --z 6 --T 900 --x-partner 0.01
```

Then ask whether the complex can move E_F^eq: compare |q_c| c_c with the largest
charged term in the neutrality condition (`codopex_tools.py fermishift`). If the
ratio r is above ~0.1, re-solve E_F^eq with the complexes included (doped
`FermiSolver`, e.g. `fixed_defects={...}` with the py-sc-fermi backend for a
fixed complex concentration) and repeat steps 4–7 at the new E_F^eq.

### Step 8. Report

Give the user: the defect set and its restrictions; the scored pair table with
tiers; which pairs were computed; E_b(E_F) of every computed configuration with
the bound ones marked; the expansion record with ΔE_b and many-body departures;
abundance in both limits; and whether E_F^eq is self-consistent. State every
assumption (supercell size, cutoffs, which tiers were skipped).

## Rules that keep the answer honest

- The scores decide the order of calculation, not the answer.
- Skipping a tier is a budget decision; say so explicitly and name what was skipped.
- A pair unselected by every mechanism was unbound in the paper's host, but that is
  an observation for one host, not a guarantee for the next.
- When a mechanism's prediction fails, report it as a failed hypothesis (which
  precondition was not met) rather than tuning the rules to fit.
