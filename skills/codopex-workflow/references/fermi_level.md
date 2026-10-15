# Equilibrium Fermi level, abundance and self-consistency (steps 1 and 7)

## Where E_F^eq enters

1. **Input.** E_F^eq fixes the charge state q_i of every isolated defect, which the
   scoring step uses (donor, acceptor, neutral). Solve it with doped at the growth
   or annealing temperature for each chemical-potential limit the user cares about:

   ```python
   EF_eq = thermo.get_equilibrium_fermi_level(bulk_dos="vasprun.xml", limit="A",
                                              temperature=900)
   ```

   If the material is quenched, doped's
   `get_fermi_level_and_concentrations(annealing_temperature=..., quenched_temperature=...)`
   gives the frozen-defect Fermi level at the operating temperature. Binding is
   decided at the growth temperature, where the complexes form.
2. **Check.** E_F^eq was solved with the isolated defects only. Once bound
   complexes are found, check that including them does not move it.

Report all binding energies across the gap, so that no conclusion rests on the
single value of E_F^eq.

## When does a complex move E_F?

Charge neutrality: p − n + Σ_i q_i c_i = 0. A complex enters only through its
charge density |q_c| c_c, and matters only when this is comparable to the largest
charged term (dominant charged defect, or free carriers). Define

    r = |q_c| c_c / (largest charged term).

For E_F pinned by a compensating donor D⁺ and acceptor A⁻ (typical of wide-gap
insulators), adding a singly charged species with fixed concentration r·c_A gives

    |ΔE_F| = k_B T ln[(r + sqrt(r² + 4)) / 2]  ≈ r k_B T / 2   (r ≪ 1).

If its concentration follows E_F like the acceptor's: |ΔE_F| = (k_B T / 2) ln(1 + r).
For a Fermi level set by free carriers n: |ΔE_F| ≈ −k_B T ln(1 − r), r = |q| c / n.

At 900 K (k_B T = 0.0776 eV):

| r | fixed c | c follows E_F |
|---|---|---|
| 1e-4 | 0.004 meV | 0.004 meV |
| 0.1 | 3.9 meV | 3.7 meV |
| 1 | 37 meV | 27 meV |
| 3.4 | 0.10 eV | 57 meV |
| 10 | 0.18 eV | 93 meV |

`python scripts/codopex_tools.py fermishift --r 1e-4 0.1 1 10 --T 900`

Rule of thumb: a shift of order k_B T needs r ≈ 1; a shift above 0.1 eV needs the
complex to exceed the dominant charged term by several times. The dependence is
logarithmic, so E_F is robust.

A simple bound often suffices: if a charged constituent (here Cu_Pb1⁻) holds most
of a dopant, the largest negative term is at least that constituent's charge
density, so r ≤ c_complex / c_constituent, without knowing which defect dominates.

## Situations in which complexes do move E_F

1. The complex becomes the main form of a dopant (large |E_b| and high doping).
2. The complex's charge differs from the sum of its constituents' charges
   (PPO: two neutral defects give a q = −1 complex, i.e. a new acceptor).
3. Passivation: a donor and an acceptor that pin E_F bind into a neutral complex,
   removing both terms (Mg–H in GaN).
4. The constituents are themselves the dominant charged defects.
5. Quenching or non-equilibrium growth freezes in concentrations formed at high T.

## Abundance

    ΔH_c = Σ_i ΔH_i + E_b

Dilute-limit estimates, both worth giving because they bracket experiment:

- **Open system** (chemical potentials fixed): fraction of defect A bound to
  partner B = z exp[−(ΔH_B + E_b)/k_B T], z = number of distinct B sites adjacent
  to A.
- **Closed system** (dopant content fixed, B occupying a fraction x_B of its sites):
  c_AB / c_A,free = z x_B exp(−E_b/k_B T).
- If A can sit on several host sites, multiply by the site partition
  c(A on the complex-forming site) / c(A on the preferred site)
  = (site ratio) exp(−Δ/k_B T).

`python scripts/codopex_tools.py abundance --eb -0.17 --z 6 --T 900 --dh-partner 1.01 --x-partner 0.01 --site-gap 0.48 --site-ratio 0.5`

## Re-solving E_F^eq with complexes (doped)

If r > ~0.1, include the complex in the neutrality condition and solve again.
With doped, add the complex as an additional defect entry of the
`DefectThermodynamics` object (formation energies from step 5) and re-solve.
For a closed-system check, also fix its total concentration with the
py-sc-fermi backend (the name must match the entry you added):

```python
from doped.thermodynamics import FermiSolver
fs = FermiSolver(thermo, bulk_dos="vasprun.xml", backend="py-sc-fermi")
df = fs.scan_temperature(temperature_range=900, limit="A",
                         fixed_defects={"Cu_Pb_D3d+S_O_Cs": c_complex_cm3})
```

Check the argument names against the installed doped version. Then repeat the
scoring and expansion at the new E_F^eq and report both values.
