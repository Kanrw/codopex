# Scoring and ranking pairs (step 4)

The scores decide the **order** in which pairs are computed. They are hypotheses
about which interactions act; the DFT binding energy decides whether a pair binds.

## Two kinds of mark

| Mark | Judged from | Question it answers |
|---|---|---|
| **+** chemical environment | chemical identity and charge states q_i of the two dopants at E_F^eq | Is there a mechanism that would make this pair bind? |
| **−** local environment | the surrounding host lattice | Does the lattice meet the precondition that mechanism needs? |

A − is always attached to a specific +. A − with no matching + has no meaning
(the pair already has nothing in its favour) and is ignored.

### + mechanisms (extend the list when the host calls for it)

- **Charge compensation**: a donor and an acceptor at E_F^eq (opposite q_i).
- **HSAB matching**: soft acid next to soft base, or hard next to hard
  (e.g. Cu⁺, Ag⁺ with S²⁻, Se²⁻; hard: high-valent small cations with O²⁻, F⁻).
  Use the oxidation state the dopant will have in the complex if it is known,
  and say which state was assumed.
- **Hole- (or electron-) mediated coupling**: two dopants that each introduce a
  hole (electron) whose states could hybridize into a bonding combination.

Other mechanisms the user may name (strain compensation, dangling-bond
passivation, spin coupling) are welcome as long as each is written down as a
stated hypothesis with the reason it applies.

### − unmet preconditions

- **Polyanion pull**: in a polyanionic host, rigid neighbouring polyanions with a
  larger charge than the dopant partner draw the dopant away from it
  (compare the effective or nominal charges of the units around the dopant).
- **Restricted relaxation**: the dopants cannot move toward each other because the
  surrounding units are stiff (compare the dopant displacement with the internal
  distortion of the neighbouring units in the relaxed single-defect cell).
- **Insufficient proximity**: the sublattice keeps the two dopants too far apart
  for their states to overlap (nearest-neighbour separation from codopex
  `dist_AB`, compared with the extent of the defect states).

Do **not** use "like-charge repulsion" as a − mark. Electrostatic repulsion of
two like-charged defects is what remains *after* a mechanism has failed, not a
precondition, and hole-mediated coupling is precisely a case where like-charged
defects may attract.

## Tiers

| Tier | Rule | Meaning |
|---|---|---|
| I | at least one + not offset by any − | compute first |
| II | every + offset by a − | compute next |
| III | no + | compute if the budget allows |

`python scripts/codopex_tools.py tiers marks.json` builds the table; an example
input for Pb3(PO4)2 is `templates/marks_ppo.json`.

## Reliability

- In Pb3(PO4)2 every tier-III pair was unbound at E_F^eq, and the only bound pair
  was in tier I. That is one host. In a new host, compute at least a sample of
  tier III before trusting it.
- The local-environment marks are judged on **isolated** defects. Complex
  formation can change oxidation states, sizes and charges (in PPO, Cu²⁺ → Cu⁺ on
  binding), so a − can be wrong. That is why a − lowers priority and never deletes
  a pair.
- Present the scoring to the user as a table with one line of justification per
  mark before any calculation is launched.
