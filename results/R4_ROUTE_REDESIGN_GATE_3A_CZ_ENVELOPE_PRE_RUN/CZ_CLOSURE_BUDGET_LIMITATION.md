# Supplementary fixture-specific structural budget limitation

This is a post-policy symbolic diagnosis of the two already registered fixtures. It changes no fixture, grid, bound, cutoff, tau or ladder. It performs no acoustic/forward/feature calculation.

The original correction policy met the requested **single unresolved-path** necessary lower bound, 37/45. A domain-midpoint reference belongs to both closed children at each axis first split, so its tree has 16 distinct terminal reference-containing boxes. With all other branches ideally rejected immediately, complete partition requests are:

| Fixed fixture | COARSE ideal necessary requests | FINE ideal necessary requests |
| --- | ---: | ---: |
| CORNER_LOW | 37 | 45 |
| DOMAIN_MIDPOINT | 451 | 579 |

The midpoint counts are derived by recursively following **every** closed child containing the prescribed point. COARSE has 225 internal retained-path boxes, 210 ideally excluded sibling boxes and 16 terminal reference boxes: 451 requests. FINE has 289 internal boxes, 274 excluded siblings and 16 terminal reference boxes: 579. Exact intermediate counts and derivation are in CZ_FIXTURE_GEOMETRIC_BUDGET_ANALYSIS.json.

The registered maxima are 439/879. Thus the coarse maximum satisfies the single-path minimum but is **12 requests below the midpoint-specific ideal minimum**. No perfectly efficient exclusion mechanism could close that midpoint's full terminal partition within 439 under these closed-cell/traversal rules. The fine maximum satisfies its optimistic geometric minimum; that is not a promise that the physical interval bounds prune all other branches sufficiently.

This restriction must be included in the final negative PRE-RUN result. It must not be presented as CZ having no information, a physics impossibility, or failure despite an adequate coarse midpoint budget. The earlier policy was insufficient for this multi-branch fixture; no successful budget can be selected. The frozen ladder still runs exactly as registered and no C4/F4 or post-result budget/fixture change is made. Development remains disabled.
