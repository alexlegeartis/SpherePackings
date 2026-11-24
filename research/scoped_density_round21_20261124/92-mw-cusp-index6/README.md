# MW92: direct sixth-jet short-witness screen

At least **2,897,035,069 of 11,013,546,097** translation-invariant
index-47^6 sublattices do not beat the parent. This consists of an
exact finite count of 2,892,049,308 primitive profiles with a literal
witness of squared norm at most 26, plus all 4,985,761 mixed
profiles by the [complete fifth-jet theorem](../92-mw-complete-index5/NONINVARIANT.md).
It does not claim that the unhit profiles attain any particular
minimum, or that the screened witnesses exhaust their shells.

For the primitive profiles use the same normalized projective
ratio over F_47[D]/D^6 as in Round 20. The replay checks rank 12
of the compiled quotient and the six-digit shift induced by D.
For each retained row it computes the integral Gram norm and
checks depth at least two. Its quotient ratio at depth d supplies
a prefix with 47^d lifts. Closing prefixes under the exact
normalizer, and retaining shortest ones first, gives the
following disjoint contributions:

| Witness depth | Shortest prefixes | Profile count |
| ---: | ---: | ---: |
| 4 | 552 | 2,693,583,912 |
| 3 | 552 | 57,310,296 |
| 2 | 63,900 | 141,155,100 |

The seed file retains one witness per contributing normalizer
orbit rather than the full search archive. It reuses the previous
short rows, the new terminal affine-CVP witnesses, and some
additional deep eigensublattice rows found during search. Every
published row's acceptance is exact and independent of discovery.
The comparison 26^46 < 16^46·47^6 is verified by integer arithmetic.

Run OPENBLAS_NUM_THREADS=1 python -B verify.py (NumPy);
the optional explicit arguments are the Round 19 index-three
verify.py, the Round 17 verify.py, and the Round 18 index-two
directory in that order.
