# Standalone scope verification
Compile all C++ files with g++ -O2 -std=c++17. Python scripts use the imports in their sources. Run from this replay directory.
- ./split parts > parts.txt and ./moments > moments.txt. python fullcert.py verifies the exact force/matrix sums; rank.cpp and rank6.cpp check the modular rank certificates.
- exact8.cpp checks the exact squared minimum and 48,600 contacts per site. ./group > group.txt and python identify.py check the attribution-preserving equivalence to the two-bank Leech construction.
- To replay the fixed original metric seven-packet obstruction for all period advances: ./braid 12 1 63 any > coarse.out; extract the positive unresolved period pairs from the OPEN lines, then ./braid 24 1 63 filter < coarse.txt. The second pass has no survivors.
- To replay the fixed recovered-record metric's discrete cases use ./fine 12 1 60708 any 9 base 8 0 1 192 (for 10 and 12 replace both the threshold and multiplicity by 60154,10 and 59208,12 respectively). Where needed, include the corresponding grouped period table whose generator is group.cpp. Status found for rad=1 can signify an unresolved conservative box, not a configuration.
- ./slice > kcost.txt then ./three 12 24 1 0 proves the separate continuous third-coordinate seven-carrier exclusion for the original fixed advances.
Only exact proof and necessary input files are included. No exploratory raw trial output or unresolved nine-site continuous box search is included. Neither equality with the known record nor any fixed-family negative claim improves a global sphere packing bound.
