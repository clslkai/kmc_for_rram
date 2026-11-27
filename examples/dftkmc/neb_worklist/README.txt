DFT/NEB WORKLIST — UNRELAXED ENDPOINTS

Each task contains POSCAR.initial, POSCAR.final and mapping.json. Atom identities and ordering are identical in both endpoint files. The target oxygen fills the source vacancy, moving by minus the actual periodic vacancy-hop displacement. Cartesian positions outside the cell preserve the intended path image.

These are candidate structures, not completed DFT or NEB calculations. Relax both endpoints at the same fixed cell and charge, verify metastability and atom mapping, then construct and converge NEB images. The conventional-cell geometry and supercell size also need convergence/relaxation assessment.

Set NELECT from the valence electron count of your own POTCAR minus q. Positive q removes electrons. Use identical POTCAR, atom count, charge, spin treatment and numerical settings throughout each path. No POTCAR is distributed and no calculation has been submitted.

catalog.pending.json deliberately contains missing numerical data and cannot run. Fill it only with completed, traceable calculations on the same structure and charge branches. Site and saddle energies require a common reference. An activation barrier alone is not an absolute saddle energy; saddle = initial site energy + forward barrier, and verify the reverse relation. Attempt prefactors require a declared calculation or approximation.

The channel list removes translation and reverse duplicates only; it is not a complete space-group symmetry reduction. A requested subset exports only those tasks, while the pending catalogue retains the complete channel worklist so uncovered rates remain blocked.

Source catalogue channel exclusions are preserved as provenance only and were NOT automatically applied. The pending catalogue regenerates all geometric candidate channels within the export cutoff; an explicit channel selection restricts endpoint tasks only. Source IDs/reasons and catalogue identity are retained in worklist.json and catalog.pending.json metadata.
