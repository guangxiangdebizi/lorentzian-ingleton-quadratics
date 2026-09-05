# Building the note

The source uses standard LaTeX packages. With `latexmk` and BibTeX installed:

```bash
cd paper
latexmk -pdf main.tex
```

The generated PDF is intentionally not tracked. The mathematical statement
does not depend on the finite computations in `src/`; those computations are
regression certificates for formulas and conventions.

The 2026-09-05 revision was also compiled with the workspace's Tectonic 0.17.0,
without network access (cached TeX bundle only):

```powershell
& '../../.tmp-tools/tectonic/bin/tectonic.exe' --only-cached --keep-logs --keep-intermediates main.tex
```

The main proof contains the full explicit `476656^(-5)` construction. The
earlier semialgebraic argument remains an alternative qualitative proof;
neither argument relies on an unpublished external note.
