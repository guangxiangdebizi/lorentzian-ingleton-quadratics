# Building the note

The source uses standard LaTeX packages. With `latexmk` and BibTeX installed:

```bash
cd paper
latexmk -pdf main.tex
```

The generated PDF is intentionally not tracked. The mathematical statement
does not depend on the finite computations in `src/`; those computations are
regression certificates for formulas and conventions.
