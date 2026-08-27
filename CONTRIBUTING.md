# Contributing

The most valuable contributions are proof audits, prior-art matches, and exact
counterexamples.

When reporting a mathematical issue, please identify the earliest failing
step and include one of:

- an explicit polynomial and four sets;
- a valuated rank-two table with the pair order fixed;
- a precise citation and theorem number;
- a minimal correction to the module-length or curve-selection argument.

Do not report a floating-point sample alone as a counterexample. Include exact
rational data or a certified interval whenever possible.

Code changes should preserve the standard-library-only reproduction path:

```bash
python src/reproduce.py --check reports/certificates.json
python -m unittest discover -s tests -v
```
