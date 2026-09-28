# Manuscript reference snapshots

`verify_numbers.py` compares Tables V and SI of the paper with the tables generated from the result files.
When the LaTeX sources are available it reads them directly (set `PAPER_DIR=/path/to/paper`, or place them in
`paper/` next to or one level above this repository); otherwise it uses these snapshots of the two tables as printed
in the submitted manuscript. Every other number the checker compares is written into `verify_numbers.py` itself,
with the section of the paper it comes from.
