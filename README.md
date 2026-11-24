# Shape optimisation of the least-modulus transmission eigenvalue

This repository accompanies the paper of the same title by Matthew J. Colbrook and Catherine Drysdale. The exact calculations bound the trial on the analytic competitor, exclude disk eigenvalues on the comparison interval, establish the subresonant Hessian estimates and verify the signs in the radial calibration.

The [mathematical supplement](mathematical_supplement.pdf) gives the trial data and the analytic error bounds. Its [LaTeX source](mathematical_supplement.tex) is included. The [certificate index](checker-map.md) gives the parameter ranges, tails and role of each calculation, and [checks.json](checks.json) specifies the commands and their inputs.

## Verification

Run from the repository root with Python 3.10 or later. The quick and full modes also require Node.js on PATH. The exact checks use the standard libraries of these runtimes.

~~~sh
python verify.py quick --output-dir verification-quick
python verify.py subresonant --output-dir verification-subresonant
python verify.py full --output-dir verification-full
~~~

| Mode | Commands | Coverage |
| --- | ---: | --- |
| `quick` | 12 | Trial, contrast interval, transition endpoints, selected subresonant estimates and coverage of the supplied radial charts |
| `subresonant` | 7 | All checks from unit contrast to the resonant endpoint |
| `full` | 18 | All checks, with fresh generation of the radial chart partitions before the coverage test |

The verifier runs the scripts in a temporary copy and compares the calculated bounds and partitions with the supplied certificates. Each `--output-dir` must be a new directory. Omit this argument to discard the detailed logs when the run finishes. Use `--jobs` to set the number of concurrent checks.

## Files

- `certificates/` contains the trial, contrast and radial calibration calculations, their exact certificates and the common interval and Bessel routines.
- `subresonant/` contains the root, plate and Hessian calculations and their rational data.
- [`figures/`](figures/) contains the figure in the paper and its plotting inputs. The numerical corner fit is included here as saved data.

## Licence

The verification and plotting code, its documentation, and the accompanying inputs and certificates are available under the [MIT License](LICENSE). The authors reserve rights in the mathematical supplement (LaTeX source and PDF) and the supplied Figure 1 PDF.
