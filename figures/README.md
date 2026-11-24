# Figure 1

[Figure 1](candidate_geometry.pdf) shows the unit disk, the analytic competitor
$\Omega_{12}$ and a fitted domain with two corners. Each domain has area $\pi$.
At this scale the analytic and cornered boundaries nearly coincide. The bound
$\lambda_{\mathrm R}(\Omega_{12},N)<5.1963$ comes from the rational trial in
the [mathematical supplement](../mathematical_supplement.pdf). The domain with
corners is a numerical comparison.

## Generate the figure

From the repository root, run the following commands in a Python environment.
The plotting packages are needed only for this figure.

~~~sh
python -m pip install numpy scipy matplotlib
python figures/make_figure.py figure1-replay.pdf
~~~

The script uses Matplotlib's noninteractive backend. At 144 dpi the replay
with Python 3.11, NumPy 2.4.4, SciPy 1.17.1 and Matplotlib 3.10.9 matched
the supplied figure pixel for pixel. PDF metadata can vary between runs.

## Curves and parameters

The orange boundary is the image of the unit circle under the polynomial map
in the paper. The twelve integer numerators and their common denominator
are in [geometry_parameters.json](geometry_parameters.json). The area
normalisation uses the exact formula for the polynomial.

For the blue boundary the fitted map satisfies

$$
f(0)=0,\qquad
f'(z)=C(1+z^2)^\beta
\exp\left(\sum_{j=1}^8 p_j z^{2j}\right).
$$

The angle $\omega=\pi(1+\beta)$ was prescribed using the root of
$\tan\omega=\omega$ in $(\pi,3\pi/2)$. Its value is approximately
$257.45^\circ$. The eight fitted parameters $p_j$ came from a numerical
search for a low positive real transmission eigenvalue at resonance.
The figure uses the saved exponent and the eight coefficients at their
original floating-point precision. To evaluate $f$, the script multiplies
the binomial and exponential series, retains 20,000 coefficients of the
derivative and integrates term by term. If $g_j$ are these coefficients
before normalisation, the conformal area formula gives the numerical
normalisation $C=(\sum_{j=0}^{19999}|g_j|^2/(2j+1))^{-1/2}$.
Both image boundaries use 8,192 samples. The disk uses 4,097 samples
including the repeated endpoint.

The parameter file contains the exact polynomial coefficients and the
saved numerical fit used for the figure in the paper.
