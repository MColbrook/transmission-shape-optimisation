# Computational inputs and their use

The [mathematical supplement](mathematical_supplement.pdf) gives the trial coefficients, contour and Hessian estimates, and the regular forms for the radial calibration checks.

The catalogue [`checks.json`](checks.json) gives the commands and their certificate paths.
The [figure directory](figures/) contains the plotting materials for Figure 1.
The inputs below are exact rational bounds used in the proofs.

## Global transition

The checker paths in this table are relative to `certificates/`.

| Checker | Inputs, range and error bound | Use in the proof |
| --- | --- | --- |
| `shape_optimization/certify_near_optimal_domain.py` | Twelve integer conformal coefficients, the clamped odd angular trial and $\tau=51963/10000$, $n_0=263447/50000$. The degree-600 reciprocal approximant $H$ has coefficient denominator $10^{24}$. If $R_*$ is the sum of absolute coefficients of $1-F'H$, then $R_* /(1-\eta)$ bounds the error using $\lVert1-F'\rVert_\infty\le\eta<1$. The checker computes $I_0$ and $I_1$ exactly and uses the reciprocal error bound only for $I_2$. | Negative trial form on the explicit smooth competitor. The contrast sweep and endpoint comparison reuse the moment bounds. |
| `continuation/certify_contrast_interval.py` | Closed interval $n\in[259/50,5281/1000]$. The form comparison uses one degree-five Bernstein interval. The order-$0,1,2$ determinants use degree-32 truncations and geometric tails from $28^k/(k!)^2$ and $k28^k/(k!)^2$ on $1\le n\le6,\ \lvert z\rvert\le8$. The 1,561 closed radial contour boxes and 29 angular intervals cover the required determinant and sign tests. A Bernstein bound gives positivity of $J_3$ on the positive squared-argument range up to $31$. | The determinant count and angular signs exclude disk eigenvalues of modulus at most $T(n)$ throughout the contrast interval. |
| `check_transition_endpoints.mjs` | Bounds on the squares of the first two $J_0$ zeros and the first $J_1$ zero, the plate root and disk roots at $5.281738,\ 5.4254619,\ 5.4254621$. For squared arguments below $31$, the ordinary Bessel series use the degree-30 and degree-31 alternating brackets. The modified series use the first omitted term times $100/99$. | The form at $n_b=5.281738$ and the slope comparison extend the certified interval to $n_b$. The checker bounds $t_D(n_b)<142$ and $5.4254619<n_{\rm pl}<5.4254621$. The plate bounds enter the product estimate for nonreal continuation. |
| `minimizer_literature/continuation/robin_denominators.py`, `global_variational/continuation/robin_normalization_bounds.py` and the `robin_roots.py` helper | Degree-32 Bernstein and degree-80 Bessel signs with geometric tails bound the Robin denominators on the two branches, using $[0,6]$ and $[14,31]$ for the outer bounds. The helper uses exact endpoint brackets, monotonicity and interval Newton on the two branches. | Normalisation and uniform derivative bounds for the radial profile in the one-signed comparison. |
| `minimizer_literature/continuation/all_contrast_center_reduced.py` | $(\eta,r^2)\in[0,1]\times[0,9/16]$, 27,167 rectangles. The checker uses degree-30 midpoint sums, the tail $\mathcal T_{30}$, the Poisson bound $1/m!$ and directed intervals for first and second derivatives. | Central chart, with the analytic cancellation of the singular powers at $r=0$. |
| `shape_optimization/continuation/all_contrast_endpoint.py A` | $s=\eta x$, $(\eta,x)\in[0,1/4]\times[0,1]$, 433 rectangles and degree 48. Exact divisions give finite enclosures at $x=0$ and $\eta=0$. The Taylor and derivative tail is bounded by $10^{30}E_C(48)(1/4)^{42}$. | Thin boundary chart for $s\le\eta$. |
| `shape_optimization/continuation/all_contrast_endpoint_fast.py E Esmall` | $s\in[0,1/4]$ and $\eta\in[1/4,1]$ or $[1/16,1/4]$, with 1,574 or 773 rectangles. Both use degree 38, exact division by the double zero and the independent exponential bound on the shortened $\mathcal W$ recurrence. The combined tail fits within $10^{22}E_C(38)(1/4)^{32}$. | Regular boundary charts, including $s=0$ and the limit $\eta=1$. |
| `shape_optimization/continuation/all_contrast_endpoint_B_covered.py` | $(s,\delta)\in[0,1/4]\times[0,1]$, $\eta=s\delta$. The driver uses `chart_B` in `all_contrast_endpoint.py` on 2,073 rectangles. It encloses the degree-38 fields separately, adds a common widening and bounds the transformed $\mathcal W$ tail separately. The 68 inherited rectangles and the portions of crossing rectangles above $\eta=1/16$ use the regular certificate. The direct calculation caps only the Robin interval at $1/16$. | Boundary chart for $\eta\le s$, including the resolved corner. |
| `continuation/audit_all_contrast_coverage.py` | Rational slab sweeps on all five closed regions, including the inherited part of B. The centre and boundary charts meet at $r=3/4$. | Together with the analytic tail bounds and chart changes, the checks give the strict signs required by the one-signed comparison at finite contrast. |
| `local_stability/check_exact_identities.py` | Seventeen exact fraction checks for disk coefficients and affine ellipse moments. Bessel signs use alternating or geometric tails. The verifier compares the transcript line by line. | Independent checks for the analytic disk and ellipse calculations. |

## Local stability below resonance

The paths are relative to `subresonant/`. The rational centre files supply
the starting boxes. Each checker verifies the closed parameter interval
in its catalogue entry, including the overlaps.

| Checker | Coverage and analytic completion |
| --- | --- |
| `unit_endpoint_check.py` | Degree-30 entire series with geometric tails at $n=1$. Isolates the first conjugate radial pair and bounds the endpoint Hessian inputs. |
| `unit_midrange/extension/check_roots.py`, `unit_midrange/extension/check_hessian.py` | $1\le n\le441/100$. Two Graeffe transforms and entire-series tails give the radial count. The 550 closed tubes control modes $2$ to $6$. The Bessel-product estimate covers higher modes. |
| `global_screen/plate_root_count_check.py`, `unit_midrange/cap/coarse_constants.py` | Rayleigh sums and Bessel zero bounds give the disk plate-product input for $\lvert\sqrt n\,\lambda\rvert\le14$. The root count and algebraic simplicity follow from the analytic argument. |
| `unit_midrange/cap/hessian_compact_check.py` | 808 closed tubes for $21/10\le\sqrt n\le459083/200000$, with entire-series derivative tails. The tubes bound the logarithmic-derivative quotients used for mode $2$. Higher modes are estimated analytically. |
| `unit_midrange/cap/hessian_endpoint_check.py` | $0\le\eta=(N-n)^{1/3}\le1/40$. The degree-12 Bessel ODE series cancels the exact $\eta^3$ factor before interval division, with a uniform determinant tail of order $\eta^{10}$. The value $\eta=0$ is the removable limit. The tube overlaps the compact check. |

The proofs use these spectrum and Hessian bounds for the subresonant
part of the classification. The nonlinear remainder and both
cube-root bounds follow from the analytic shape estimates.
