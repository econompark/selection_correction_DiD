// ------------------------------------------------------------- description ---
* In this example data, the parallel trends assumption is violated due to unobserved factors.

* The true average treatment effects are equal to the relative time. For instance, (g, t) = (2, 5) has a relative time of 3 (i.e., t - g = 3), so ATT(2, 5) = 3.

* The following example demonstrates that the method proposed by Park (2025) successfully estimates the treatment effects under this setup, while the method by Callaway and Sant'Anna (2021), which relies on the conditional parallel trends assumption, does not.

* The dataset is simulated data with the same structure as the mpdta dataset by Fernando Rios-Avila (friosavila.github.io).

// -------------------------------------------------------------------- data ---
use example_data, clear

// ----------------------------------------------------------- certification ---
// ordered probit
recode first_treat (2 = 1) (4 = 2) (5 = 3) (0 = 4), gen(G)

qui oprobit G X
mat b = e(b)
matrix list b

// generalized inverse Mills' ratio
predict xpi, xb
gen a1xpi = b[1, 2] - xpi
gen a2xpi = b[1, 3] - xpi
gen a3xpi = b[1, 4] - xpi

forv l = 1/3 {
	gen pdf`l' = normalden(a`l'xpi)
	gen cdf`l' = normal(a`l'xpi)
}

gen lam_1 = -(pdf1 / cdf1)
gen lam_2 = (pdf1 - pdf2) / (cdf2 - cdf1)
gen lam_3 = (pdf2 - pdf3) / (cdf3 - cdf2)
gen lam_4 = (pdf3) / (1 - cdf3)

gen lam = .
forv l = 1/4 {
	replace lam = lam_`l' if G == `l'
}
drop lam_*

// trend variable
xtset pid tid
gen Y21 = Y - L.Y  if tid == 2
gen Y31 = Y - L2.Y if tid == 3
gen Y41 = Y - L3.Y if tid == 4
gen Y51 = Y - L4.Y if tid == 5

gen Y43 = Y - L.Y  if tid == 4
gen Y53 = Y - L2.Y if tid == 5

gen Y54 = Y - L.Y  if tid == 5

// estimate: (g, t) = (2, 2) & true effect = relative time
qui reg Y21 X lam if first_treat == 0 | first_treat > 2, nocons
qui predict Y21hat, xb
qui gen teff22 = Y21 - Y21hat
su teff22 if first_treat == 2

// estimate: (g, t) = (2, 5) & true effect = relative time
qui reg Y51 X lam if first_treat == 0, nocons
qui predict Y51hat, xb
qui gen teff25 = Y51 - Y51hat
su teff25 if first_treat == 2

// ------------------------------------------------------------------- csdid ---
csdid Y X, ivar(pid) time(tid) gvar(first_treat) notyet
