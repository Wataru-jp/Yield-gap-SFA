*** 
* Author: Wataru Kodama
* Purpose: efficiency analysis for SEA rice bowls FINAL
* Date: Nov 5, 2021
* Update: Nov 5, 2021
* Note:
***


*** Myanmar 
cd /Users/kodam1/Documents/document/研究/YieldGap/Data_final/dta
use Myanmar_SFA2.dta, clear
gen sowing = abs(planting_date - 3)
replace sowing = abs(planting_date - 37) if season == 1
replace ln_N=asinh(N)
replace ln_P=asinh(P)
replace ln_K=asinh(K)
gen ln_exp=asinh(farm_experience)
replace variety2=1 if variety1==1 
gen income = rincome_season/tincome*100
* Dry season 
eststo model1d: sfcross ln_yield ln_N variety2 variety3 sowing herbicide if season==0, ///
emean( fertsplits ln_exp income ) distribution(tnormal) nolog
estadd scalar sigma =  e(sigma_u)+e(sigma_v)
estadd scalar lambda2 =  e(sigma_u)/(e(sigma_u)+e(sigma_v))
predict exp_u_d, jlms
su exp_u_d if season==0
estadd scalar tscore =  r(mean)
* Wet season 
eststo model1w: sfcross ln_yield ln_N variety2 variety3 sowing if season==1 & ln_yield>6.5, ///
emean( fertsplits ln_exp income) distribution(tnormal) nolog
estadd scalar sigma =  e(sigma_u)+e(sigma_v)
estadd scalar lambda2 =  e(sigma_u)/(e(sigma_u)+e(sigma_v))
predict exp_u_w, jlms
su exp_u_w if season==1
estadd scalar tscore =  r(mean)
* save temp 
gen exp_u = exp_u_d
replace exp_u = exp_u_w if season==1
replace country = "Myanmar"
tempfile Myanmar
save `Myanmar', replace


*** Vietnam 
cd /Users/kodam1/Documents/document/研究/YieldGap/Data_final/dta
use Vietnam_SFA2.dta, clear
gen sowing = abs(planting_date - 48)
replace sowing = abs(planting_date - 18) if season == 1
replace ln_N=asinh(N)
replace ln_P=asinh(P)
replace ln_K=asinh(K)
gen ln_exp=asinh(farm_experience)
gen income = rincome_season/tincome*100
* Dry season 
eststo model2d: sfcross ln_yield ln_N ln_P ln_K variety1 sowing herbicide if season==0, ///
emean( fertsplits ln_exp income) distribution(tnormal) nolog
estadd scalar sigma =  e(sigma_u)+e(sigma_v)
estadd scalar lambda2 =  e(sigma_u)/(e(sigma_u)+e(sigma_v))
predict exp_u_d, jlms
su exp_u_d if season==0
estadd scalar tscore =  r(mean)
* Wet season 
eststo model2w: sfcross ln_yield ln_N ln_P ln_K variety1 sowing if season==1, ///
emean( fertsplits ln_exp income) distribution(tnormal) nolog
estadd scalar sigma = e(sigma_u)+e(sigma_v)
estadd scalar lambda2 = e(sigma_u)/(e(sigma_u)+e(sigma_v))
predict exp_u_w, jlms
su exp_u_d if season==1
estadd scalar tscore =  r(mean)
* save temp 
gen exp_u = exp_u_d
replace exp_u = exp_u_w if season==1
replace country = "Vietnam"
tempfile Vietnam
save `Vietnam', replace


*** Thailand
cd /Users/kodam1/Documents/document/研究/YieldGap/Data_final/dta
use Thailand_SFA2.dta, clear
gen sowing = abs(planting_date - 51)
replace sowing = abs(planting_date - 22) if season == 1
//drop if yield<200
replace ln_N=asinh(N)
replace ln_P=asinh(P)
replace ln_K=asinh(K)
gen ln_exp=asinh(farm_experience)
gen income = rincome_season/tincome*100
* Dry season 
eststo model3d: sfcross ln_yield ln_N
* Wet season 
eststo model3w: sfcross ln_yield ln_N ln_P ln_K variety2 sowing herbicide if season==1, ///
emean( fertsplits ln_exp income) distribution(tnormal) nolog
estadd scalar sigma =  e(sigma_u)+e(sigma_v)
estadd scalar lambda2 =  e(sigma_u)/(e(sigma_u)+e(sigma_v))
predict exp_u_w, jlms
predict exp_u_w2, u
su exp_u_w if season==1
estadd scalar tscore =  r(mean)
* save temp 
gen exp_u = .
replace exp_u = exp_u_w if season==1
replace country = "Thailand"
tempfile Thailand
save `Thailand', replace


*** Indonesia
cd /Users/kodam1/Documents/document/研究/YieldGap/Data_final/dta
use Indonesia_SFA2.dta, clear
gen sowing = abs(planting_date - 23)
replace sowing = abs(planting_date - 51) if season == 1
replace ln_N=asinh(N)
replace ln_P=asinh(P)
replace ln_K=asinh(K)
gen ln_exp=asinh(farm_experience)
gen income = rincome_season/tincome*100
* Dry season 
eststo model4d: sfcross ln_yield ln_N ln_P variety1 sowing if season==0 & ln_yield>7, ///
emean( fertsplits ln_exp income) distribution(tnormal) nolog
estadd scalar sigma =  e(sigma_u)+e(sigma_v)
estadd scalar lambda2 =  e(sigma_u)/(e(sigma_u)+e(sigma_v))
predict exp_u_d, jlms
su exp_u_d if season==0
estadd scalar tscore =  r(mean)
* Wet season 
eststo model4w: sfcross ln_yield ln_N ln_P variety1 sowing if season==1, ///
emean( fertsplits ln_exp income) distribution(tnormal) nolog 
estadd scalar sigma =  e(sigma_u)+e(sigma_v)
estadd scalar lambda2 =  e(sigma_u)/(e(sigma_u)+e(sigma_v))
predict exp_u_w, jlms
su exp_u_w if season==1
estadd scalar tscore =  r(mean)
* save temp 
gen exp_u = exp_u_d
replace exp_u = exp_u_w if season==1
replace country = "Indonesia"
tempfile Indonesia
save `Indonesia', replace



** output: combine data
use `Myanmar', clear
drop seedling_age
append using `Vietnam'
drop seedling_age
append using `Thailand'
drop seedling_age
append using `Indonesia'
drop seedling_age
drop if yield == .
label define season0 0 "DRY SEASON" 1 "WET SEASON" 
label values season season0
* variable
replace yield = yield/1000
gen y = yield/exp_u
bysort country: su exp_u
gen ygap = 1-exp_u

twoway kdensity ygap if country=="Myanmar", color(gs4) graphregion(color(white)) ///
xlabel(0(0.2)1) ylabel(0(2)10) yscale(range(0 4))	///
|| kdensity ygap if country=="Vietnam", color(midblue) ///
|| kdensity ygap if country=="Thailand"& ygap<0.8, color(dkorange) ///
|| kdensity ygap if country=="Indonesia", color(midgreen) ///
by(season, note("") graphregion(fcolor(white))) ///
ytitle("Density") xtitle("Efficiency yield gap")  ///
legend(order(1 "Myanmar" 2 "Vietnam" 3 "Thailand" 4 "Indonesia")) 

stop 
* Fig 5
twoway kdensity ygap if country=="Myanmar", color(gs4) graphregion(color(white)) ///
xlabel(0(0.2)1) ylabel(0 " ") yscale(range(0 5))	///
|| kdensity ygap if country=="Vietnam", color(midblue) ///
|| kdensity ygap if country=="Thailand", color(dkorange) ///
|| kdensity ygap if country=="Indonesia", color(midgreen) ///
by(season, note("") graphregion(fcolor(white))) ///
ytitle("Density") xtitle("Efficiency yield gap")  ///
legend(order(1 "Myanmar" 2 "Vietnam" 3 "Thailand" 4 "Indonesia")) 
graph export "/Users/kodam1/Desktop/YieldGap/SFA/latex/fig1.jpg", as(jpg) name("Graph") quality(100) replace
 
* Fig 6
twoway scatter y yield if country=="Myanmar", msize(tiny) mcolor(gs4) xlabel(0(2)11) ylabel(0(2)13) graphregion(color(white)) ///
|| scatter y yield if country=="Vietnam" & y<12.5, msize(tiny) mcolor(midblue) msymbol(0) ///
|| scatter y yield if country=="Thailand", msize(tiny) mcolor(dkorange) msymbol(0) ///
|| scatter y yield if country=="Indonesia" & y<12.5, msize(tiny) mcolor(midgreen) msymbol(0) ///
|| function y=0.99*x, range(0 11) color(gs5) lpattern(solid) /// 
by(season, note("") graphregion(fcolor(white))) ///
ytitle("Technical efficient yield (t/ha)") xtitle("Actual yield (t/ha)") ///
legend(order(1 "Myanmar" 2 "Vietnam" 3 "Thailand" 4 "Indonesia")) 
graph export "/Users/kodam1/Desktop/YieldGap/SFA/latex/fig2.jpg", as(jpg) name("Graph") quality(100) replace

twoway scatter y yield if country=="Myanmar" & season==0, msize(tiny) mcolor(gs4) xlabel(0(2)11) ylabel(0(2)15) graphregion(color(white)) xsize(5) ysize(4.5) ///
|| scatter y yield if country=="Vietnam" & season==0 & y<14, msize(tiny) mcolor(midblue) msymbol(0) ///
|| scatter y yield if country=="Thailand" & season==0, msize(tiny) mcolor(dkorange) msymbol(0) ///
|| scatter y yield if country=="Indonesia" & season==0, msize(tiny) mcolor(midgreen) msymbol(0) ///
|| function y=0.99*x, range(0 11) color(gs5) lpattern(solid) /// 
ytitle("Technical efficient yield (t/ha)") xtitle("Actual yield (t/ha)") ///
legend(order(1 "Myanmar" 2 "Vietnam" 3 "Thailand" 4 "Indonesia")) 
graph export "/Users/kodam1/Desktop/YieldGap/SFA/latex/fig2a.jpg", as(jpg) name("Graph") quality(100) replace

twoway scatter y yield if country=="Myanmar" & season==1, msize(tiny) mcolor(gs4) xlabel(0(2)11) ylabel(0(2)13) graphregion(color(white)) xsize(5) ysize(4.5) ///
|| scatter y yield if country=="Vietnam" & season==1 & y<13, msize(tiny) mcolor(midblue) msymbol(0) ///
|| scatter y yield if country=="Thailand" & season==1, msize(tiny) mcolor(dkorange) msymbol(0) ///
|| scatter y yield if country=="Indonesia" & season==1, msize(tiny) mcolor(midgreen) msymbol(0) ///
|| function y=0.99*x, range(0 11) color(gs5) lpattern(solid) /// 
ytitle("Technical efficient yield (t/ha)") xtitle("Actual yield (t/ha)") ///
legend(order(1 "Myanmar" 2 "Vietnam" 3 "Thailand" 4 "Indonesia")) 
graph export "/Users/kodam1/Desktop/YieldGap/SFA/latex/fig2b.jpg", as(jpg) name("Graph") quality(100) replace




* Table
replace actual_yield = actual_yield/1000
label var actual_yield "Yield (t ha$^{-1}$)"
label var farm_size "Farm size (ha)"
label var N "N applied (kg N ha$^{-1}$)"
label var P "P applied (kg P ha$^{-1}$)"
label var K "K applied (kg P ha$^{-1}$)"
label var rice_ecosystem "Irrigated lowland"
label var sowing "Sowing week"
label var fertsplits "Fertilizer splits (\#)"
label var herbicide "Herbicide use"
label var ce "Transplanting"
label var threshing "Threshing machine use"
label var harvest "Mechanized harvest"
label var education "Education"
label var farm_experience "Farm experience (year)"
label var land_ownership "Land ownership"
label var rincome_season "Rice income (USD ha$^{-1}$)"
label var ln_exp "Farm experience, log"
replace rice_ecosystem=1 if country!="Myanmar"
label var variety1 "Short duration"
label var variety2 "Medium-short duration"
label var variety3 "Medium-long duration"
label var variety4 "Long duration"
label var farm_experience "Farm experience"
replace variety3=0 if country=="Indonesia"
replace variety4=0 if country=="Indonesia"
		
** Table 1 
cd "/Users/kodam1/Documents/document/研究/YieldGap/SFA/latex"
gen grp = "Myanmar DS"
replace grp = "Vietnam DS" if country == "Vietnam" & season==0
replace grp = "Thailand DS" if country == "Thailand" & season==0
replace grp = "Indonesia DS" if country == "Indonesia" & season==0
replace grp = "Myanmar WS" if country == "Myanmar" & season==1
replace grp = "Vietnam WS" if country == "Vietnam" & season==1
replace grp = "Thailand WS" if country == "Thailand" & season==1
replace grp = "Indonesia WS" if country == "Indonesia" & season==1

eststo dstat: estpost tabstat actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4 herbicide fertsplits threshing education farm_experience, ///
by(grp) statistics(mean sd) columns(statistics) listwise nototal
esttab dstat using table1_1.tex, replace label ///
	unstack nostar nonote nomtitle nonumber ///
    cells(mean(fmt(2)) sd(par)) nogaps compress par ///
    collabels(none)

eststo M0: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income if country == "Myanmar" & season == 0
eststo M1: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income  if country == "Myanmar" & season == 1
eststo V0: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income if country == "Vietnam" & season == 0
eststo V1: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income if country == "Vietnam" & season == 1
eststo T0: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income if country == "Thailand" & season == 0
eststo T1: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income if country == "Thailand" & season == 1
eststo I0: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income if country == "Indonesia" & season == 0
eststo I1: qui estpost su actual_yield farm_size N P K herbicide sowing variety1 variety2 variety3 variety4  fertsplits threshing harvest education farm_experience land_ownership income if country == "Indonesia" & season == 1

cd /Users/kodam1/Documents/document/研究/YieldGap/SFA/latex
esttab M0 M1 V0 V1 T0 T1 I0 I1 using table1.tex, replace ///
label cells(mean(pattern(1 1 1 1 1 1 1 1) fmt(2)) sd(par)) nogap nonotes ///
nonumber nomtitles ///
mgroups("Myanmar" "Vietnam" "Thailand" "Indonesia", ///
		pattern(1 0 1 0 1 0 1 0) ///
		span prefix(\multicolumn{@span}{c}{) suffix(}) ///
		erepeat(\cmidrule(lr){@span})) 

cd /Users/kodam1/Documents/document/研究/YieldGap/MSword
esttab M0 M1 V0 V1 T0 T1 I0 I1 using table1.rtf, replace ///
label cells(mean(pattern(1 1 1 1 1 1 1 1) fmt(2)) sd(par)) nogap nonotes ///
nonumber nomtitles ///
mgroups("Myanmar" "Vietnam" "Thailand" "Indonesia", ///
		pattern(1 0 1 0 1 0 1 0)) 




** Table 3
cd /Users/kodam1/Documents/document/研究/YieldGap/SFA/latex
label var Area "Farm size log"
label var ln_exp "Farm experience log"
label var herbicide "Herbicide use"
label var income "Ratio of rice income"
label var variety1 "Short duration"
label var variety2 "Medium-short duration"
label var variety3 "Medium-long duration"
label var variety4 "Long duration"
* Table 
esttab model1d model1w model2d model2w model3w model4d model4w using table3.tex, replace ///
	se label nogap nonotes b(%4.3f) ///
	star(* 0.1 ** 0.05 *** 0.01)  ///
	order( ln_N ln_P ln_K herbicide sowing variety1 variety2 variety3 fertsplits herbicide ln_exp) ///
	s(chi2 tscore sigma lambda2 N, fmt(%9.2f %9.3f %9.3f %9.3g %9.0g) labels("Wald Chi-squared" "TE score" "$\sigma^2=\sigma^2_u+\sigma^2_v$" "$\lambda=\sigma^2_u/\sigma^2$" "Observations")) ///
mgroups("Myanmar" "Vietnam" "Thailand" "Indonesia", ///
		pattern(1 0 1 0 1 1 0) ///
		span prefix(\multicolumn{@span}{c}{) suffix(}))

cd /Users/kodam1/Documents/document/研究/YieldGap/MSword
esttab model1d model1w model2d model2w model3w model4d model4w using table3.rtf, replace ///
	se label nogap nonotes b(%4.3f) ///
	star(* 0.1 ** 0.05 *** 0.01)  ///
	order( ln_N ln_P ln_K herbicide sowing variety1 variety2 variety3 fertsplits herbicide ln_exp) ///
	s(chi2 tscore sigma lambda2 N, fmt(%9.2f %9.3f %9.3f %9.3g %9.0g) labels("Wald Chi-squared" "TE score" "$\sigma^2=\sigma^2_u+\sigma^2_v$" "$\lambda=\sigma^2_u/\sigma^2$" "Observations")) ///
mgroups("Myanmar" "Vietnam" "Thailand" "Indonesia", ///
		pattern(1 0 1 0 1 1 0))

		
keep(ln_N ln_P ln_K variety1 variety2 variety3 planting_RANK fertsplits herbicide threshing education ln_exp _cons) ///
