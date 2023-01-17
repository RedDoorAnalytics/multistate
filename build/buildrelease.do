// 
// build new version of multistate package
// --> run whole do file

//!!
// --> build in 15
//!!

local drive /Users/Michael/Documents/reddooranalytics/products/multistate
cd `drive'

local includemata = 0

//===========================================================================//

//build SSC release

//current version up is 4_4_0
local newversion 4_4_1
if `includemata' {
	local newversion `newversion'_mata
}
cap mkdir ./release/version_`newversion'
local fdir `drive'/release/version_`newversion'/

//===========================================================================//

//pkg description
copy ./build/multistate_details.txt `fdir', replace
	

copy ./multistate.sthlp `fdir', replace

//predictms, msset, stms
copy ./predictms/predictms.ado `fdir', replace
copy ./predictms/predictms.sthlp `fdir', replace

//msset
copy ./msset/msset.ado `fdir', replace
copy ./msset/msset.sthlp `fdir', replace

//msgraph
copy ./graphms/graphms.ado `fdir', replace
copy ./graphms/graphms.sthlp `fdir', replace

//msaj
copy ./msaj/msaj.ado `fdir', replace
copy ./msaj/msaj.sthlp `fdir', replace

//msboxes
copy ./msboxes/msboxes.ado `fdir', replace
copy ./msboxes/msboxes.sthlp `fdir', replace
copy ./msboxes/msboxes_examples.ado `fdir', replace
copy ./msboxes/freq_total.ado `fdir', replace

//mlib
cap erase `fdir'/lmultistate.mlib
copy ./lmultistate.mlib `fdir', replace

//source mata files
if `includemata' {
	copy ./predictms/predictms.mata `fdir', replace
	copy ./predictms/predictms_aj.mata `fdir', replace
	copy ./predictms/predictms_analytic.mata `fdir', replace
	copy ./predictms/predictms_analytic_cr.mata `fdir', replace
	copy ./predictms/predictms_analytic_extilld.mata `fdir', replace
	copy ./predictms/predictms_analytic_illd.mata `fdir', replace
	copy ./predictms/predictms_analytic_survival.mata `fdir', replace
	copy ./predictms/predictms_functions.mata `fdir', replace
	copy ./predictms/predictms_merlin.mata `fdir', replace
	copy ./predictms/predictms_model_predict.mata `fdir', replace
	copy ./predictms/predictms_pred.mata `fdir', replace
	copy ./predictms/predictms_setup.mata `fdir', replace
	copy ./predictms/predictms_sim.mata `fdir', replace
	copy ./predictms/predictms_sim_cr.mata `fdir', replace
	copy ./predictms/predictms_sim_cr_cox.mata `fdir', replace
	copy ./predictms/predictms_sim_cr_latent.mata `fdir', replace
	copy ./predictms/predictms_sim_root.mata `fdir', replace
	copy ./predictms/predictms_writejson2.mata `fdir', replace
	copy ./msaj/msaj.mata `fdir', replace
}

//===========================================================================//
