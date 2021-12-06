
	capture erase lmultistate.mlib
	mata: mata set matastrict off
	qui {
			do "./predictms/predictms_setup.mata"
			do "./predictms/predictms_functions.mata"
			do "./predictms/predictms.mata"
			do "./predictms/predictms_aj.mata"
			do "./predictms/predictms_analytic.mata"
			do "./predictms/predictms_analytic_survival.mata"
			do "./predictms/predictms_analytic_cr.mata"
			do "./predictms/predictms_analytic_illd.mata"
			do "./predictms/predictms_analytic_extilld.mata"
			do "./predictms/predictms_merlin.mata"
			do "./predictms/predictms_model_predict.mata"
			do "./predictms/predictms_post.mata"
			do "./predictms/predictms_post_dmcis.mata"
			do "./predictms/predictms_post_bootcis.mata"
			do "./predictms/predictms_sim.mata"
			do "./predictms/predictms_sim_cr.mata"
			do "./predictms/predictms_sim_cr_cox.mata"
			do "./predictms/predictms_sim_cr_latent.mata"
			do "./predictms/predictms_sim_root.mata"
// 			do "./predictms/predictms_writejson.mata"
			do "./predictms/predictms_writejson2.mata"
			do "./msaj/msaj.mata"
			
			mata: mata mlib create lmultistate, dir(.)
			mata: mata mlib add    lmultistate *(), dir(.)
			mata: mata d *()  
			mata mata clear
			mata mata mlib index

	}
