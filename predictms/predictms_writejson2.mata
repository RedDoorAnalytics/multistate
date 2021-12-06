version 15.1

mata:
function predictms_writejson2() {
	filename = st_local("jsonfile")
	
	Nstates = strtoreal(st_local("Nstates"))
	Ntrans = strtoreal(st_local("Ntrans"))
	Nats = strtoreal(st_local("Nats"))
	from = strtoreal(tokens(st_local("from")))

	refat=1
	//mint=strtoreal(st_local("mint"))
	Natscomp=Nats-1
	hasprob=st_local("probability")!=""
	haslos=st_local("los")!=""
	hashaz=st_local("hazard")!=""
    hasvisit=st_local("visit")!=""
	hasdiff=st_local("difference")!=""
	hasratio=st_local("ratio")!=""
	hasuser=st_local("userfunction")!=""
	hasci=st_local("ci")!=""
	touse = st_local("touse")
	
	hastimevar=st_local("timevar")!=""
    hasmint=st_local("mint")!=""
	hasmaxt=st_local("maxt")!=""

	tmat = st_matrix(st_local("transmatrix"))
	ncoltmat=cols(tmat)
	nrowtmat=rows(tmat)
		 
	if (hasmint==1 | hasmaxt==1 ) {
		timevar = st_data(.,"_time",touse)'
		timevar = invtokens(strofreal(timevar,"%9.5f"),",")
	}
	
	if (hastimevar==1) {
		timevar = st_data(.,st_local("timevar"),touse)'
		timevar = invtokens(strofreal(timevar,"%9.5f"),",")
	}
	
	atlist = J(1,Nats,"")
	
	for(i=1;i<=Nats;i++) {
		atlist[1,i] = `"""' + st_local("at"+strofreal(i)) + `"""'
	}
	atlist = "[" + invtokens(atlist,",") + "]"

	
// open file
	fh = fopen(filename, "w")
	
	fput(fh,"{")
	
	fput(fh,`""timevar": ["' + timevar + "],")
	
	fput(fh,`""Nats": "' + strofreal(Nats) +",")
	
	fput(fh,`""Ntransitions": "' + strofreal(Ntrans) +",")
	
    tmptmat  = `""tmat":["'  
        for (j=1; j<=nrowtmat; j++) {  
        	tmptmat = tmptmat + "[" + invtokens(strofreal(tmat[j,]),",") + "]"
            if(j!=nrowtmat) tmptmat = tmptmat +", "
        }
        tmptmat  = tmptmat + "],"
		tmptmat = subinstr(tmptmat,".", `"""'+"NA"+`"""')
        fput(fh,tmptmat)

	
	/************Main estimates*********************/
	
for (g=1; g<=cols(from);g++) {	

end_states=select( 1..cols(tmat[from[1,g],]), tmat[from[1,g],] :!= .)
	
	
if (hashaz==1) {
	
		for(j=1;j<=length(end_states);j++) {
		   Hreturn = J(1,Nats,"")
		   
		
		   
		         for(ati=1;ati<=Nats;ati++) {
			         tmpH = st_data(.,"_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j]),touse)'
			         Hreturn[1,ati] = "[" + invtokens(strofreal(tmpH,"%9.5f"),",") + "]"
		         }
		   
			    Hname = "Haz_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])
			    Hreturn = invtokens(Hreturn,",")
			    Hreturn = subinstr(Hreturn,".,","null,")
			    Hreturn = subinstr(Hreturn,".]","null]")
			    fput(fh,`"""' + Hname + `"": ["' + Hreturn + "],")
	      
		}

	if (hasdiff==1) {
		for(j=1;j<=length(end_states);j++) {
		

		
				Hdreturn = J(1,Natscomp,"")
				for(ati=2;ati<=Nats;ati++) {

						tmpHd = st_data(.,"_diff_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j]),touse)'
						Hdreturn[1,ati-1] = "[" + invtokens(strofreal(tmpHd,"%9.5f"),",") + "]"
				}
				
				Hdname = "Haz_diff_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])
				Hdreturn = invtokens(Hdreturn,",")
				Hdreturn = subinstr(Hdreturn,".,","null,")
				Hdreturn = subinstr(Hdreturn,".]","null]")
				fput(fh,`"""' + Hdname + `"": ["' + Hdreturn + "],")
		 
	   }
	}
	
	if (hasratio==1) {
	
		 for(j=1;j<=length(end_states);j++) {
		 
				  
						Hrreturn = J(1,Natscomp,"")
						
						for(ati=2;ati<=Nats;ati++) {
						tmpHr = st_data(.,"_ratio_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j]),touse)'
						Hrreturn[1,ati-1] = "[" + invtokens(strofreal(tmpHr,"%9.5f"),",") + "]"
				        }
						
						Hrname = "Haz_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])
						Hrreturn = invtokens(Hrreturn,",")
						Hrreturn = subinstr(Hrreturn,".,","null,")
						Hrreturn = subinstr(Hrreturn,".]","null]")
						fput(fh,`"""' + Hrname + `"": ["' + Hrreturn + "],")
			 
		}						
  }

}

if (hasprob==1) {
	

	
			for(j=1;j<=Nstates;j++) {
				Preturn = J(1,Nats,"")
				for(ati=1;ati<=Nats;ati++) {
					tmpP = st_data(.,"_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
					Preturn[1,ati] = "[" + invtokens(strofreal(tmpP,"%9.5f"),",") + "]"
				}
			Pname = "P_"+strofreal(from[1,g])+"_to_"+strofreal(j)
			Preturn = invtokens(Preturn,",")
			Preturn = subinstr(Preturn,".,","null,")
			Preturn = subinstr(Preturn,".]","null]")
			fput(fh,`"""' + Pname + `"": ["' + Preturn + "],")
			}
	
	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		  Pdreturn = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {

			  tmpPd = st_data(.,"_diff_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			  Pdreturn[1,ati-1] = "[" + invtokens(strofreal(tmpPd,"%9.5f"),",") + "]"
		  }
		  Pdname = "P_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		  Pdreturn = invtokens(Pdreturn,",")
		  Pdreturn = subinstr(Pdreturn,".,","null,")
		  Pdreturn = subinstr(Pdreturn,".]","null]")
		  fput(fh,`"""' + Pdname + `"": ["' + Pdreturn + "],")
	    }
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {
				  
		  Prreturn = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {
			  tmpPr = st_data(.,"_ratio_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			  Prreturn[1,ati-1] = "[" + invtokens(strofreal(tmpPr,"%9.5f"),",") + "]"
		  }
		  Prname = "P_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		  Prreturn = invtokens(Prreturn,",")
		  Prreturn = subinstr(Prreturn,".,","null,")
		  Prreturn = subinstr(Prreturn,".]","null]")
		  fput(fh,`"""' + Prname + `"": ["' + Prreturn + "],")
	    }
	}	

}
	
if (hasvisit==1) {
	
		for(j=1;j<=Nstates;j++) {
		   Vreturn = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpV = st_data(.,"_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			   Vreturn[1,ati] = "[" + invtokens(strofreal(tmpV,"%9.5f"),",") + "]"
		   }
		Vname = "Visit_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		Vreturn = invtokens(Vreturn,",")
		Vreturn = subinstr(Vreturn,".,","null,")
		Vreturn = subinstr(Vreturn,".]","null]")
		fput(fh,`"""' + Vname + `"": ["' + Vreturn + "],")
	    }

	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		  Vdreturn = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {

			  tmpVd = st_data(.,"_diff_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			  Vdreturn[1,ati-1] = "[" + invtokens(strofreal(tmpVd,"%9.5f"),",") + "]"
		  }
		  Vdname = "Visit_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		  Vdreturn = invtokens(Vdreturn,",")
		  Vdreturn = subinstr(Vdreturn,".,","null,")
		  Vdreturn = subinstr(Vdreturn,".]","null]")
		  fput(fh,`"""' + Vdname + `"": ["' + Vdreturn + "],")
	    }
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {
				  
		  Vrreturn = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {
			  tmpVr = st_data(.,"_ratio_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			  Vrreturn[1,ati-1] = "[" + invtokens(strofreal(tmpVr,"%9.5f"),",") + "]"
		  }
		  Vrname = "Visit_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		  Vrreturn = invtokens(Vrreturn,",")
		  Vrreturn = subinstr(Vrreturn,".,","null,")
		  Vrreturn = subinstr(Vrreturn,".]","null]")
		  fput(fh,`"""' + Vrname + `"": ["' + Vrreturn + "],")
	    }
	}		
			
		
}

	
if (haslos==1) {	
	
		for(j=1;j<=Nstates;j++) {
		   Lreturn = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpL = st_data(.,"_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			   Lreturn[1,ati] = "[" + invtokens(strofreal(tmpL,"%9.5f"),",") + "]"
		   }
	    Lname = "Los_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		Lreturn = invtokens(Lreturn,",")
		Lreturn = subinstr(Lreturn,".,","null,")
		Lreturn = subinstr(Lreturn,".]","null]")
		fput(fh,`"""' + Lname + `"": ["' + Lreturn + "],")
	    }
	
	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		
		   Ldreturn = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpLd = st_data(.,"_diff_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			   Ldreturn[1,ati-1] = "[" + invtokens(strofreal(tmpLd,"%9.5f"),",") + "]"
		   }
	    Ldname = "Los_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		Ldreturn = invtokens(Ldreturn,",")
	    Ldreturn = subinstr(Ldreturn,".,","null,")
		Ldreturn = subinstr(Ldreturn,".]","null]")
		fput(fh,`"""' + Ldname + `"": ["' + Ldreturn + "],")
	    }	
	
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {

		   Lrreturn = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpLr = st_data(.,"_ratio_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j),touse)'
			   Lrreturn[1,ati-1] = "[" + invtokens(strofreal(tmpLr,"%9.5f"),",") + "]"
		   }
	    Lrname = "Los_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)
		Lrreturn = invtokens(Lrreturn,",")
		Lrreturn = subinstr(Lrreturn,".,","null,")
	    Lrreturn = subinstr(Lrreturn,".]","null]")
		Lrreturn = subinstr(Lrreturn,".]","null]")
		fput(fh,`"""' + Lrname + `"": ["' + Lrreturn + "],")
	    }	
	
	}
}

if (hasuser==1) {	
	
		
		   Ureturn = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpU = st_data(.,"_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1",touse)'
			   Ureturn[1,ati] = "[" + invtokens(strofreal(tmpU,"%9.5f"),",") + "]"
		   }
	    Uname = "User_"+strofreal(from[1,g])
		Ureturn = invtokens(Ureturn,",")
		Ureturn = subinstr(Ureturn,".,","null,")
		Ureturn = subinstr(Ureturn,".]","null]")
		fput(fh,`"""' + Uname + `"": ["' + Ureturn + "],")
	   
	
	if (hasdiff==1) {
	
		   Udreturn = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpUd = st_data(.,"_diff_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1",touse)'
			   Udreturn[1,ati-1] = "[" + invtokens(strofreal(tmpUd,"%9.5f"),",") + "]"
		   }
	    Udname = "User_diff_"+strofreal(from[1,g])
		Udreturn = invtokens(Udreturn,",")
	    Udreturn = subinstr(Udreturn,".,","null,")
		Udreturn = subinstr(Udreturn,".]","null]")
		fput(fh,`"""' + Udname + `"": ["' + Udreturn + "],")
	    	
	
	}
	
	if (hasratio==1) {

		   Urreturn = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpUr = st_data(.,"_ratio_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1",touse)'
			   Urreturn[1,ati-1] = "[" + invtokens(strofreal(tmpUr,"%9.5f"),",") + "]"
		   }
	    Urname = "User_ratio_"+strofreal(from[1,g])
		Urreturn = invtokens(Urreturn,",")
		Urreturn = subinstr(Urreturn,".,","null,")
		Urreturn = subinstr(Urreturn,".]","null]")
		fput(fh,`"""' + Urname + `"": ["' + Urreturn + "],")
	    	
	
	}
}

/*********************************************************************************************/

if (hasci==1) {

/************Upper limit ci*********************/


if (hashaz==1) {
	
		for(j=1;j<=length(end_states);j++) {
		   Hreturn_uci = J(1,Nats,"")
		   

		   
		         for(ati=1;ati<=Nats;ati++) {
			         tmpH_uci = st_data(.,"_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j])+"_uci",touse)'
			         Hreturn_uci[1,ati] = "[" + invtokens(strofreal(tmpH_uci,"%9.5f"),",") + "]"
		         }
		   
			    Hname_uci = "Haz_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])+"_uci"
			    Hreturn_uci = invtokens(Hreturn_uci,",")
			    Hreturn_uci = subinstr(Hreturn_uci,".,","null,")
			    Hreturn_uci = subinstr(Hreturn_uci,".]","null]")
			    fput(fh,`"""' + Hname_uci + `"": ["' + Hreturn_uci + "],")
	      
		}

	if (hasdiff==1) {
		for(j=1;j<=length(end_states);j++) {
		

		
				Hdreturn_uci = J(1,Natscomp,"")
				for(ati=2;ati<=Nats;ati++) {

						tmpHd_uci = st_data(.,"_diff_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j])+"_uci",touse)'
						Hdreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpHd_uci,"%9.5f"),",") + "]"
				}
				
				Hdname_uci = "Haz_diff_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])+"_uci"
				Hdreturn_uci = invtokens(Hdreturn_uci,",")
				Hdreturn_uci = subinstr(Hdreturn_uci,".,","null,")
				Hdreturn_uci = subinstr(Hdreturn_uci,".]","null]")
				fput(fh,`"""' + Hdname_uci + `"": ["' + Hdreturn_uci + "],")
		 
	   }
	}
	
	if (hasratio==1) {
	
		 for(j=1;j<=length(end_states);j++) {
		 

				  
						Hrreturn_uci = J(1,Natscomp,"")
						
						for(ati=2;ati<=Nats;ati++) {
						tmpHr_uci = st_data(.,"_ratio_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j])+"_uci",touse)'
						Hrreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpHr_uci,"%9.5f"),",") + "]"
				        }
						
						Hrname_uci = "Haz_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])+"_uci"
						Hrreturn_uci = invtokens(Hrreturn_uci,",")
						Hrreturn_uci = subinstr(Hrreturn_uci,".,","null,")
						Hrreturn_uci = subinstr(Hrreturn_uci,".]","null]")
						fput(fh,`"""' + Hrname_uci + `"": ["' + Hrreturn_uci + "],")
			 
		}						
  }

}
	
	
if (hasprob==1) {
	for(j=1;j<=Nstates;j++) {
		Preturn_uci = J(1,Nats,"")
		for(ati=1;ati<=Nats;ati++) {
			tmpP_uci = st_data(.,"_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			Preturn_uci[1,ati] = "[" + invtokens(strofreal(tmpP_uci,"%9.5f"),",") + "]"
		}
	    Pname_uci = "P_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		Preturn_uci = invtokens(Preturn_uci,",")
		Preturn_uci = subinstr(Preturn_uci,".,","null,")
		Preturn_uci = subinstr(Preturn_uci,".]","null]")
		fput(fh,`"""' + Pname_uci + `"": ["' + Preturn_uci + "],")
	}

	
	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		  Pdreturn_uci = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {

			  tmpPd_uci = st_data(.,"_diff_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			  Pdreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpPd_uci,"%9.5f"),",") + "]"
		  }
		  Pdname_uci = "P_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		  Pdreturn_uci = invtokens(Pdreturn_uci,",")
		  Pdreturn_uci = subinstr(Pdreturn_uci,".,","null,")
		  Pdreturn_uci = subinstr(Pdreturn_uci,".]","null]")
		  fput(fh,`"""' + Pdname_uci + `"": ["' + Pdreturn_uci + "],")
	    }
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {
				  
		  Prreturn_uci = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {
			  tmpPr_uci = st_data(.,"_ratio_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			  Prreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpPr_uci,"%9.5f"),",") + "]"
		  }
		  Prname_uci = "P_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		  Prreturn_uci = invtokens(Prreturn_uci,",")
		  Prreturn_uci = subinstr(Prreturn_uci,".,","null,")
		  Prreturn_uci = subinstr(Prreturn_uci,".]","null]")
		  fput(fh,`"""' + Prname_uci + `"": ["' + Prreturn_uci + "],")
	    }
	}
	

}

if (hasvisit==1) {
	
		for(j=1;j<=Nstates;j++) {
		   Vreturn_uci = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpV_uci = st_data(.,"_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			   Vreturn_uci[1,ati] = "[" + invtokens(strofreal(tmpV_uci,"%9.5f"),",") + "]"
		   }
		Vname_uci = "Visit_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		Vreturn_uci = invtokens(Vreturn_uci,",")
		Vreturn_uci = subinstr(Vreturn_uci,".,","null,")
		Vreturn_uci = subinstr(Vreturn_uci,".]","null]")
		fput(fh,`"""' + Vname_uci + `"": ["' + Vreturn_uci + "],")
	    }
		
		
	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
			Vdreturn_uci = J(1,Natscomp,"")
			for(ati=2;ati<=Nats;ati++) {
			   tmpVd_uci = st_data(.,"_diff_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			   Vdreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpVd_uci,"%9.5f"),",") + "]"
			}
			Vdname_uci 		= "Visit_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
			Vdreturn_uci 	= invtokens(Vdreturn_uci,",")
			Vdreturn_uci 	= subinstr(Vdreturn_uci,".,","null,")
			Vdreturn_uci 	= subinstr(Vdreturn_uci,".]","null]")
			fput(fh,`"""' + Vdname_uci + `"": ["' + Vdreturn_uci + "],")
	    }	
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {

		   Vrreturn_uci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpVr_uci = st_data(.,"_ratio_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			   Vrreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpVr_uci,"%9.5f"),",") + "]"
		   }
	    Vrname_uci = "Visit_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		Vrreturn_uci = invtokens(Vrreturn_uci,",")
		Vrreturn_uci = subinstr(Vrreturn_uci,".,","null,")
		Vrreturn_uci = subinstr(Vrreturn_uci,".]","null]")
		fput(fh,`"""' + Vrname_uci + `"": ["' + Vrreturn_uci + "],")
	    }	
	
	}	
			
}

	
if (haslos==1) {	
	
		for(j=1;j<=Nstates;j++) {
		   Lreturn_uci = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpL_uci = st_data(.,"_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			   Lreturn_uci[1,ati] = "[" + invtokens(strofreal(tmpL_uci,"%9.5f"),",") + "]"
		   }
	    Lname_uci = "Los_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		Lreturn_uci = invtokens(Lreturn_uci,",")
		Lreturn_uci = subinstr(Lreturn_uci,".,","null,")
		Lreturn_uci = subinstr(Lreturn_uci,".]","null]")
		fput(fh,`"""' + Lname_uci + `"": ["' + Lreturn_uci + "],")
	    }
	
	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		
		   Ldreturn_uci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpLd_uci = st_data(.,"_diff_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			   Ldreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpLd_uci,"%9.5f"),",") + "]"
		   }
	    Ldname_uci = "Los_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		Ldreturn_uci = invtokens(Ldreturn_uci,",")
	    Ldreturn_uci = subinstr(Ldreturn_uci,".,","null,")
		Ldreturn_uci = subinstr(Ldreturn_uci,".]","null]")
		fput(fh,`"""' + Ldname_uci + `"": ["' + Ldreturn_uci + "],")
	    }	
	
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {

		   Lrreturn_uci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpLr_uci = st_data(.,"_ratio_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_uci",touse)'
			   Lrreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpLr_uci,"%9.5f"),",") + "]"
		   }
	    Lrname_uci = "Los_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_uci"
		Lrreturn_uci = invtokens(Lrreturn_uci,",")
		Lrreturn_uci = subinstr(Lrreturn_uci,".,","null,")
		Lrreturn_uci = subinstr(Lrreturn_uci,".]","null]")
		fput(fh,`"""' + Lrname_uci + `"": ["' + Lrreturn_uci + "],")
	    }	
	
	}
}

if (hasuser==1) {	
	
		
		   Ureturn_uci = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpU_uci = st_data(.,"_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1"+"_uci",touse)'
			   Ureturn_uci[1,ati] = "[" + invtokens(strofreal(tmpU_uci,"%9.5f"),",") + "]"
		   }
	    Uname_uci = "User_"+strofreal(from[1,g])+"_uci"
		Ureturn_uci = invtokens(Ureturn_uci,",")
		Ureturn_uci = subinstr(Ureturn_uci,".,","null,")
	    Ureturn_uci = subinstr(Ureturn_uci,".]","null]")
		fput(fh,`"""' + Uname_uci + `"": ["' + Ureturn_uci + "],")
	   
	
	if (hasdiff==1) {
	
		   Udreturn_uci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpUd_uci = st_data(.,"_diff_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1"+"_uci",touse)'
			   Udreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpUd_uci,"%9.5f"),",") + "]"
		   }
	    Udname_uci = "User_diff_"+strofreal(from[1,g])+"_uci"
		Udreturn_uci = invtokens(Udreturn_uci,",")
	    Udreturn_uci = subinstr(Udreturn_uci,".,","null,")
	    Udreturn_uci = subinstr(Udreturn_uci,".]","null]")
		fput(fh,`"""' + Udname_uci + `"": ["' + Udreturn_uci + "],")
	    	
	
	}
	
	if (hasratio==1) {

		   Urreturn_uci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpUr_uci = st_data(.,"_ratio_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1"+"_uci",touse)'
			   Urreturn_uci[1,ati-1] = "[" + invtokens(strofreal(tmpUr_uci,"%9.5f"),",") + "]"
		   }
	    Urname_uci = "User_ratio_"+strofreal(from[1,g])+"_uci"
		Urreturn_uci = invtokens(Urreturn_uci,",")
		Urreturn_uci = subinstr(Urreturn_uci,".,","null,")
	    Urreturn_uci = subinstr(Urreturn_uci,".]","null]")
		fput(fh,`"""' + Urname_uci + `"": ["' + Urreturn_uci + "],")
	    	
	
	}
}


/*********************************************************************************************/

	


/************Lower limit ci*********************/




if (hashaz==1) {
	
		for(j=1;j<=length(end_states);j++) {
		   Hreturn_lci = J(1,Nats,"")
		   

		   
		         for(ati=1;ati<=Nats;ati++) {
			         tmpH_lci = st_data(.,"_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j])+"_lci",touse)'
			         Hreturn_lci[1,ati] = "[" + invtokens(strofreal(tmpH_lci,"%9.5f"),",") + "]"
		         }
		   
			    Hname_lci = "Haz_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])+"_lci"
			    Hreturn_lci = invtokens(Hreturn_lci,",")
			    Hreturn_lci = subinstr(Hreturn_lci,".,","null,")
			    Hreturn_lci = subinstr(Hreturn_lci,".]","null]")
			    fput(fh,`"""' + Hname_lci + `"": ["' + Hreturn_lci + "],")
	      
		}

	if (hasdiff==1) {
		for(j=1;j<=length(end_states);j++) {
		

		
				Hdreturn_lci = J(1,Natscomp,"")
				for(ati=2;ati<=Nats;ati++) {

						tmpHd_lci = st_data(.,"_diff_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j])+"_lci",touse)'
						Hdreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpHd_lci,"%9.5f"),",") + "]"
				}
				
				Hdname_lci = "Haz_diff_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])+"_lci"
				Hdreturn_lci = invtokens(Hdreturn_lci,",")
				Hdreturn_lci = subinstr(Hdreturn_lci,".,","null,")
				Hdreturn_lci = subinstr(Hdreturn_lci,".]","null]")
				fput(fh,`"""' + Hdname_lci + `"": ["' + Hdreturn_lci + "],")
		 
	   }
	}
	
	if (hasratio==1) {
	
		 for(j=1;j<=length(end_states);j++) {
		 

				  
						Hrreturn_lci = J(1,Natscomp,"")
						
						for(ati=2;ati<=Nats;ati++) {
						tmpHr_lci = st_data(.,"_ratio_hazard_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(end_states[1,j])+"_lci",touse)'
						Hrreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpHr_lci,"%9.5f"),",") + "]"
				        }
						
						Hrname_lci = "Haz_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(end_states[1,j])+"_lci"
						Hrreturn_lci = invtokens(Hrreturn_lci,",")
						Hrreturn_lci = subinstr(Hrreturn_lci,".,","null,")
						Hrreturn_lci = subinstr(Hrreturn_lci,".]","null]")
						fput(fh,`"""' + Hrname_lci + `"": ["' + Hrreturn_lci + "],")
			 
		}						
  }

}
	


	
if (hasprob==1) {

	
	for(j=1;j<=Nstates;j++) {
	
		Preturn_lci = J(1,Nats,"")
		for(ati=1;ati<=Nats;ati++) {
			tmpP_lci = st_data(.,"_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			Preturn_lci[1,ati] = "[" + invtokens(strofreal(tmpP_lci,"%9.5f"),",") + "]"
		}
	    Pname_lci = "P_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		Preturn_lci = invtokens(Preturn_lci,",")
		Preturn_lci = subinstr(Preturn_lci,".,","null,")
		Preturn_lci = subinstr(Preturn_lci,".]","null]")
		fput(fh,`"""' + Pname_lci + `"": ["' + Preturn_lci + "],")
	}
	
	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		  Pdreturn_lci = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {

			  tmpPd_lci = st_data(.,"_diff_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			  Pdreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpPd_lci,"%9.5f"),",") + "]"
		  }
		  Pdname_lci = "P_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		  Pdreturn_lci = invtokens(Pdreturn_lci,",")
		  Pdreturn_lci = subinstr(Pdreturn_lci,".,","null,")
		  Pdreturn_lci = subinstr(Pdreturn_lci,".]","null]")
		  fput(fh,`"""' + Pdname_lci + `"": ["' + Pdreturn_lci + "],")
	    }
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {
				  
		  Prreturn_lci = J(1,Natscomp,"")
		  for(ati=2;ati<=Nats;ati++) {
			  tmpPr_lci = st_data(.,"_ratio_prob_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			  Prreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpPr_lci,"%9.5f"),",") + "]"
		  }
		  Prname_lci = "P_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		  Prreturn_lci = invtokens(Prreturn_lci,",")
		  Prreturn_lci = subinstr(Prreturn_lci,".,","null,")
		  Prreturn_lci = subinstr(Prreturn_lci,".]","null]")
		  fput(fh,`"""' + Prname_lci + `"": ["' + Prreturn_lci + "],")
	    }
	}	

}
	
if (hasvisit==1) {
	
		for(j=1;j<=Nstates;j++) {
		   Vreturn_lci = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpV_lci = st_data(.,"_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			   Vreturn_lci[1,ati] = "[" + invtokens(strofreal(tmpV_lci,"%9.5f"),",") + "]"
		   }
		Vname_lci = "Visit_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		Vreturn_lci = invtokens(Vreturn_lci,",")
		Vreturn_lci = subinstr(Vreturn_lci,".,","null,")
		Vreturn_lci = subinstr(Vreturn_lci,".]","null]")
		fput(fh,`"""' + Vname_lci + `"": ["' + Vreturn_lci + "],")
	    }
		
			if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		
		   Vdreturn_lci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpVd_lci = st_data(.,"_diff_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			   Vdreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpVd_lci,"%9.5f"),",") + "]"
		   }
	    Vdname_lci = "Visit_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		Vdreturn_lci = invtokens(Vdreturn_lci,",")
	    Vdreturn_lci = subinstr(Vdreturn_lci,".,","null,")
		Vdreturn_lci = subinstr(Vdreturn_lci,".]","null]")
		fput(fh,`"""' + Vdname_lci + `"": ["' + Vdreturn_lci + "],")
	    }	
	
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {

		   Vrreturn_lci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpVr_lci = st_data(.,"_ratio_visit_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			   Vrreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpVr_lci,"%9.5f"),",") + "]"
		   }
	    Vrname_lci = "Visit_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		Vrreturn_lci = invtokens(Vrreturn_lci,",")
		Vrreturn_lci = subinstr(Vrreturn_lci,".,","null,")
		Vrreturn_lci = subinstr(Vrreturn_lci,".]","null]")
		fput(fh,`"""' + Vrname_lci + `"": ["' + Vrreturn_lci + "],")
	    }	
	
	}	
}

	
if (haslos==1) {	
	
		for(j=1;j<=Nstates;j++) {
		   Lreturn_lci = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpL_lci = st_data(.,"_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			   Lreturn_lci[1,ati] = "[" + invtokens(strofreal(tmpL_lci,"%9.5f"),",") + "]"
		   }
	    Lname_lci = "Los_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		Lreturn_lci = invtokens(Lreturn_lci,",")
		Lreturn_lci = subinstr(Lreturn_lci,".,","null,")
		Lreturn_lci = subinstr(Lreturn_lci,".]","null]")
		fput(fh,`"""' + Lname_lci + `"": ["' + Lreturn_lci + "],")
	    }
	
	if (hasdiff==1) {
		for(j=1;j<=Nstates;j++) {
		
		
		   Ldreturn_lci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpLd_lci = st_data(.,"_diff_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			   Ldreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpLd_lci,"%9.5f"),",") + "]"
		   }
	    Ldname_lci = "Los_diff_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		Ldreturn_lci = invtokens(Ldreturn_lci,",")
	    Ldreturn_lci = subinstr(Ldreturn_lci,".,","null,")
		Ldreturn_lci = subinstr(Ldreturn_lci,".]","null]")
		fput(fh,`"""' + Ldname_lci + `"": ["' + Ldreturn_lci + "],")
	    }	
	
	}
	
	if (hasratio==1) {
		for(j=1;j<=Nstates;j++) {

		   Lrreturn_lci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpLr_lci = st_data(.,"_ratio_los_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+strofreal(j)+"_lci",touse)'
			   Lrreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpLr_lci,"%9.5f"),",") + "]"
		   }
	    Lrname_lci = "Los_ratio_"+strofreal(from[1,g])+"_to_"+strofreal(j)+"_lci"
		Lrreturn_lci = invtokens(Lrreturn_lci,",")
		Lrreturn_lci = subinstr(Lrreturn_lci,".,","null,")
		Lrreturn_lci = subinstr(Lrreturn_lci,".]","null]")
		fput(fh,`"""' + Lrname_lci + `"": ["' + Lrreturn_lci + "],")
	    }	
	
	}
}

if (hasuser==1) {	
	
		
		   Ureturn_lci = J(1,Nats,"")
		   for(ati=1;ati<=Nats;ati++) {
			   tmpU_lci = st_data(.,"_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1"+"_lci",touse)'
			   Ureturn_lci[1,ati] = "[" + invtokens(strofreal(tmpU_lci,"%9.5f"),",") + "]"
		   }
	    Uname_lci = "User_"+strofreal(from[1,g])+"_to_"+"lci"
		Ureturn_lci = invtokens(Ureturn_lci,",")
		Ureturn_lci = subinstr(Ureturn_lci,".,","null,")
		Ureturn_lci = subinstr(Ureturn_lci,".]","null]")
		fput(fh,`"""' + Uname_lci + `"": ["' + Ureturn_lci + "],")
	   
	
	if (hasdiff==1) {
	
		   Udreturn_lci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpUd_lci = st_data(.,"_diff_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1"+"_lci",touse)'
			   Udreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpUd_lci,"%9.5f"),",") + "]"
		   }
	    Udname_lci = "User_diff_"+strofreal(from[1,g])+"_lci"
		Udreturn_lci = invtokens(Udreturn_lci,",")
	    Udreturn_lci = subinstr(Udreturn_lci,".,","null,")
		Udreturn_lci = subinstr(Udreturn_lci,".]","null]")
		fput(fh,`"""' + Udname_lci + `"": ["' + Udreturn_lci + "],")
	    	
	
	}
	
	if (hasratio==1) {

		   Urreturn_lci = J(1,Natscomp,"")
		   for(ati=2;ati<=Nats;ati++) {
			   tmpUr_lci = st_data(.,"_ratio_user_at"+strofreal(ati)+"_"+strofreal(from[1,g])+"_"+"1"+"_lci",touse)'
			   Urreturn_lci[1,ati-1] = "[" + invtokens(strofreal(tmpUr_lci,"%9.5f"),",") + "]"
		   }
	    Urname_lci = "User_ratio_"+strofreal(from[1,g])+"_lci"
		Urreturn_lci = invtokens(Urreturn_lci,",")
		Urreturn_lci = subinstr(Urreturn_lci,".,","null,")
		Urreturn_lci = subinstr(Urreturn_lci,".]","null]")
		fput(fh,`"""' + Urname_lci + `"": ["' + Urreturn_lci + "],")
	    	
	
	}
}



/*********************************************************************************************/
}  //hasci end

}  // for from loop end

	fput(fh,`""atlist": "' + atlist)
	fput(fh,"}")
	fclose(fh)
	
} //mata command end

end
