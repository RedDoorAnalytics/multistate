version 14.2
mata:

// Main program
void function AJ()
{

	// Import data into mata	
	touse 				= st_local("touse")
						st_view(t,.,"_t",touse)  
						st_view(d,.,"_d",touse)
						st_view(trans,.,"_trans",touse) 
						st_view(Nrisk,.,st_local("Nrisk"), touse)
						st_view(Nevents,.,st_local("Nevents"), touse)
	
	// Read in by data
	if(st_local("by") != "") {
						st_view(by,.,st_local("by"),touse)
		bylevels 		= strtoreal(tokens(st_local("bylevels")))
	}
	else {
		bylevels 		= 1
		by 				= J(rows(t),1,1)
	}

	// Get unique rows for time, trans, haz, from, to, Nrisk and Nevents (only at event times)
	haz 				= Nevents:/Nrisk
	d_index 			= selectindex(d:==1)
	data_all			= uniqrows((t[d_index,], by[d_index], trans[d_index], haz[d_index]))

	// Get ci information
	ci 					= st_local("ci") != ""
	se					= st_local("se") != ""
	if (ci | se) {
						st_view(from,.,"_from",touse)
						st_view(to,.,"_to",touse)
		dataVar_all		= uniqrows((t[d_index,], by[d_index], from[d_index],to[d_index],Nrisk[d_index],Nevents[d_index]))
	}
	
	// Get transmatrix, states, transitions
	transmat 			= st_matrix(st_local("transmatrix"))
	transmat_index 		= transRowCol(transmat)	
	Ntrans 				= check_transmatrix(transmat)
	Nstates 			= rows(transmat)	
	
	// Get enter and exit time
	entertime 			= strtoreal(st_local("ltruncated"))
	exittime 			= strtoreal(st_local("exit"))
	
	// Get from state & check not absorbing
	fromS 				= strtoreal(st_local("from"))
	if (max(transmat[fromS,]) == .) {
		errprintf("From state is an absorbing state\n")
		exit(198)
	}
		
	// Get information for LOS 
	los 				= st_local("los") != ""	
	
	// Create variables to store results
	if (ci) P_state 	= J(rows(t), Nstates*3,.)
	else P_state 		= J(rows(t), Nstates,.)
	if (se) se_state	= J(rows(t), Nstates,.)
	if (los) LOS_state 	= J(rows(t), Nstates,.)
	
	// Get useful matrices
	Imat 				= I(Nstates)
	if (ci | se) zeros 	= J(Nstates, Nstates, 0)
	
	
	// Main loop for calculations, loop over by
	for(k=1;k<=cols(bylevels);k++) {
	
		b 						= bylevels[k]
	
		// Set initial probabilities
		P 						= I(Nstates)
		if (ci | se) varP		= J(Nstates:^2, Nstates:^2, 0)

		// Select data for by level
		by_d_index 				= selectindex(data_all[,2]:==b)
		data					= data_all[by_d_index,]
		if (ci | se) dataVar	= dataVar_all[by_d_index,]	
		t_d_unique 				= uniqrows(data[,1])
		Ntd		 				= rows(t_d_unique)

		if (los) {
			by_index			= selectindex(by:==b)			
			t_unique 			= uniqrows(t[by_index])
			Nt 					= rows(t_unique)
			P_unique 			= J(Nt, Nstates,.)
			P_t_unique 			= J(Nt, 1, .)
			v 					= 1
			past_t				= -1
			LOS 				= J(1, Nstates, 0)
			t0		 			= entertime
			P1 					= J(Nt, Nstates, .) 
			l 					= 1
		}
		
		// Work out the probabilities
		for(u=1; u<=Ntd; u++) {
		
			// Get current t and check before exit time and after enter time
			current_t = t_d_unique[u]
			if (current_t < entertime | current_t > exittime) continue
			
			// If current t = enter time then this will be first current t and P=identity already

			// If current t > enter time then work out the probability
			if (current_t > entertime) {			
				
				// Build hazard matrix for current time
				H = J(Nstates, Nstates,0)
				
				// Fill in the off-diagonals
				for(j=1;j<=Ntrans;j++) {

					if(max((data[,1]:==current_t) :& (data[,3]:==j))==1) {
						H[transmat_index[j,1],transmat_index[j,2]] = select(data[,4],((data[,1]:==current_t) :& (data[,3]:==j)) :==1)			
					}
				}

				// Work out the diagonals
				for(i=1;i<=Nstates;i++) {
					H[i,i]=-quadsum(H[i,])
				}
				
				
				// Work out variance
				if (ci | se) {
					
					// Work out Y_k and N_k., N_kn and N_kl before main loop for efficiency
					RiskEvents		= select(dataVar,dataVar[,1]:==current_t)
					atrisk_k 		= J(1, Nstates, 0)
					events_k		= J(1, Nstates, 0)
					events_kn 		= J(Nstates, Nstates, 0)
					events_kl 		= J(Nstates, Nstates, 0)
		
					for(vk = 1;vk<=Nstates;vk++) {
					
						any_atrisk_k = select(RiskEvents[,5],RiskEvents[,3]:==vk)
						if (rows(any_atrisk_k) != 0) {
							
							atrisk_k[vk] = any_atrisk_k[1]
							events_k[vk] = sum(select(RiskEvents[,6],RiskEvents[,3]:==vk))
							
							for(vl = 1;vl<=Nstates;vl++) {
								any_events_kl = select(RiskEvents[,6],RiskEvents[,3]:==vk :& RiskEvents[,4]:==vl)
								if(rows(any_events_kl) != 0) events_kl[vk, vl] = any_events_kl
							}
							
							for(vn = 1;vn<=Nstates;vn++) {
								any_events_kn = select(RiskEvents[,6],RiskEvents[,3]:==vk :& RiskEvents[,4]:==vn)
								if(rows(any_events_kn) != 0) events_kn[vk, vn] = any_events_kn
							}
						}
					}
					atrisk_k3 = atrisk_k:^(-3)
					
					// Calculate the covariance matrix of (A_kl, A_mn) 
					// Stored in (l,n)th block with row k and column m
					VarBlock = asarray_create("real",2)
					for(vl = 1;vl<=Nstates;vl++) {
					
						for(vn = 1;vn<=Nstates;vn++) {
						
							tempBlock = zeros
							delta_ln = (vl :== vn)
							
							for(vk = 1;vk<=Nstates;vk++) {
							
								if (atrisk_k[vk] != 0) {
								
									if(vk == vl & vk == vn & events_k[vk] != 0) {
										tempBlock[vk,vk] = (atrisk_k[vk] - events_k[vk])*events_k[vk]*atrisk_k3[vk]                                          
									}
																						
									if(vk == vl & vk != vn & events_kn[vk,vn] != 0) {
										tempBlock[vk,vk] = -(atrisk_k[vk] - events_k[vk])*events_kn[vk,vn]*atrisk_k3[vk]
									}
										
									if(vk != vl & vk != vn & events_kn[vk,vn] != 0) {
										tempBlock[vk,vk] = (delta_ln*atrisk_k[vk] - events_kl[vk,vl])*events_kn[vk,vn]*atrisk_k3[vk]
									}
										
								}
							}
							
							asarray(VarBlock,(vl,vn),tempBlock)
							
						}
					}
					
					// Format covariance to be a huge martix rather than array
					VarAd = J(Nstates:^2,Nstates:^2,0)
						for(vl = 1;vl<=Nstates;vl++) {
							for(vn = 1;vn<=Nstates;vn++) {
								trows = ((vl-1)*Nstates+1)..((vl-1)*Nstates+Nstates)
								tcols = ((vn-1)*Nstates+1)..((vn-1)*Nstates+Nstates)
								VarAd[trows,tcols] = asarray(VarBlock,(vl,vn))
							}
						}

					// Make sure it is symmetrical
					for(i = 1; i<= Nstates:^2; i++) {
						for (j = 1; j<= Nstates:^2; j++) {
							if (VarAd[i,j] == 0 & VarAd[j,i] != 0) VarAd[i,j] = VarAd[j,i]
						}
					}
					
					// Calculate the variance
					tmp1 		= ((Imat + H)' # Imat) * varP * ((Imat + H) # Imat)
					tmp2 		= (Imat # P) * VarAd * (Imat # P')
					varP 		= tmp1 + tmp2 
				}
							
				
				// Matrix multiplication - calculate probability
				P = P*(Imat + H)
			}
				
			// Store results
			res_index = selectindex(t:==current_t :& d:==1 :& by:==b)
				
			if (ci | se) {
				seP = (rowshape(sqrt(diagonal(varP)),Nstates)[,fromS])'
				if (se) se_state[res_index,] = J(rows(res_index),1,seP[1,])
				if (ci) {			
					for(j=1;j<=Nstates;j++) {			
						P_colindex = (3*(j-1)+1)
						P_state[res_index,P_colindex..P_colindex+2] = J(rows(res_index),1,(P[fromS,j], P[fromS,j] :-1.96:*seP[1,j], P[fromS,j] :+1.96:*seP[1,j]))
					}
				}
			}
			if (ci != 1) P_state[res_index,] = J(rows(res_index),1,P[fromS,])
			
			if (los) {
				if (current_t != past_t) {
					P_t_unique[v] = current_t
					P_unique[v,] = P[fromS,]
					v = v + 1
				}
				past_t = current_t
			}
		}

		
		// LOS
		if (los) {
			
			for(u=1; u<=Nt; u++) {
				
				// Get current t and check before exit time and after enter time
				current_t = t_unique[u]
				if (current_t < entertime | current_t > exittime) continue
				
				// Work out expanded prob matrix
				if (P_t_unique[l] != .) {
					if (current_t == P_t_unique[l]) {
						P1[u,] = P_unique[l,]
						l = l + 1
					}
					else {
						if (l == 1) P1[u, ] = Imat[fromS,]
						else P1[u, ] = P_unique[l-1,]
					}
				}
				else {
					if (l == 1) P1[u, ] = Imat[fromS,]
					else P1[u, ] = P_unique[l-1,]
				}
				
				// Calculate LOS
				if (t0 == entertime) LOS = LOS :+ (current_t - t0):*Imat[fromS,]
				else LOS = LOS :+ (current_t - t0):*P1[u-1,]
				
				// Store LOS
				los_index = selectindex(t:==current_t :& by:==b)
				LOS_state[los_index,] = J(rows(los_index),1,LOS)		
				
				// Get t0 for next iteration
				t0 = current_t
			}
		}
	}
	
	
	// Check all probabilities within [0,1]
	P_state = 0:*(P_state:<0) :+ 1:*(P_state:>1) :+ P_state:*(P_state:>=0 :& P_state:<=1)

	// Output results
	newvars = tokens(st_local("newvars"))
	st_store(.,newvars,touse,P_state)
		
	if (se) {
		newvars_se = tokens(st_local("newvars_se"))
		st_store(.,newvars_se,touse,se_state)
	}
		
	if (los) {
		newvars_LOS = tokens(st_local("newvars_LOS"))
		st_store(.,newvars_LOS,touse,LOS_state)
	}
}


// Does checks and returns No. of transitions
function check_transmatrix(tmat)
{
	tmat_ind = tmat:!=.                                                  
	if (max(diagonal(tmat_ind))>0) {
		errprintf("All elements on the diagonal of transmatrix() must be coded missing = .\n")
		exit(198)
	}
	
	row = 1
	rtmat = rows(tmat)
	trans = 1
	while (row<rtmat) {
		for (i=1;i<=rtmat;i++) {
		
			if (sum(tmat:==tmat[row,i])>1 & tmat[row,i]!=.) {
				errprintf("Elements of transmatrix() are not unique\n")
				exit(198)
			}
			
			if (tmat[row,i]!=. & tmat[row,i]!=trans){
				errprintf("Elements of transmatrix() must be sequentially numbered from 1,...,K, where K = number of transitions\n")
				exit(198)
			}               
                        
			if (tmat[row,i]!=.) trans++
			}
                
			row++
        }
		
	return(trans-1)
}

		
// Creates to/from states for each transition
function transRowCol(tmat)
{
	tmat_index = J(max(tmat),2,.)
	row = 1
	rtmat = rows(tmat)
	trans = 1

	while (row<rtmat) {
		for (i=1;i<=rtmat;i++) {
			if(tmat[row,i] == trans) {
				tmat_index[trans,] = (row,i)
				trans++
			}
		}
		row++
	}
        return(tmat_index)
}
end	
