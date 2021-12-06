version 12.1

local SS 	string scalar
local RS	real scalar
local RM 	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct predictms_struct scalar) scalar

mata:

`RM' predictms_prob(	`PS' p,		///
						`RS' from,	///
						`RM' res,	///
						`RM' res2)
{
	struct predictms_struct scalar S	
	S = *p
	
	//for each state need to extract time at which that state was reached, or max cens time.
	//do for both res and res2, then compare
	ncols = cols(res)/2
	states1 = res[,1..ncols]
	stimes1 = res[,(ncols+1)..cols(res)]
	states2 = res2[,1..ncols]
	stimes2 = res2[,(ncols+1)..cols(res)]
	
	nvars = (S.Nstates-1)*2
	retres = J(S.obs,nvars,0)
	/*postind = 2

	for (i=from;i<=S.Nstates;i++) {
		//which column to extract times from
		ind = i+1					
		//if there's next states
		if (S.Nnextstates[i]) {
			posnextstates = asarray(S.posnextstatesj,i)		//extract state number of next states and loop
			for (j=1;j<=S.Nnextstates[i];j++) {
				//extract next state of interest, also using to post outcome
				col = posnextstates[j]										
				postind = (col-1) * 2 - 1

				//first set
				index = selectindex((states1[,ind]:==col) :* (states1[,ind-1]:!=col))
				temp1 = (stimes1[index,ind] :< stimes2[index,ind])
				for (t=1;t<=S.obs;t++) retres[t,postind] = retres[t,postind] :+ sum(temp1 :* (stimes1[index,ind] :< S.predtime[t])):/S.N

				//inverse probs
				postind++
				index = selectindex((states2[,ind]:==col) :* (states2[,ind-1]:!=col))
				temp1 = (stimes2[index,ind] :< stimes1[index,ind]) 
				for (t=1;t<=S.obs;t++) retres[t,postind] = retres[t,postind] :+ sum(temp1 :* (stimes2[index,ind] :< S.predtime[t])):/S.N
			}
		}
	}
	*/
	//need to extract the earliest time each state is reached
	//1,2,3
	estime1 = estime2 = J(S.N,S.Nstates,.)

	//loop over each state, extract first time that state is entered, then that ob is done
	for (i=2;i<=S.Nstates;i++) {
		index1 = index2 = S.coreindex
		for (j=2;j<=cols(stimes1);j++) {

			flag1 = states1[index1,j]:==i
			flag2 = 1:-flag1
			index12 = select(index1,flag1)	//these are in state i for first time

			flag3 = states2[index2,j]:==i
			flag4 = 1:-flag3

			index22 = select(index2,flag3)	//these are in state i for first time
			//post their times
			if (rows(index12) & cols(index12)) {
				estime1[index12,i] = stimes1[index12,j]
				//remove them from searching index	
				index1 = select(index1,flag2)
			}
			if (rows(index22) & cols(index22)) {
				estime2[index22,i] = stimes2[index22,j]
				//remove them from searching index	
				index2 = select(index2,flag4)
			}		
		}
	}

	ind1 = 1
	ind2 = 2
	for (i=2;i<=S.Nstates;i++) {
		temp1 = (estime1[,i] :< estime2[,i]) 
		temp2 = (estime2[,i] :< estime1[,i]) 
		for (t=1;t<=S.obs;t++) {
			temp3 = (estime1[,i] :< S.predtime[t]) :* (estime2[,i] :< S.predtime[t])
			retres[t,ind1] =  sum(temp1 :* temp3):/S.N
			retres[t,ind2] =  sum(temp2 :* temp3):/S.N
		}	
		ind1++
		ind2++
	}
	
	//retres[,1] = 1:-rowsum(retres[,2..nvars])
	return(retres)
}
end
