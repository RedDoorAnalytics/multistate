local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

version 14.2

local ss 	string scalar
local RS	real scalar
local NM	numeric matrix
local RM 	real matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct predictms_struct scalar) scalar
local SS	struct predictms_struct scalar
local gml 	struct merlin_struct scalar

mata:

`RC' predictms_sim_root(transmorphic x,      						/// bj: will be replaced by solution
        pointer(real matrix function) scalar f,		/// Address of the function whose zero will be sought for
        `RS' bx,									/// Root will be sought for within a range [t0,bx]
        `RS' tol,   								/// Acceptable tolerance for the root value (default 0)
        `RS' maxit,									/// maximum # of iterations (default: 1000)
        `RS' Nsim,|									///
        `opts')            							//  additional args to pass on to f
{
    transmorphic  fs    // setup for f
    `RC' a, b, c       	// Abscissae, descr. see above
    `RS' fa, fb, fc  	// f(a), f(b), f(c)
    `RS' prev_step     	// Distance from the last but one
    `RS' tol_act       	// Actual tolerance
    `RS' p             	// Interpolation step is calculated in the form p/q; division operations is delayed until the last moment
    `RS' q             	// 
    `RS' new_step      	// Step at this iteration
    `RS' t1, cb, t2
    `RS' itr
    `RC' index

    fs = mm_callf_setup(f, args()-6, `opts') 	// bj: prepare function call
	index	= o3
	result 	= J(Nsim,1,.)	

	//inputs to update with indexing
	tempo1 	= o1					//t0
	tempo2 	= o2					//logU
	
	//limits
        a = o1
	b = J(Nsim,1,bx)

	fa = mm_callf(fs, a);  fb = mm_callf(fs, b)
        c = a;  fc = fa

    //if ( fa==. ) return(0)      // bj: abort if fa missing
	tempindex = selectindex(fa:==.)
	nti = rows(tempindex)
	if (nti & cols(tempindex)) result[tempindex,] = J(nti,1,0)

	//remove tempindex as they are done 
	index = select(index,fa:!=.)

	if (rows(index) & cols(index)) { //not done

		tempindex = select(index,((fa[index]:>0) :* (fb[index]:>0))) 
		if (rows(tempindex) & cols(tempindex)) {			
			
			flag1 = abs(fa[tempindex]) :< abs(fb[tempindex])
			flag2 = 1:-flag1
			tempindex2 = select(tempindex,flag1)
			nti = rows(tempindex2)
			if (nti & cols(tempindex2)) {
				result[tempindex2] = J(nti,1,2)
				x[tempindex2] = a[tempindex2]
			}
			tempindex2 = select(tempindex,flag2)
			nti = rows(tempindex2)
			if (nti & cols(tempindex2)) {
				result[tempindex2] = J(nti,1,3)
				x[tempindex2] = b[tempindex2]
			}
			//update index
			index = select(index,x[index]:==.)	
			if (!rows(index) & !cols(index)) return(result)
		}

		tempindex = select(index,((fa[index]:<0) :* (fb[index]:<0)))
		if (rows(tempindex) & cols(tempindex)) {			

			flag1 = abs(fa[tempindex]) :< abs(fb[tempindex])
			flag2 = 1:-flag1
			tempindex2 = select(tempindex,flag1)
			nti = rows(tempindex2)
			if (nti & cols(tempindex2)) {
				result[tempindex2] = J(nti,1,2)
				x[tempindex2] = a[tempindex2]
			}
			
			tempindex2 = select(tempindex,flag2)
			nti = rows(tempindex2)
			if (nti & cols(tempindex2)) {
				result[tempindex2] = J(nti,1,3)
				x[tempindex2] = b[tempindex2]
			}
			//update index
			index = select(index,x[index]:==.)	
			if (!rows(index) & !cols(index)) return(result)
		}
	}
	else return(result)

	for (itr=1; itr<=maxit; itr++) {

		tempindex = index[selectindex(fb[index]:==.)]

		if (cols(tempindex)) result[tempindex] = J(rows(tempindex),1,0)

		//remove tempindex as they are done 
		index = select(index,fb[index]:!=.)

		if (!cols(index)) return(result)

		tempindex = select(index,abs(fc[index]) :< abs(fb[index]))

		if (cols(tempindex)) {
			a[tempindex] = b[tempindex];  b[tempindex] = c[tempindex];  c[tempindex] = a[tempindex];         // best approximation
            fa[tempindex] = fb[tempindex];  fb[tempindex] = fc[tempindex];  fc[tempindex] = fa[tempindex]
		}

		tol_act = 2:*epsilon_vec(b[index]) :+ tol:/2
                new_step = (c[index]:-b[index]):/2

		flag1 = (abs(new_step):<=tol_act) :+ (fb[index]:==0)
		flag2 = (flag1:==0)
		tempindex = select(index,flag1)

		if (cols(tempindex)) {
			x[tempindex] = b[tempindex]
			result[tempindex] = J(rows(tempindex),1,0)
		}

		index = select(index,flag2)  
		if (!cols(index) | !rows(index)) return(result)

		//update stuff
		tol_act = select(tol_act,flag2)
		new_step = select(new_step,flag2)

        // Decide if the interpolation can be tried
		prev_step = b[index]:-a[index]

		tempindex11 = (abs(prev_step) :>= tol_act) :* (abs(fa[index]) :> abs(fb[index]))
		tempindex = select(index,tempindex11)

		if (cols(tempindex)) {
		
			cb = c[tempindex] :- b[tempindex]
			
			p = q  = cb:*0					//fix

			flag1 = a[tempindex] :== c[tempindex]
			flag2 = 1:-flag1
			tempindex2 = select(tempindex,flag1)
			if (cols(tempindex2)) {
				t1 = fb[tempindex2]:/fa[tempindex2]
				p[selectindex(flag1)] = select(cb,flag1) :* t1
				q[selectindex(flag1)] = 1:- t1
			}
			tempindex2 = select(tempindex,flag2)
			if (cols(tempindex2)) {			
				q[selectindex(flag2)] = fa[tempindex2]:/fc[tempindex2]; t1 = fb[tempindex2]:/fc[tempindex2]; t2 = fb[tempindex2]:/fa[tempindex2]
				p[selectindex(flag2)] = t2 :* ( select(cb,flag2) :* q[selectindex(flag2)] :* (q[selectindex(flag2)] :- t1) :- (b[tempindex2]:-a[tempindex2]):*(t1:-1) )
                q[selectindex(flag2)] = (q[selectindex(flag2)]:-1) :* (t1:-1) :* (t2:-1)
			}
			flag1 = p:>0
			flag2 = 1:-flag1
			tempindex = selectindex(flag1)
			if (cols(tempindex)) q[tempindex] = -q[tempindex]
			tempindex = selectindex(flag2)
			if (cols(tempindex)) p[tempindex] = -p[tempindex]

			tempindex = (p :< (0.75:*cb:*q:-abs(select(tol_act,tempindex11):*q):/2))  :* (p :< abs(select(prev_step,tempindex11):*q:/2))
			if (cols(tempindex)) {
				//update tempindex11
				tempindex22 = select(selectindex(tempindex11),tempindex)
				if (cols(tempindex22) & rows(tempindex22)) new_step[tempindex22] = p[selectindex(tempindex)]:/q[selectindex(tempindex)]
			}
			
		}

		tempindex = selectindex(abs(new_step) :< tol_act)
		if (rows(tempindex)) {
			flag1 = new_step[tempindex] :> 0
			flag2 = 1:-flag1
			tempindex2 = select(tempindex,flag1)
			if (rows(tempindex2)) new_step[tempindex2] = tol_act[tempindex2]
			tempindex2 = select(tempindex,flag2)
			if (rows(tempindex2)) new_step[tempindex2] = -tol_act[tempindex2]
        }

        a[index] = b[index];  fa[index] = fb[index]                   // Save the previous approx.
        b[index] = b[index] :+ new_step

		o1 = tempo1[index]
		o2 = tempo2[index]	
		//update indexing
		o3 = index

		fb[index] = mm_callf(fs, b[index]) // Do step to a new approxim.

		tempindex1 = select(index,((fb[index]:>0) :* (fc[index]:>0)))
		tempindex2 = select(index,((fb[index]:<0) :* (fc[index]:<0)))

		if (cols(tempindex1)) {			
			c[tempindex1] = a[tempindex1]
			fc[tempindex1] = fa[tempindex1]
		}
		if (cols(tempindex2)) {			
			c[tempindex2] = a[tempindex2]
			fc[tempindex2] = fa[tempindex2]
		}
    }

	x[index] = b[index]
	result[index] = J(rows(index),1,0)
    return(result)                             // bj: convergence not reached
}

real colvector epsilon_vec(real colvector y)
{
	res = J(rows(y),1,.)
	for (i=1;i<=rows(y);i++) res[i] = epsilon(y[i])
	return(res)
}

end
