local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"
local RS real scalar
local RC real colvector

mata:


real colvector mm_root_vec(
						 transmorphic x,      // bj: will be replaced by solution
						 pointer(real scalar function) scalar f,
											  // Address of the function whose zero will be sought for
						 `RS' ax,      // Root will be sought for within a range [ax,bx]
						 `RS' bx,      //
						 | real scalar tol,   // Acceptable tolerance for the root value (default 0)
						   real scalar maxit, // bj: maximum # of iterations (default: 1000)
						   `opts')            // bj: additional args to pass on to f
{
    transmorphic  fs            // setup for f
    `RC'   a, b, c       // Abscissae, descr. see above
    //real scalar   fa, fb, fc    // f(a), f(b), f(c)
    real scalar   prev_step     // Distance from the last but one
    real scalar   tol_act       // Actual tolerance
    real scalar   p             // Interpolation step is calcu-
    real scalar   q             // lated in the form p/q; divi-
                                // sion operations is delayed
                                // until the last moment
    real scalar   new_step      // Step at this iteration
    real scalar   t1, cb, t2
    real scalar   itr

	`RS' nobs 
	`RC' index
	real colvector   fa, fb, fc    // f(a), f(b), f(c)
	
    if (args()<5) tol = 0       // bj: set tolerance
    if (args()<6) maxit = 1000  // bj: maximum # of iterations

    fs = mm_callf_setup(f, args()-6, `opts') // bj: prepare function call

    //x = .                       // bj: initialize output
	
	tempo1 = o1
	tempo6 = o6
	
	nobs = rows(x)
	index = 1::nobs
	result = J(nobs,1,.)	
	x = J(nobs,1,.)

    a = J(nobs,1,ax);  b = J(nobs,1,bx);  fa = mm_callf(fs, a);  fb = mm_callf(fs, b)
    c = a;  fc = fa

    //if ( fa==. ) return(0)      // bj: abort if fa missing
	tempindex = selectindex(fa:==.)
	nti = rows(tempindex)
	if (nti & cols(tempindex)) result[tempindex,] = J(nti,1,0)

	//remove tempindex as they are done 
	index = select(index,fa:!=.)

	if (rows(index)) { //not done
		tempindex = select(index,((fa[index]:>0) :* (fb[index]:>0))) 
		if (cols(tempindex)) {			
			
			flag1 = abs(fa[tempindex]) :< abs(fb[tempindex])
			flag2 = 1:-flag1
			tempindex2 = select(tempindex,flag1)
			nti = rows(tempindex2)
			if (nti) {
				result[tempindex2] = J(nti,1,2)
				x[tempindex2] = a[tempindex2]
			}
			tempindex2 = select(tempindex,flag2)
			nti = rows(tempindex2)
			if (nti) {
				result[tempindex2] = J(nti,1,3)
				x[tempindex2] = b[tempindex2]
			}
			//update index
			index = select(index,x:==.)					//better way to do this?
			if (rows(index)==0) return(result)
		}
		tempindex = select(index,((fa[index]:<0) :* (fb[index]:<0)))
		if (cols(tempindex)) {			

			flag1 = abs(fa[tempindex]) :< abs(fb[tempindex])
			flag2 = 1:-flag1
			tempindex2 = select(tempindex,flag1)
			nti = rows(tempindex2)
			if (nti) {
				result[tempindex2] = J(nti,1,2)
				x[tempindex2] = a[tempindex2]
			}
			tempindex2 = select(tempindex,flag2)
			nti = rows(tempindex2)
			if (nti) {
				result[tempindex2] = J(nti,1,3)
				x[tempindex2] = b[tempindex2]
			}
			//update index
			index = select(index,x:==.)					//better way to do this?
			if (rows(index)==0) return(result)
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
        b[index] = b[index] + new_step

		o1 = tempo1[index]
		o6 = tempo6[index]

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

real scalar mm_root2(
 transmorphic x,      // bj: will be replaced by solution
 pointer(real scalar function) scalar f,
                      // Address of the function whose zero will be sought for
 real scalar ax,      // Root will be sought for within a range [ax,bx]
 real scalar bx,      //
 | real scalar tol,   // Acceptable tolerance for the root value (default 0)
   real scalar maxit, // bj: maximum # of iterations (default: 1000)
   `opts')            // bj: additional args to pass on to f
{
    transmorphic  fs            // setup for f
    real scalar   a, b, c       // Abscissae, descr. see above
    real scalar   fa, fb, fc    // f(a), f(b), f(c)
    real scalar   prev_step     // Distance from the last but one
    real scalar   tol_act       // Actual tolerance
    real scalar   p             // Interpolation step is calcu-
    real scalar   q             // lated in the form p/q; divi-
								// sion operations is delayed
                                // until the last moment
    real scalar   new_step      // Step at this iteration
    real scalar   t1, cb, t2
    real scalar   itr

    if (args()<5) tol = 0       // bj: set tolerance
    if (args()<6) maxit = 1000  // bj: maximum # of iterations

    fs = mm_callf_setup(f, args()-6, `opts') // bj: prepare function call

    x = .                       // bj: initialize output

    a = ax;  b = bx;  
	timer_on(10)
	fa = mm_callf(fs, a);  
	fb = mm_callf(fs, b)
	timer_off(10)
    c = a;  fc = fa

    if ( fa==. ) return(0)      // bj: abort if fa missing

    if ( (fa > 0 & fb > 0) |    // bj: f(ax) and f(bx) do not
         (fa < 0 & fb < 0) ) {  //     have opposite signs
        if ( abs(fa) < abs(fb) ) {
            x = a; return(2)    // bj: fa closer to zero than fb
        }
        x = b; return(3)        // bj: fb closer to zero than fa
    }

    for (itr=1; itr<=maxit; itr++) {
        if ( fb==. ) return(0)            // bj: abort if fb missing

        prev_step = b-a

        if( abs(fc) < abs(fb) ) {         // Swap data for b to be the
            a = b;  b = c;  c = a;        // best approximation
            fa = fb;  fb = fc;  fc = fa
        }

        tol_act = 2*epsilon(b) + tol/2
        new_step = (c-b)/2

        if( abs(new_step) <= tol_act | fb == 0 ) {
             x = b                        // Acceptable approx. is found
             return(0)
        }

        // Decide if the interpolation can be tried
        if( abs(prev_step) >= tol_act     // If prev_step was large enough
             & abs(fa) > abs(fb) ) {      // and was in true direction,
                                          // Interpolation may be tried
            cb = c-b
            if( a==c ) {                  // If we have only two distinct
                t1 = fb/fa                // points linear interpolation
                p = cb*t1                 // can only be applied
                q = 1.0 - t1
            }
            else {                        // Quadric inverse interpolation
                q = fa/fc;  t1 = fb/fc;  t2 = fb/fa
                p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) )
                q = (q-1.0) * (t1-1.0) * (t2-1.0)
            }
            if( p>0 )                     // p was calculated with the op-
              q = -q                      // posite sign; make p positive
            else                          // and assign possible minus to
              p = -p                      // q
            if( p < (0.75*cb*q-abs(tol_act*q)/2) // If b+p/q falls in [b,c]
               & p < abs(prev_step*q/2) ) // and isn't too large
             new_step = p/q               // it is accepted
                                          // If p/q is too large then the
                                          // bisection procedure can
                                          // reduce [b,c] range to more
                                          // extent
        }

        if( abs(new_step) < tol_act ) {   // Adjust the step to be not less
            if( new_step > 0 )            // than tolerance
             new_step = tol_act
            else
             new_step = -tol_act
        }

        a = b;  fa = fb                   // Save the previous approx.
        b = b + new_step;  
		timer_on(10)
		fb = mm_callf(fs, b) // Do step to a new approxim.
		timer_off(10)
        if( (fb > 0 & fc > 0) | (fb < 0 & fc < 0) )  {
                      // Adjust c for it to have a sign opposite to that of b
            c = a;  fc = fa
        }
    }
    x = b
    return(1)                             // bj: convergence not reached
}

real colvector epsilon_vec(real colvector y)
{
	res = J(rows(y),1,.)
	for (i=1;i<=rows(y);i++) res[i] = epsilon(y[i])
	return(res)
}


end

