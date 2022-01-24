*! version 1.0.0 21Dec2020

program define freq_total
        version 14.2
        syntax /*[if][in]*/,                    [ID(varname)   ///
                                                Tmatrix(name)  ///
                                                STATESNames(string asis) ///
											    TRANSNames(string asis) ///
												TIMEPoints(varname)  ///
												SCale_freq(real 1)]
												
		
			
// check tranmatrix exists- create tmpstates, tmptrans

   if "`tmatrix'" != "" {
         
                qui summ _to, meanonly
                local tmpNstates `r(max)'
                tempname tmatrix
                matrix `tmatrix'= tmat
                       
   }       
   else {
                cap confirm matrix `tmatrix'
                if _rc>0 {
                        di as error "tmatrix(`tmatrix') not found"
                        exit 198
                }
   }
		
	
   local Nstates = rowsof(`tmatrix')

				
/**Find starting, intermediate and absorbing states***/

   forvalues i = 1/`Nstates' {
              local Absorbing`i' 1
			  local Starting`i'  1
     }

   qui summ _to
   
   local min_interm `r(min)'
		
   forvalues i = 1/`Nstates' {
			  local Starting`i'  0 
			  if `i'<`min_interm' {
					local Starting`i'  1
			  }
   }
		
   qui summ _trans
   
   local Ntransitions `r(max)'

   forvalues i = 1/`Nstates' {
	
			local Absorbing`i' 1
			forvalues j = 1/`Nstates' {
			
				if el(`tmatrix',`i',`j')!=.	{
					local Absorbing`i' 0
					continue, break
				}
			
	    }
   }


   forvalues i = 1/`Nstates' {
			local Intermediate`i'  0 
			
			if `Starting`i''==0 & `Absorbing`i''==0 {
					local Intermediate`i'  1 
	        }
   }
	

	/***********************************************************/
	
/*Looping over the timepoints computing the frequencies of states */

 // the timepoints for which we compute the freuence of people in each state and transition

   tempvar id2
   gen `id2'=`id'
   

   qui levelsof `timepoints', local(levels) 
   global times "`levels'"

	 
/*****Store the frequency results in a matrix- Creating the dimensions of the matrix*****/

   qui summ `timepoints'
   local tmpNtimepoints =`r(N)'

	
   qui summ _trans
   local Ntransitions `r(max)'


   local tmpNcol= 1+`Nstates'+`Ntransitions'+1
	
   tempvar time_label
	
   egen `time_label'= seq() , f(1) t(`tmpNtimepoints')
	

   tempname freq_matrix
   matrix `freq_matrix' = J(`tmpNtimepoints',`tmpNcol',.)
		
   local	p=1
   
/******************************************************************************************/
// Computing the frequencies in each state and transition for up to each specified timepoint
/******************************************************************************************/

   foreach k in $times { 

		local namek : subinstr local k "." ""

        forvalues i = 1/`Nstates' { 

	  
		     if `Starting`i''==1 {
		     	   tempvar x1
		     	   gen `x1'=0
		     	   qui bysort `id' (_trans) :  replace `x1'= 1 if _from==`i' & _status==1 & _stop<=`k'/`scale_freq'
		     	   tempvar x2
		     	   qui bysort `id' :  egen `x2'=total (`x1') 
		     	   tempvar firstobs
		     	   qui bysort `id': gen `firstobs' = _n==1
		     	   qui count if `firstobs' & `x2'==0
			
				   local freq_`i'_`namek' = `r(N)'
				   drop `x1' `x2'
			
	        }



		    if `Intermediate`i''==1 {
		     	   tempvar x1
		     	   gen `x1'=0
		            qui bysort `id' (_trans) :  replace `x1'= 1 if _from==`i' & _status==1 & _stop<=`k'/`scale_freq'
		     	   tempvar x2
		     	   qui bysort `id' :  egen `x2'=total (`x1') 
		     	   
		     	   tempvar x3
		     	   gen `x3'=0	
		     	   qui bysort `id' (_trans) :  replace `x3'= 1 if _from!=`i'  & _to==`i' & _status==1 & _stop<=`k'/`scale_freq'
		     	   tempvar x4
		     	   qui bysort `id' :  egen `x4'=total (`x3') 
		     	   
		     	   tempvar firstobs
		     	   qui bysort `id': gen `firstobs' = _n==1
		     	   qui count if `firstobs' & `x2'==0 & `x4'!=0
		     	   
		     	   
		     	   local freq_`i'_`namek' = `r(N)'
        
		     	   drop `x1' `x2' `x3' `x4'
	    
		   }

	
		   if `Absorbing`i''==1 {
  
  			      tempvar x1
			      gen `x1'=0
			      qui bysort `id' (_trans):  replace `x1'= 1 if  _to==`i' & _status==1 & _stop<=`k'/`scale_freq'
			      tempvar x2
			      qui bysort `id' :  egen `x2'=total (`x1') 
			      tempvar firstobs
			      qui bysort `id': gen `firstobs' = _n==1
			      qui count if `firstobs' & `x2'!=0

				  local freq_`i'_`namek' = `r(N)'
				  drop `x1' `x2'

		  }
	
       }



	/*Looping over the timepoints computing the frequencies of  transitions*/

       qui summ _trans
       local Ntransitions `r(max)'

       forvalues i = 1/`Ntransitions' {
  
    		tempvar x1
			gen `x1'=0
			qui bysort `id' :  replace `x1'= 1 if _trans==`i' & _status==1 & _stop<=`k'/`scale_freq'
			tempvar x2
			qui bysort `id' :  egen `x2'=total(`x1') 
			tempvar firstobs
			qui bysort `id': gen `firstobs' = _n==1
			qui count if `firstobs' & `x2'!=0
			
			local freq_trans_`i'_`namek' = `r(N)'

			drop `x1' `x2'
			
			
	   }

	   
// Building of the frequency matrix

		matrix `freq_matrix'[`p', 1] = `p'
			
	    forvalues i = 1/`Nstates' {	
                        matrix `freq_matrix'[`p', 1+`i'] =`freq_`i'_`namek''
        }
			  
		forvalues j = 1/`Ntransitions' {
                            matrix `freq_matrix'[`p', 1+`Nstates' +`j'] =  `freq_trans_`j'_`namek''
        }	
			
		matrix `freq_matrix'[`p',  1+`Nstates' +`Ntransitions'+1] = `timepoints'[`p']	
			
		local p= `p'+1
		
		
   }
   
/***************************************************************/



// Saving the frequancy matrix 
   matrix def frequencies=  `freq_matrix'	
   
   
// Give names to the columns of the frequency matrix

   matname frequencies time_label , columns(1) explicit

   local l=2
   forvalues i = 1/`Nstates' {
			matname frequencies  State`i'  ,  columns(`l') explicit
			local l=`l'+1
   }
  
    

   if `"`transnames'"' == "" {	
			local l=1+`Nstates'+1
			forvalues i = 1/`Ntransitions' {
				matname frequencies  h`i'  ,  columns(`l') explicit
				local l=`l'+1
			}
   }
   
   
    if `"`transnames'"' != "" {
			local L0=1+`Nstates'+1
			local l=1+`Nstates'+`Ntransitions'
			matname frequencies  "`transnames'"   ,  columns(`L0'..`l') explicit
	}
	
	local L_final=1+`Nstates'+`Ntransitions'+1
	matname frequencies timevar   ,  columns(`L_final') explicit
		
end
