*! version 1.0.0 26jan2021 MJC

/*
History
NS 26jan2021 version 1.0.0 - lines 30, 34,35, 50, 41 adding extra options to the msboxes command
						   - line 145 Make sure that if someone defines the transnames they give length equal to Ntransitions
						   - line 187 Default transition Names for the frequency matrix if none are defined by the user during the msboxes command
						   - line 242 create locals Trans`i' with the  transition names
						   - lines 291, 310, 335, 353, 363  inserting the information about the coordinates of the transition names 
						   - line 398 addition of the local arrowtext_trans in the scatterplot (coordinates local)
						   - line 419-432 Error messages: If freqat option is not defined or not a variable or negative we do not invoke freq_total command and thus do not ceate the frequency matrix
						   - line 436  If freqat option is defined then the freq_total command is run and the frequency matrix is created
						   - update WriteJson() mata function so that it puts the frequency matrix in the json file only when this matrix exists (a.k.a use of freqat option)
PL 30dec2017 version 0.1.0
*/

program define msboxes, rclass
        version 14.2
        syntax /*[if][in]*/, 			Xvalues(numlist)          	///
										Yvalues(numlist)          	///
										ID(varname)               	///
									[								///
										TRANSMATrix(name)        	///
										STATENames(string asis)  	///
										TRANSNames(string asis)   	/// transition names option
										BOXHeight(real 0.3)       	///
										BOXWidth(real 0.2)        	///
										YRANge(string)            	///
										SCale(real 1)             	/// scale
										FREQat(varname)           	/// time points we want frequency for
										XRANge(string)            	///
										YSIZE(string)             	///
										XSIZE(string)             	///
										GRid                      	///
										INTERactive               	///  interactive option    
										JSONPath(string)          	///  jsonpath
										CR                        	///
									]								//
        //marksample touse
        
        tempvar endstate
	
        // check tranmatrix exists
        if "`cr'" != "" {
                if "`transmatrix'" != "" {
                        di "do not specify both the transmatrix and cr option"
                        exit 198
                }
                summ _to, meanonly
                local tmpNstates `r(max)'
                tempname transmatrix
                matrix `transmatrix' = J(`tmpNstates',`tmpNstates',.)
                forvalues i = 2/`tmpNstates' {
                        local tmptrans = `i' - 1
                        matrix `transmatrix'[1,`i'] = `tmptrans'
                }

        }       
        else {
                cap confirm matrix `transmatrix'
                if _rc>0 {
                        di as error "transmatrix(`transmatrix') not found"
                        exit 198
                }
        }
        
        local Nstates = rowsof(`transmatrix')
        
        // Checks for interactive options
        if "`jsonpath'" != "" & "`interactive'" == "" {
                di as error "You have used the jsonpath option without using the interactive option."
                exit 198
        }
        if "`interactive'" != "" {
                if "`jsonpath'" != "" {
                        mata st_local("direxists",strofreal(direxists("`jsonpath'")))
                        if !`direxists' {
                                di as error "Folder `jsonpath' does not exist."
                                exit 198
                        }
                        mata st_local("jsonfile",pathjoin("`jsonpath'","msboxes.json"))
                }
                else {
                        local jsonfile msboxes.json
                }
                capture confirm file "`jsonfile'"
                if !_rc {
                        capture erase "`jsonfile'"
                        if _rc {
                                display as error "`jsonfile' cannot be deleted'"
                        }
                }
        }
         
        if "`yrange'" == "" {
                local ymin 0
                local ymax 1
        }
        else {
                numlist "`yrange'", ascending min(2) max(2)
                local ymin = word("`r(numlist)'",1)
                local ymax = word("`r(numlist)'",2)     
        }
        if "`xrange'" == "" {
                local xmin 0
                local xmax 1
        }
        else {
                numlist "`xrange'", ascending min(2) max(2)
                local xmin = word("`r(numlist)'",1)
                local x max = word("`r(numlist)'",2)    
        }       
        if "`ysize'" != "" local ysize ysize(`ysize')
        if "`xsize'" != "" local xsize xsize(`xsize')
        

        // Check data seems to be stset correctly
        confirm var _trans 
        confirm var _status
        confirm var _start
        confirm var _stop
        
        // Xvalues and Yvalues
        if (wordcount("`xvalues'") != `Nstates') | (wordcount("`yvalues'") != `Nstates') {
                di as error "xvalues and yvalues must be of length `Nstates' (the number of states)"
                exit 198
        }
        local wc : word count `statenames'
        if `"`statenames'"' != "" & `wc' != `Nstates' {
                di as error "Length of statenames option should be equal to `Nstates'" ///
                                        "(The number of states)"
                exit 198
        }
		
		
		
		
		
		//Change 16/12/2020
		 // Make sure that if someone defines the transnames they give length equal to Ntransitions
		qui summ _trans
        local Ntransitions `r(max)'
		local wc_trans : word count `transnames'
        if `"`transnames'"' != "" & `wc_trans' != `Ntransitions' {
                di as error "Length of transnames option should be equal to `Ntransitions'" ///
                                        "(The number of transitions)"
                exit 198
        }

// ============================================================================ 
//      set up 
// ============================================================================ 

    //    qui summ _trans
    //    local Ntransitions `r(max)'

        forvalues i = 1/`Nstates' {
                local Absorbing`i' 1
                forvalues j = 1/`Nstates' {
                        if el(`transmatrix',`i',`j') != . {
                                local Absorbing`i' 0
                                continue, break
                        }
                }
        }
        
// end state (assumes some ordering)
        tempvar maxeventtime totalmiss
        qui bysort `id' (_trans): egen `maxeventtime' = max(_stop*_status )
        qui bysort `id' (_trans): gen `endstate' = _to if _stop == `maxeventtime' & _status==1 
        qui bysort `id' (_trans): egen `totalmiss' = total(`endstate'==.)
        qui bysort `id' (_trans): replace `endstate' = _from if `totalmiss' == _N & _n==1

// Default State Names
        if `"`statenames'"' == "" {
                forvalues i = 1/`Nstates' {
                        local statenames `statenames' State_`i'
                }
        }
		
 //Change 16/12/2020
// Default transition Names for json file frequency matrix
/*check it out*/ 

        if `"`transnames'"' == "" {
                forvalues i = 1/`Ntransitions' {
                        local transnames `transnames' h`i'
                }
        }
//		

        
// Start and end state for each transition
        forvalues i=1/`Nstates' {
                forvalues j=1/`Nstates' {
                        if el(`transmatrix',`i',`j') != . {
                                local tmptrans = el(`transmatrix',`i',`j')
                                local trans`tmptrans'_start `i'
                                local trans`tmptrans'_stop `j'
                        }
                }
        }
        

// ============================================================================
// Do calculations
// ============================================================================ 

// Number at risk at start
        tempvar fromfirst
        bysort `id' _from: gen `fromfirst' = _n==1
        forvalues i = 1/`Nstates' {
                qui count if _from == `i' & `fromfirst' & !`Absorbing`i''
                local Nstart`i' `r(N)'
                qui count if `endstate' == `i' 
                local Nend`i' `r(N)'
        }
// Number of subjects transitioning
        forvalues i = 1/`Ntransitions' {
                qui count if _trans==`i' & _status == 1 
                local Nevents`i' `r(N)'
        }
        
// ============================================================================
// Draw Graph
// ============================================================================ 

// state names
        forvalues i = 1/`Nstates' {
                local x = word("`xvalues'",`i')
                local y = real(word("`yvalues'",`i')) + `boxheight'/3 
                local text: word `i' of `statenames'
                local statetext `statetext' text(`y' `x' "`text'", placement(c))
        }

// Change 16/12/2020 create locals Trans`i' with the  transition names
        if `"`transnames'"' != "" {
		    forvalues i = 1/`Ntransitions' {
		
				   local text_trans: word `i' of `transnames'
                   local Trans`i'  "`text_trans'"		
            }
        }
//

// boxes
        forvalues i = 1/`Nstates' {
                local x1 = real(word("`xvalues'",`i')) - `boxwidth'/2
                local x2 = real(word("`xvalues'",`i')) + `boxwidth'/2
                local y1 = real(word("`yvalues'",`i')) - `boxheight'/2          
                local y2 = real(word("`yvalues'",`i')) + `boxheight'/2          
                local boxes `boxes' (pci `y1' `x1' `y1' `x2' `y1' `x1' `y2' `x1' `y2' `x1' `y2' `x2' `y2' `x2' `y1' `x2', lcolor(black))
        }
        
// Number at Start
        forvalues i = 1/`Nstates' {
                local x = real(word("`xvalues'",`i')) - `boxwidth'/2 + 0.01
                local y = real(word("`yvalues'",`i')) - `boxheight'/3 
                local text = "`Nstart`i''"
                local Nstarttext `Nstarttext' text(`y' `x' "`text'", placement(e))
        }
// Number at End
        forvalues i = 1/`Nstates' {
                local x = real(word("`xvalues'",`i')) + `boxwidth'/2 -0.01
                local y = real(word("`yvalues'",`i')) - `boxheight'/3 
                local text = "`Nend`i''"
                local Nendtext `Nendtext' text(`y' `x' "`text'", placement(w))
        }       
                
// arrows & arrow text
        forvalues i = 1/`Ntransitions' {
                // horizontal
                if      word("`yvalues'",`trans`i'_start') == word("`yvalues'",`trans`i'_stop') {
                        local lefttoright = real(word("`xvalues'",`trans`i'_start')) < real(word("`xvalues'",`trans`i'_stop'))
                        local y1 = real(word("`yvalues'",`trans`i'_start'))
                        local x1 = real(word("`xvalues'",`trans`i'_start')) + (cond(`lefttoright',1,-1)*`boxwidth'/2)
                        local y2 = real(word("`yvalues'",`trans`i'_stop')) 
                        local x2 = real(word("`xvalues'",`trans`i'_stop')) + (cond(`lefttoright',-1,1)*`boxwidth'/2)
                        
                        local ytext = real(word("`yvalues'",`trans`i'_start')) + 0.01
                        local xtext = (real(word("`xvalues'",`trans`i'_start')) + real(word("`xvalues'",`trans`i'_stop')))/2
        
                        local arrowtext `arrowtext' text(`ytext' `xtext' "`Nevents`i''", placement(n))
						
						// Change 16/12/2020 add transition names
						local ytext_trans = real(word("`yvalues'",`trans`i'_start')) - 0.01
                        local xtext_trans = (real(word("`xvalues'",`trans`i'_start')) + real(word("`xvalues'",`trans`i'_stop')))/2
                        local arrowtext_trans `arrowtext_trans' text(`ytext_trans' `xtext_trans' "`Trans`i''", placement(s) )
						//

                } 
                else if word("`xvalues'",`trans`i'_start') == word("`xvalues'",`trans`i'_stop') {
                        local toptobottom = real(word("`yvalues'",`trans`i'_start')) > real(word("`yvalues'",`trans`i'_stop'))
                        local y1 = real(word("`yvalues'",`trans`i'_start')) + (cond(`toptobottom',-1,1)*`boxheight'/2)
                        local x1 = real(word("`xvalues'",`trans`i'_start'))
                        local y2 = real(word("`yvalues'",`trans`i'_stop')) + (cond(`toptobottom',1,-1)*`boxheight'/2)
                        local x2 = real(word("`xvalues'",`trans`i'_stop')) 
                        
                        local ytext = (real(word("`yvalues'",`trans`i'_start')) + real(word("`yvalues'",`trans`i'_stop')))/2
                        local xtext = (real(word("`xvalues'",`trans`i'_start')) -0.01)
        
                        local arrowtext `arrowtext' text(`ytext' `xtext' "`Nevents`i''", placement(w))
						
						// Change 16/12/2020 add transition names
						local ytext_trans = (real(word("`yvalues'",`trans`i'_start')) + real(word("`yvalues'",`trans`i'_stop')))/2
                        local xtext_trans = (real(word("`xvalues'",`trans`i'_start')) + 0.01)
                        local arrowtext_trans `arrowtext_trans' text(`ytext_trans' `xtext_trans' "`Trans`i''", placement(e))
						//
                        
                }
                else {
                        local cutoff 0.5
                        local gradient =  abs(real(word("`yvalues'",`trans`i'_start')) - real(word("`yvalues'",`trans`i'_stop'))) ///
                                                        / abs(real(word("`xvalues'",`trans`i'_start')) - real(word("`xvalues'",`trans`i'_stop')))
                        
                        local textleft = (real(word("`xvalues'",`trans`i'_start')) > real(word("`xvalues'",`trans`i'_stop')))                           
                        if      `gradient' < `cutoff' {
                                local lefttoright = real(word("`xvalues'",`trans`i'_start')) < real(word("`xvalues'",`trans`i'_stop'))
                                local y1 = real(word("`yvalues'",`trans`i'_start'))
                                local x1 = real(word("`xvalues'",`trans`i'_start')) + (cond(`lefttoright',1,-1)*`boxwidth'/2)
                                local y2 = real(word("`yvalues'",`trans`i'_stop')) 
                                local x2 = real(word("`xvalues'",`trans`i'_stop')) + (cond(`lefttoright',-1,1)*`boxwidth'/2)
                                
                                local ytext = (real(word("`yvalues'",`trans`i'_start')) + real(word("`yvalues'",`trans`i'_stop')))/2
                                local xtext = (real(word("`xvalues'",`trans`i'_start')) + real(word("`xvalues'",`trans`i'_stop')))/2 + cond(`textleft',-1,1)*0.02
        
                                local arrowtext `arrowtext' text(`ytext' `xtext' "`Nevents`i''", placement(`=cond(`textleft',"w","e")'))
								
								// Change 16/12/2020 add transition names
								local ytext_trans = (real(word("`yvalues'",`trans`i'_start')) + real(word("`yvalues'",`trans`i'_stop')))/2
                                local xtext_trans = (real(word("`xvalues'",`trans`i'_start')) + real(word("`xvalues'",`trans`i'_stop')))/2 + cond(`textleft',-1,1)*(-0.02)
                                local arrowtext_trans `arrowtext_trans' text(`ytext_trans' `xtext_trans' "`Trans`i''", placement(`=cond(`textleft',"e","w")'))
								//
                        }
                        else {
                                local toptobottom = real(word("`yvalues'",`trans`i'_start')) > real(word("`yvalues'",`trans`i'_stop'))
                                local y1 = real(word("`yvalues'",`trans`i'_start')) + (cond(`toptobottom',-1,1)*`boxheight'/2)
                                local x1 = real(word("`xvalues'",`trans`i'_start'))
                                local y2 = real(word("`yvalues'",`trans`i'_stop')) + (cond(`toptobottom',1,-1)*`boxheight'/2)
                                local x2 = real(word("`xvalues'",`trans`i'_stop')) 

                                local ytext = (real(word("`yvalues'",`trans`i'_start')) + real(word("`yvalues'",`trans`i'_stop')))/2
                                local xtext = (real(word("`xvalues'",`trans`i'_start')) + real(word("`xvalues'",`trans`i'_stop')))/2 + (cond(`textleft',-1,1)*0.02)
        
                                local arrowtext `arrowtext' text(`ytext' `xtext' "`Nevents`i''", placement(`=cond(`textleft',"w","e")'))
							    
								// Change 16/12/2020 add transition names
                                local ytext_trans = (real(word("`yvalues'",`trans`i'_start')) + real(word("`yvalues'",`trans`i'_stop')))/2
                                local xtext_trans = (real(word("`xvalues'",`trans`i'_start')) + real(word("`xvalues'",`trans`i'_stop')))/2 + (cond(`textleft',-1,1)*(-0.02))
                                local arrowtext_trans `arrowtext_trans' text(`ytext_trans' `xtext_trans' "`Trans`i''", placement(`=cond(`textleft',"e","w")'))
                                //
                        }
                }
                local arrowstextx `arrowstextx' `xtext'
                local arrowstexty `arrowstexty' `ytext'
				
				// Change 16/12/2020 add transition names
				local arrowtext_transx  `arrowtext_transx' `xtext'
				local arrowtext_transy  `arrowtext_transy' `ytext'
                //
				
                local arrows `arrows' (pcarrowi `y1' `x1' `y2' `x2', color(red) barbsize(1) msize(2))
                local arrowsx1 `arrowsx1' `x1' 
                local arrowsy1 `arrowsy1' `y1' 
                local arrowsx2 `arrowsx2' `x2' 
                local arrowsy2 `arrowsy2' `y2' 
        }
// arrow text


// Grid option
        if "`grid'" == "" {
             local xlab xlab(minmax,nolabels noticks nogrid) 
             local ylab ylab(minmax,nolabels noticks nogrid) 
        }
        else {
             local xlab xlab(,grid  glcolor(gs7)) 
             local ylab ylab(,grid  glcolor(gs7))
        }


// Now plot everything
        twoway (scatteri `ymin' `xmin' `ymax' `xmax', msymbol(none)) ///
                        `boxes' ///
                        `arrows', ///
                        xscale(off range(`xmin' `xmax')) ///
                        yscale(off range(`ymin' `ymax')) ///
                        `statetext' ///
                        `Nstarttext' ///
                        `Nendtext' ///
                        `arrowtext' ///
						`arrowtext_trans' ///  Change 16/12/2020 add transition names
                        `ylab' ///
                        `xlab' ///
                        graphregion(color(white) margin(0 0 0 0)) bgcolor(white) ///
                        plotregion(margin(0 0 0 0)) ///
                        `ysize' `xsize' ///
                        legend(off)


	   tempvar id_inner
	   gen `id_inner'=`id'

	   tempname trmat
	   mat define `trmat'=`transmatrix'

	   local states `"`statenames'"'
     
	   local trans `"`transnames'"'
	  
	   local scale_inner =`scale'

//Change 16/12/2020 : Only if freqat is defined then we invoke the freq_total command and ceate the frequency matrix
// Use varname in the option 
      if "`freqat'"!="" {
                        
                cap confirm numeric variable `freqat'
                if _rc {
                         di as error "Invalid freqat()"
                         exit 198
                }
               
                if `freqat'<0 {
                         di as error "freqat() must be >0"
                         exit 198
                }
              
                        

				freq_total, id(`id_inner') tmatrix(`trmat') ///
						statesnames(`states') transnames(`trans') ///
						timepoints(`freqat') scale_freq(`scale_inner')
	
				local nrow=rowsof(frequencies)
	 }      
//
	
	mat def tmat2=`trmat'     

				
	// save to JSON file
	if `"`interactive'"' != "" {
			mata: WriteJson()
	}

	if "`freqat'"!="" {
		return matrix frequencies = frequencies
	}

end
/***************************************************************************************************/

program closeallfiles
        forvalues i=0(1)50 {
                capture mata: fclose(`i')
        }
end

// Updated 16/12/2020 WriteJson() to not produce the frequency matrix when freqat is not defined
mata
function WriteJson() {
        filename = st_local("jsonfile")
		
		hasfreq=st_local("freqat")!=""
		
		Nstates = st_local("Nstates")
		Ntransitions =st_local("Ntransitions")
		boxwidth = strtoreal(st_local("boxwidth"))
        boxheight = strtoreal(st_local("boxheight"))
        xvalues = invtokens(strofreal(strtoreal(tokens(st_local("xvalues"))) :- 0.5:*boxwidth, "%9.5f"),",")
        yvalues = invtokens(strofreal(strtoreal(tokens(st_local("yvalues"))) :+ 0.5:*boxheight,"%9.5f"),",")
        ymin = st_local("ymin")
        ymax = st_local("ymax")
        xmin = st_local("xmin")
        xmax = st_local("xmax")
        statenames = tokens(st_local("statenames"))
        statenames = invtokens(`"""' :+ statenames :+ `"""',",")
		transnames = tokens(st_local("transnames"))
        transnames = invtokens(`"""' :+ transnames :+ `"""',",")
		tmat = st_matrix("tmat2")
		
		
	    if (hasfreq==1) {
				frequencies = st_matrix("frequencies")
				clabels = st_matrixcolstripe("frequencies")
				nrow = strtoreal(st_local("nrow"))
		}
        arrowsx1 = "["+invtokens(strofreal(strtoreal(tokens(st_local("arrowsx1"))),"%9.5f"),",")+"]"
        arrowsy1 = "["+invtokens(strofreal(strtoreal(tokens(st_local("arrowsx2"))),"%9.5f"),",")+"]"
        arrowsx2 = "["+invtokens(strofreal(strtoreal(tokens(st_local("arrowsy1"))),"%9.5f"),",")+"]"
        arrowsy2 = "["+invtokens(strofreal(strtoreal(tokens(st_local("arrowsy2"))),"%9.5f"),",")+"]"
        arrows = `"{"x1":"' + arrowsx1 + `","y1":"'+arrowsy1+`","x2":"'+arrowsx2 +`","y2":"'+arrowsy2+"}"
        arrowstextx = "["+invtokens(strofreal(strtoreal(tokens(st_local("arrowstextx"))),"%9.5f"),",")+"]"
        arrowstexty = "["+invtokens(strofreal(strtoreal(tokens(st_local("arrowstexty"))),"%9.5f"),",")+"]"
        arrowstext = `"{"x":"' + arrowstextx + `","y":"'+arrowstexty+"}"
		//filename
		
		NstatesR=strtoreal(Nstates)
		NtransitionsR=strtoreal(Ntransitions)
		ncol=NstatesR+NtransitionsR+2
		//filename
        // open file
        fh = fopen(filename, "w")
        fput(fh,"{")
        fput(fh,`""Nstates":"' + Nstates + ",")
        fput(fh,`""Ntransitions":"' + Ntransitions + ",")
        fput(fh, `""xvalues": ["' + xvalues + "],")
        fput(fh, `""yvalues": ["' + yvalues + "],")
        fput(fh,`""ymin":"' + ymin + ",")
        fput(fh,`""ymax":"' + ymax + ",")
        fput(fh,`""xmin":"' + xmin + ",")
        fput(fh,`""xmax":"' + xmax + ",")
        fput(fh,`""boxwidth":"' + strofreal(boxwidth,"%9.5f") + ",")
        fput(fh,`""boxheight":"' + strofreal(boxheight,"%9.5f") + ",")
        fput(fh,`""statenames": ["' + statenames + `"],"')
		fput(fh,`""transnames": ["' + transnames + `"],"')
        fput(fh,`""arrows": "' + arrows + ",")
        fput(fh,`""arrowstext": "' + arrowstext+ ",")
		

		tmptmat  = `""tmat":["'  
        for (j=1; j<=NstatesR; j++) {  
        	tmptmat = tmptmat + "[" + invtokens(strofreal(tmat[j,]),",") + "]"
          if(j!=NstatesR) tmptmat = tmptmat +", "
        }
		tmptmat = subinstr(tmptmat,".", `"""'+"NA"+`"""')
		
		if (hasfreq==0) {	
				tmptmat  = tmptmat + "]"
				fput(fh,tmptmat)
				fput(fh,"}")
		}
	
	
		if (hasfreq==1) {	
				tmptmat  = tmptmat + "],"
	
	
				fput(fh,tmptmat)
	

  /********************************************/
	
				tmptmat2  =`""frequencies":["'
 
				for (i=1; i<=nrow-1; i++) {
				
					for (j=1; j<=1; j++) {
					
							tmptmat2 = tmptmat2 + "{" + `"""'+ clabels[j,2]+`"""'+`":"'+strofreal(frequencies[i,j])+`","'	
					}
					
					for (j=2; j<=ncol-1; j++) { 
					
							tmptmat2 = tmptmat2 +  `"""'+ clabels[j,2]+`"""'+`":"'+strofreal(frequencies[i,j])+`","'	
					}
					
					tmptmat2 = tmptmat2 +   `"""'+ clabels[j,2]+`"""'+`":"'+strofreal(frequencies[i,ncol],"%9.5f") 
					tmptmat2 = tmptmat2+ "},"
				}
     
				for (i=nrow; i<=nrow; i++) {
				
						for (j=1; j<=1; j++) {
								tmptmat2 = tmptmat2 + "{" + `"""'+ clabels[j,2]+`"""'+`":"'+strofreal(frequencies[i,j])+`","'	
						}
						for (j=2; j<=ncol-1; j++) { 
								tmptmat2 = tmptmat2 + `"""'+ clabels[j,2]+`"""'+`":"'+strofreal(frequencies[i,j])+`","'	
						}
						tmptmat2 = tmptmat2 +   `"""'+ clabels[j,2]+`"""'+`":"'+strofreal(frequencies[i,ncol],"%9.5f") 
						tmptmat2 = tmptmat2+ "}"
				}
 

				fput(fh,tmptmat2)

   
				fput(fh,`"]"')
		
			
  
				fput(fh,"}")
		
	  }

	  fclose(fh)

}
end
        
