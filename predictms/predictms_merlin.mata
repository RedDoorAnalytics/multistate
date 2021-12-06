version 14.2

local gml 	struct merlin_struct scalar
local SS 	string scalar
local RS	real scalar
local RM	numeric matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct merlin_struct scalar) scalar
local SGS	struct predictms_struct scalar

mata:

`PS' predictms_merlin_setup(`SGS' S, `RR' b, `RS' N, `RS' trans, | `RC' t0)
{
	strtrans 	= strofreal(trans)
	at 			= S.at
	std 		= S.std
	
	//start by recalling merlin 
	
	stata("qui preserve")
	stata("local _tempobs = int("+strofreal(N)+")")
	stata("cap set obs "+st_local("_tempobs"))
	
	if (S.hasmodels) {
		stata("cap estimates restore "+st_local("modelests"+strtrans))
	}

	if (st_local("survsim")=="") {

		//start by replacing allvars with zeroes (if not in stdvars)
		allvars 	= tokens(st_global("e(allvars)"))
		Nallvars	= cols(allvars)
		zerovars 	= J(1,0,"")
		
		if (S.standardise) {
			stdvars = tokens(st_local("at"+strofreal(at)+"stdvars"+strtrans))
			for (i=1; i<=Nallvars; i++) {
				flag = 0
				for (j=1;j<=cols(stdvars);j++) {
					if (stdvars[j]==allvars[i]) flag = 1
				}
				if (!flag) zerovars = zerovars,allvars[i]
			}
		}
		else zerovars = allvars

		for (i=1;i<=cols(zerovars);i++) {
			stata("qui replace "+zerovars[i]+" = 0 if _n<="+strofreal(N))
		}

		//need to replace any variables with their at#()
		ats 		= tokens(st_local("at"+strofreal(at)))'
		Natstodo 	= rows(ats)/2
		i = 1
		j = 2
		for (a=1;a<=Natstodo;a++) {
			stata("qui replace "+ats[i]+" = "+ats[j]+" if _n<="+strofreal(N))
			i = i+2
			j = j+2
		}

		//standardisation
		if (S.standardise) {
			st_view(poststdvars=.,.,stdvars,st_local("stdtouse"))
			for (i=1;i<=cols(stdvars);i++) {
				stata("qui replace "+stdvars[i]+" = "+strofreal(poststdvars[std,i])+" if _n<="+strofreal(N))
			}
		}

		//if stacked model, update *_trans# variables
		if (!S.hasmodels) {

			for (i=1;i<=S.Ntrans;i++) {
				ti = "_trans"+strofreal(i)
				
				if (i==trans) 	stata("qui replace "+ti+" = 1 if _n<="+strofreal(N))
				else 			stata("qui replace "+ti+" = 0 if _n<="+strofreal(N))
				
				//rebuild any *_transi variables
				if (st_local("toupdate")!="") {
					for (k=1;k<=rows(S.toupdate);k++) {
						stata("cap replace "+S.toupdate[k]+ti+"= "+S.toupdate[k]+"*"+ti+" if _n<="+strofreal(N))
					}
				}
			}
			
		}

	}	

	//if I replace ltruncated() timevar with newly simulated entry time
	//--> everthing will work for multiple timescales (only time since entry + main one)
	//--> this must be after zeros and at, so it's replaced appropriately
	if (args()==5 & st_global("e(ltruncated1)")!="") {
		st_view(lt=.,.,st_global("e(ltruncated1)"))
		lt[|1,1\N,1|] = t0									//!! need to check any sorting doesn't screw this up
	}

	//post coefficients
	stata("tempname bmat")
	st_matrix(st_local("bmat"),b)

	//remove any options
	stata("local cmd "+st_global("e(cmdline)"))
	stata("gettoken merlin cmd : cmd")
	stata(`"gettoken cmd rhs : cmd, parse(",") bind"')

	//strip of any ifs/ins
		//global
		cmd = strtrim(st_local("cmd"))
		lastbracket = strrpos(cmd,")")
		cmd = substr(cmd,1,lastbracket)

		//model specific
		ifpos = strpos(cmd," if ")
		inpos = strpos(cmd," in ")
		comma = strpos(cmd,",")
		if (ifpos) {
			rest = substr(cmd,ifpos,.)
			comma = strpos(rest,",")
			cmd = substr(cmd,1,ifpos-1)+substr(rest,comma,.)
			
		}
		if (inpos) {
			rest = substr(cmd,inpos,.)
			comma = strpos(rest,",")
			cmd = substr(cmd,1,inpos-1)+substr(rest,comma,.)
		}
	
	//add global if
	cmd = cmd + " if _n<=" + st_local("_tempobs")
	st_local("cmd",cmd)

	//recall merlin
	stata("tempname tousem")
	stata("tempname GML"+strtrans)
	rmexternal(st_local("GML"+strtrans))

	cmd = ""
	cmd = cmd + " merlin_parse"
	cmd = cmd + " " + st_local("GML"+strtrans)
	cmd = cmd + " " + ", touse("+st_local("tousem")+") : " + st_local("cmd")
	cmd = cmd + " " + ", predict galahad nogen"
	cmd = cmd + " " + "from("+st_local("bmat")+")"
	cmd = cmd + " " + "npredict("+strofreal(N)+")"
	cmd = cmd + " " + "devcode1("+st_local("devcode1")+")"		//offsets i.e. multiple timescales
	cmd = cmd + " " + "devcode9("+st_local("devcode9")+")"		//Cox

	//run
	rc = _stata(cmd)
	if (rc) exit(1986)

	//get filled up struct
	struct merlin_struct gml
	swap(gml,*findexternal(st_local("GML"+strofreal(trans))))
	gml.fixedonly 	= st_local("marginal")==""

	stata("qui restore")
	return(&gml)
}

end
