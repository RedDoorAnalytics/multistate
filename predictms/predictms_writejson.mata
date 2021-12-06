version 12.1

mata:
function predictms_writejson() {
	filename = st_local("jsonfile")
	Nstates = strtoreal(st_local("Nstates"))
	Ntrans = strtoreal(st_local("Ntrans"))
	Nats = strtoreal(st_local("Nats"))
	touse = st_local("touse")
	timevar = st_data(.,st_local("timevar"),touse)'
	timevar = invtokens(strofreal(timevar),",")
	atlist = J(1,Nats,"")
	for(i=1;i<=Nats;i++) {
		atlist[1,i] = `"""' + st_local("at"+strofreal(i)) + `"""'
	}
	atlist = "[" + invtokens(atlist,",") + "]"

	
// open file
	fh = fopen(filename, "w")
	fput(fh,"mspred = {")
	for(j=1;j<=Nstates;j++) {
		Preturn = J(1,Nats,"")
		for(ati=1;ati<=Nats;ati++) {
			tmpP = st_data(.,"_prob_at"+strofreal(ati)+"_1_"+strofreal(j),touse)'
			Preturn[1,ati] = "[" + invtokens(strofreal(tmpP),",") + "]"
			Pname = "P"+strofreal(j)
		}
		Preturn = invtokens(Preturn,",")
		fput(fh,`"""' + Pname + `"": ["' + Preturn + "],")
	}
	fput(fh,`""timevar": ["' + timevar + "],")
	fput(fh,`""Nats": "' + strofreal(Nats) +",")
// 	for(i=1;i<=Ntrans;i++) {
// 		hreturn = J(1,Nats,"")
// 		for(j=1;j<=Nats;j++) {
// 			tmph = st_data(.,st_local("_hazard"+"_at"+strofreal(j)),touse)'
// 			hreturn[1,j] = "[" + invtokens(strofreal(tmph),",") + "]"
// 		}
// 		hname = "h"+strofreal(i)
// 		hreturn = invtokens(hreturn,",")
// 		hreturn = subinstr(hreturn,".,","null,")
// 		fput(fh,`"""' + hname + `"": ["' + hreturn + "],")
// 	}	
	fput(fh,`""atlist": "' + atlist)
	fput(fh,"}")
	fclose(fh)
}
end
	
