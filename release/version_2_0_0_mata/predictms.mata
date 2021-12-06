/*
predictms Mata sourcecode
*/

version 14.2

local RM 	real matrix
local ss 	string scalar
local RS	real scalar
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local SS	struct predictms_struct scalar
local PS	pointer(struct predictms_struct scalar) scalar
local PC	pointer(struct predictms_struct scalar) colvector
local Ps	pointer scalar

mata:

void predictms()
{
	`SS' S
	`RC' from
	`RS' Nstarts
	
	predictms_setup(S)
	
	from 	= strtoreal(tokens(st_local("from")))'
	Nstarts = rows(from)

	for (fr=1;fr<=Nstarts;fr++) {
		
		if (S.usesim) {
			for (at=1;at<=S.Nats;at++) predictms_core(S,from[fr],at)
		}
		else {
			for (at=1;at<=S.Nats;at++) predictms_aj(S,from[fr],at)
		}

		if (S.getprobs) predictms_post_predictions(S,from[fr],0)
	
		if (S.getlos) 	predictms_post_predictions(S,from[fr],1)

		if (S.hasuser) 	predictms_post_predictions(S,from[fr],2)
		
		if (S.getvisit) predictms_post_predictions(S,from[fr],3)
	}
	
}

`RM' ms_user_prob(`SS' S, `RS' b)
{
	return(S.pt[,b])
}

`RM' ms_user_los(`SS' S, `RS' b)
{
	return(S.los[,b])
}

end

