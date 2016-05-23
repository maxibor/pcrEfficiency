nb_cycle = 30
nb_amplicon_final = 160000
efficiency = 1.73


def nb_init_templates (nb_cycle, nb_amplicon_final, efficiency) :
	'''
	nb_cycle : number of cycles of the PCR
	nb_amplicon_final : number of amplicon at the end of PCR
	efficiency : calculated with PCR efficiency
	'''
	nb_init = (nb_amplicon_final) / (efficiency**nb_cycle)
	if nb_init < 1 :
		return 1
	else :
		return nb_init

print nb_init_templates(nb_cycle,nb_amplicon_final,efficiency) , "number of DNA templates at the beginning of PCR"
