function Kinetics(t,x,PROBLEM_DICTIONARY)
	num_rates = 24;
	idx = findall(x->(x<0),x);
	if(size(idx, 1)>0)
		x[idx] = zeros(size(idx))
	end
	
	# Alias the species -
	FII	 = x[1]
	FIIa	= x[2]
	PC	  = x[3]
	APC	 = x[4]
	ATIII   = x[5]
	TM	  = x[6]
	TRIGGER = x[7]
	Eps = x[8] #frac platelets acativated
	FV_FX = x[9]
	FV_FXA = x[10]
	PROTHOMBINASE_PLATELETS = x[11]

	Fibrin = x[12]
	Plasmin = x[13]
	Fibrinogen = x[14]
	Plasminogen = x[15]
	tPA = x[16]
	uPA = x[17]
	Fibrin_monomer = x[18]
	Protofibril = x[19]
	antiplasmin =x[20]
	PAI_1 = x[21]
	Fiber =x[22] 

	# Grab the kinetic parameetrs from the problem dictionary -
	kinetic_parameter_vector = PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"]
	control_parameter_vector = PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"]
	qualitative_factor_level_vector = PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"]
	platelet_parameter_vector = PROBLEM_DICTIONARY["PLATELET_PARAMS"]
	timing = PROBLEM_DICTIONARY["TIME_DELAY"]
	fibrin_kinetic_parameter_vector = PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]
	fibrin_control_vector = PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]

	# Alias the qualitative factors -
	TFPI = qualitative_factor_level_vector[1]
	#FV = qualitative_factor_level_vector[2]
	FVIII = qualitative_factor_level_vector[3]
	FIX = qualitative_factor_level_vector[4]
	#FX = qualitative_factor_level_vector[5]
	Platelets = qualitative_factor_level_vector[6]
	TAFI=qualitative_factor_level_vector[7]
	FXIII = qualitative_factor_level_vector[8]

	alpha_trigger_activation = control_parameter_vector[1]
	order_trigger_activation = control_parameter_vector[2]
	alpha_trigger_inhibition_APC = control_parameter_vector[3]
	order_trigger_inhibition_APC = control_parameter_vector[4]
	alpha_trigger_inhibition_TFPI = control_parameter_vector[5]
	order_trigger_inhibition_TFPI = control_parameter_vector[6]
	
	# Amplification -
	alpha_amplification_FIIa = control_parameter_vector[7]
	order_amplification_FIIa = control_parameter_vector[8]
	alpha_amplification_APC = control_parameter_vector[9]
	order_amplification_APC = control_parameter_vector[10]
	alpha_amplification_TFPI = control_parameter_vector[11]
	order_amplification_TFPI = control_parameter_vector[12]
	
	# APC generation -
	alpha_shutdown_APC = control_parameter_vector[13]
	order_shutdown_APC = control_parameter_vector[14]

	#prothomibase complex
	alpha_FV_activation=control_parameter_vector[15]
	order_FV_activation=control_parameter_vector[16]
 	alpha_FX_activation=control_parameter_vector[17]
    	order_FX_activation=control_parameter_vector[18]
        alpha_FX_inhibition=control_parameter_vector[19]
        order_FX_inhibition=control_parameter_vector[20]
	
	#platelets
	kplatelts = platelet_parameter_vector[1]#1 rate constant
	platelet_pwr = platelet_parameter_vector[2] #2 power for control function
	platelet_denom = platelet_parameter_vector[3] #3 adjustment in denominator
	EpsMax0 = platelet_parameter_vector[4] #4 Epsmax0
	aida = platelet_parameter_vector[5] #5 aida
	koffplatelets = platelet_parameter_vector[6]

	#fibrin controls
	alpha_fib_inh_fXIII=fibrin_control_vector[1]
	order_fib_inh_fXIII=fibrin_control_vector[2]
	alpha_fib_inh_TAFI=fibrin_control_vector[3]
	order_fib_inh_TAFI=fibrin_control_vector[4]
	alpha_tPA_inh_PAI_1=fibrin_control_vector[5]
	order_tPA_inh_PAI_1=fibrin_control_vector[6]
	alpha_uPA_inh_PAI_1=fibrin_control_vector[7]
	order_uPA_inh_PAI_1=fibrin_control_vector[8]

	#FVIII function
	FVIII_control = PROBLEM_DICTIONARY["FVIII_CONTROL"]

	#platelet control
	#update aleph so that it holds the maximum value of FIIa
	if(FIIa>PROBLEM_DICTIONARY["ALEPH"])
		PROBLEM_DICTIONARY["ALEPH"]=FIIa
	end

	#alias kinetic parameters
	k_trigger = kinetic_parameter_vector[1]
	K_trigger = kinetic_parameter_vector[2]
	k_amplification = kinetic_parameter_vector[3]
	K_FII_amplification = kinetic_parameter_vector[4]
	k_APC_formation = kinetic_parameter_vector[5]
	K_PC_formation = kinetic_parameter_vector[6]
	k_inhibition = kinetic_parameter_vector[7]
	K_FIIa_inhibition = kinetic_parameter_vector[8]
	k_inhibition_ATIII = kinetic_parameter_vector[9]
	k_FV_X_activation = kinetic_parameter_vector[10]
	K_FV_X_actiation = kinetic_parameter_vector[11]
	k_complex = kinetic_parameter_vector[14]
	k_amp_prothombinase = kinetic_parameter_vector[15]
	K_FII_amp_prothombinase=kinetic_parameter_vector[16]
	k_amp_active_factors= kinetic_parameter_vector[17]
	K_amp_active_factors = kinetic_parameter_vector[18]

	k_cat_Fibrinogen = fibrin_kinetic_parameter_vector[1]
	Km_Fibrinogen=fibrin_kinetic_parameter_vector[2]
	k_fibrin_monomer_association=fibrin_kinetic_parameter_vector[3]
	k_protofibril_association=fibrin_kinetic_parameter_vector[4]
	k_protofibril_monomer_association=fibrin_kinetic_parameter_vector[5]
	k_cat_plasminogen_Fibrin_tPA=fibrin_kinetic_parameter_vector[6]
	Km_plasminogen_Fibrin_tPA=fibrin_kinetic_parameter_vector[7]
	k_cat_plasminogen_Fibrin_uPA=fibrin_kinetic_parameter_vector[8]
	Km_plasminogen_Fibrin_uPA=fibrin_kinetic_parameter_vector[9]
	k_cat_Fibrin=fibrin_kinetic_parameter_vector[10]
	Km_Fibrin=fibrin_kinetic_parameter_vector[11]
	k_inhibit_plasmin=fibrin_kinetic_parameter_vector[12]
	k_inhibit_PAI_1_APC=fibrin_kinetic_parameter_vector[13]
	k_inhibit_tPA_PAI_1=fibrin_kinetic_parameter_vector[14]
	k_inhibit_uPA_PAI_1=fibrin_kinetic_parameter_vector[15]
	k_fibrin_formation=fibrin_kinetic_parameter_vector[16]
	k_cat_fiber=fibrin_kinetic_parameter_vector[17]
	Km_fiber=fibrin_kinetic_parameter_vector[18]
	Ka_plasminogen_Fibrin_tPA=fibrin_kinetic_parameter_vector[19]
	K_plasminogen_Fibrin_tPA = fibrin_kinetic_parameter_vector[20]
	k_cat_fibrin_monomer_deg=fibrin_kinetic_parameter_vector[21]
	Km_fibrin_monomer_deg=fibrin_kinetic_parameter_vector[22]

	aleph = PROBLEM_DICTIONARY["ALEPH"]
	faleph = aleph^platelet_pwr/(aleph^platelet_pwr + platelet_denom^platelet_pwr)
	EpsMax = EpsMax0+(1-EpsMax0)*faleph
	#@show aleph, falpeh, EpsMax

	#timing
	time_delay = timing[1]
	time_coeff = timing[2]

    	rate_vector = similar(x, num_rates)#zeros(1,num_rates)
	rate_vector[1] = k_trigger*TRIGGER*(FV_FX/(K_trigger + FV_FX))
	rate_vector[2] = k_amplification*FIIa*(FII/(K_FII_amplification + FII))
	rate_vector[3] = k_APC_formation*TM*((PC)/(K_PC_formation + PC))
	rate_vector[4] = (k_inhibition_ATIII)*(ATIII)*((FIIa)^1.26)
	rate_vector[5] = k_complex*FV_FXA*aida/Eps
	rate_vector[6] = k_amp_prothombinase*PROTHOMBINASE_PLATELETS*FII/(K_FII_amp_prothombinase + FII)
	rate_vector[7] = k_amp_active_factors*FV_FXA*FII/(K_amp_active_factors+FII)
	rate_vector[8] = (FIIa*k_cat_Fibrinogen*Fibrinogen)/(Km_Fibrinogen+Fibrinogen)                             # Cleavage of fibrinopeptides A and/or B to form fibrin monomer
	rate_vector[9] = k_fibrin_monomer_association*(Fibrin_monomer^2);                                          # Protofibril formation through association of fibrin monomers
	rate_vector[10] = k_protofibril_association*(Protofibril^2);                                                # Association of protofibril-protofibril to form fibers
	rate_vector[11] = k_protofibril_monomer_association*(Protofibril)*(Fibrin_monomer);                         # Fibril association with monomer forms protofibril again
	rate_vector[12] = k_fibrin_formation*(Fiber);                                                               # Fibrin growth
  	rate_vector[13] = (tPA*k_cat_plasminogen_Fibrin_tPA*Plasminogen)/(Km_plasminogen_Fibrin_tPA+Plasminogen);
  	rate_vector[14] = (uPA*k_cat_plasminogen_Fibrin_uPA*Plasminogen)/(Km_plasminogen_Fibrin_uPA+Plasminogen);
  	rate_vector[15] = (Plasmin*k_cat_Fibrin*Fibrin)/(Km_Fibrin+Fibrin)
  	rate_vector[16] = k_inhibit_plasmin*Plasmin*(antiplasmin);
  	rate_vector[17] = k_inhibit_PAI_1_APC*(APC)*(PAI_1);
  	rate_vector[18] = k_inhibit_tPA_PAI_1*(tPA)*PAI_1;
 	rate_vector[19] = k_inhibit_uPA_PAI_1*(uPA)*PAI_1;
  	rate_vector[20] = (Plasmin*k_cat_fiber*Fiber)/(Km_fiber+Fiber)                                           # Plasmin degrading fiber
  	rate_vector[21] = (Plasmin*k_cat_fibrin_monomer_deg*Fibrin_monomer)/(Km_fibrin_monomer_deg+Fibrin_monomer)	#plasmin degrading fibrin monomer
	rate_vector[22]= kplatelts*(EpsMax-Eps) #platelet activiation
	rate_vector[23]= koffplatelets*Eps  #platelet degredation
	#rate_vector[24]=kincorportation*Eps*Fiber

	for j in collect(1:length(rate_vector))
		if(isnan(rate_vector[j]) || !isfinite(rate_vector[j]))
		#	println("Found NAN at rate_vector $j")
			rate_vector[j] = 0.0
		end
	end

	return rate_vector
end
