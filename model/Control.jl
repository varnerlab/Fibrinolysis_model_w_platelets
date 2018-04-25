@everywhere function Control(t,x,PROBLEM_DICTIONARY)	
	num_rates = 24;
	idx = find(x->(x<0),x);
	x[idx] = 0.0;

	idx = find(x->(abs(x)<1E-15),x); #lets make anything less than 1E-9 zero
	#@show idx
	x[idx] = 0.0;
	
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
	#k_FX_activation = kinetic_parameter_vector[12]
	#K_FX_activation = kinetic_parameter_vector[13]
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
	k_cat_fibrinogen_deg=fibrin_kinetic_parameter_vector[21]
	Km_fibrinogen_deg=fibrin_kinetic_parameter_vector[22]

	# Initiation model -
	initiation_trigger_term = ((alpha_trigger_activation*TRIGGER)^order_trigger_activation)/(1 + ((alpha_trigger_activation*TRIGGER)^order_trigger_activation))
	initiation_TFPI_term = 1 - ((alpha_trigger_inhibition_TFPI*TFPI)^order_trigger_inhibition_TFPI)/(1 + ((alpha_trigger_inhibition_TFPI*TFPI)^order_trigger_inhibition_TFPI))

	# Amplification model -
	activation_term = ((alpha_amplification_FIIa*FIIa)^order_amplification_FIIa)/(1 + ((alpha_amplification_FIIa*FIIa)^order_amplification_FIIa))
	inhibition_term =  ((alpha_amplification_APC*APC)^order_amplification_APC)/(1 + ((alpha_amplification_APC*APC)^order_amplification_APC))
	inhibition_term_TFPI =  ((alpha_amplification_TFPI*TFPI)^order_amplification_TFPI)/(1 + ((alpha_amplification_TFPI*TFPI)^order_amplification_TFPI))
	#factor_product = FIX*FVIII*FV_FX/(nominal_FIX*nominal_FVIII*nominal_FV_X)
	#factor_amplification_term = ((0.1*factor_product)^2)/(1+((0.1*factor_product)^2))

	# Shutdown phase -
	shutdown_term = ((alpha_shutdown_APC*FIIa)^order_shutdown_APC)/(1 + ((alpha_shutdown_APC*FIIa)^order_shutdown_APC))


	#prothominabse complex formation
	activation_FV_by_thrombin = (alpha_FV_activation*FIIa)^order_FV_activation/(1+(alpha_FV_activation*FIIa)^order_FV_activation)
	activation_FX_by_trigger = (alpha_FX_activation*TRIGGER)^order_FX_activation/(1+(alpha_FX_activation*TRIGGER)^order_FX_activation)
	inhibition_of_FX_by_ATIII=1-(alpha_FX_inhibition*ATIII)^order_FX_inhibition/(1+(alpha_FX_inhibition*ATIII)^order_FX_inhibition)

	# control_term_fibrinolysis_TAFI -
	control_term_fibrinolysis_TAFI  = 1 - ((alpha_fib_inh_TAFI.*TAFI).^order_fib_inh_TAFI)/(1+((alpha_fib_inh_TAFI.*TAFI).^order_fib_inh_TAFI))
	# control_term_fibrinolysis_fXIII -
	control_term_fibrinolysis_fXIII  = 1 - ((alpha_fib_inh_fXIII.*FXIII).^order_fib_inh_fXIII)./(1+((alpha_fib_inh_fXIII.*FXIII).^order_fib_inh_fXIII))

	control_vector = ones(1,num_rates)
	control_vector[1] = min(initiation_trigger_term,initiation_TFPI_term)
	control_vector[2] = min(inhibition_term,inhibition_term_TFPI)
	control_vector[3] = shutdown_term
	control_vector[4] = 1
	control_vector[5] = 1
	#control_vector[6] = min(inhibition_term,inhibition_term_TFPI)
	control_vector[6] = FVIII_control
	control_vector[7] = min(1-max(inhibition_term,inhibition_term_TFPI))
	if (control_term_fibrinolysis_TAFI>0.001)
		control_vector[13] = control_term_fibrinolysis_TAFI
		control_vector[14] = control_term_fibrinolysis_TAFI
	end

	if (control_term_fibrinolysis_fXIII>0.001)
		control_vector[15] = control_term_fibrinolysis_fXIII
	end
	for j in collect(1:length(control_vector))
		if(isnan(control_vector[j])|| !isfinite(control_vector[j]))
		#	println("Found NAN at control_vector $j")
			control_vector[j] = 1.0
		end
	end
	return control_vector
end
