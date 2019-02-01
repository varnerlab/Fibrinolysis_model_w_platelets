 function buildCoagulationModelDictionary(platelet_count)
    
    # Initialize -
    PROBLEM_DICTIONARY = Dict()
    # Initial condition -
	#from simulateFig3Butenas2002
    initial_condition_vector = Float64[]
    push!(initial_condition_vector,1400)       # 1 FII
    push!(initial_condition_vector,0.001)        # 2 FIIa
    push!(initial_condition_vector,60.0)        # 3 PC
    push!(initial_condition_vector,0.0)        # 4 APC
    push!(initial_condition_vector,3400)       # 5 ATIII
    push!(initial_condition_vector,12)         # 6 TM
    push!(initial_condition_vector,5.0E-3)          # 7 TRIGGER
    push!(initial_condition_vector, 0.01) #8 Fraction of platelets activated
    push!(initial_condition_vector,22.0)          #9  FV+FX
    push!(initial_condition_vector,0.0)          #10  FVa+FXa
    push!(initial_condition_vector,0.000)         #11  prothombinase complex
    push!(initial_condition_vector, 0.0)	#12 Fibrin
    push!(initial_condition_vector,0.0)		#13 plasmin
    push!(initial_condition_vector,7400)	#14 Fibrinogen
    push!(initial_condition_vector, 2000)	#15 plasminogen
    push!(initial_condition_vector,0.0)		#16 tPA
    push!(initial_condition_vector,0.0)		#17 uPA
    push!(initial_condition_vector,0.0)		#18 fibrin monomer
    push!(initial_condition_vector,0.0)		#19 protofibril
    push!(initial_condition_vector,1180)	#20 antiplasmin
    push!(initial_condition_vector,.56)		#21 PAI_1
    push!(initial_condition_vector,0.0)		#22 Fiber
    PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"] = initial_condition_vector
    
    # Kinetic parameters -
    kinetic_parameter_vector = Float64[]
    push!(kinetic_parameter_vector,7200/1.2)       #  k_trigger
    push!(kinetic_parameter_vector,100)       # 1 K_trigger
    push!(kinetic_parameter_vector,5.0)        # 2 k_amplification
    push!(kinetic_parameter_vector,12)       # 3 K_FII_amplification
    push!(kinetic_parameter_vector,.5)        # 4 k_APC_formation
    push!(kinetic_parameter_vector,3/100.0)         # 5 K_PC_formation
    push!(kinetic_parameter_vector,0.2*10)        # 6 k_inhibition
    push!(kinetic_parameter_vector,120)       # 7 K_FIIa_inhibition
    push!(kinetic_parameter_vector,3E-4)     # 8 k_inhibition_ATIII
    #push!(kinetic_parameter_vector,0.001)      # 9 K_inhibition_ATIII
    #push!(kinetic_parameter_vector,100.0)      # 10 K_inhibition_FIIa
    push!(kinetic_parameter_vector, 2E7*60*10.0^-4) #9 k_FV_activation, from reaction 16 in Diamond 2010 paper
    push!(kinetic_parameter_vector, 1E8*80*10.0^-6/100) #10 K_FV_activation 
    push!(kinetic_parameter_vector, 600.0) #11 k_FX_activation from reaction 6 in Diamond 2010 paper
    push!(kinetic_parameter_vector, .28) #12 K_FX_activation
    push!(kinetic_parameter_vector, 24) #13 k_complex
    push!(kinetic_parameter_vector, 63.5*60*2.5 )#14 k_amp_prothombinase from reaction 18 in Diamond 2010
    push!(kinetic_parameter_vector,10) #13 K_FII_amp_prothombinase
    push!(kinetic_parameter_vector, 1.0) #k_amp_active_factors
    push!(kinetic_parameter_vector, 1.0) #K_amp_active_factor


    PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"] = kinetic_parameter_vector
    
    # Control parameters -
    control_parameter_vector =Float64[]
    # Trigger -
    push!(control_parameter_vector,1.0)      # 0 9 alpha_trigger_activation = control_parameter_vector[0]
    push!(control_parameter_vector,1.0)        # 1 10 order_trigger_activation = control_parameter_vector[1]
    push!(control_parameter_vector,1.0)        # 2 11 alpha_trigger_inhibition_APC = control_parameter_vector[2]
    push!(control_parameter_vector,1.0)        # 3 12 order_trigger_inhibition_APC = control_parameter_vector[3]
    push!(control_parameter_vector,0.1)        # 4 13 alpha_trigger_inhibition_TFPI = control_parameter_vector[4]
    push!(control_parameter_vector,1.0)        # 5 14 order_trigger_inhibition_TFPI = control_parameter_vector[5]
    
    # Amplification -
    push!(control_parameter_vector,0.1)        # 6 15 alpha_amplification_FIIa = control_parameter_vector[6]
    push!(control_parameter_vector,2.0)        # 7 16 order_amplification_FIIa = control_parameter_vector[7]
    push!(control_parameter_vector,0.4)        # 8 17 alpha_amplification_APC = control_parameter_vector[8]
    push!(control_parameter_vector,2.0)        # 9 18 order_amplification_APC = control_parameter_vector[9]
    push!(control_parameter_vector,0.01*100)       # 10 19 alpha_amplification_TFPI = control_parameter_vector[10]
    push!(control_parameter_vector,2.0)        # 11 20 order_amplification_TFPI = control_parameter_vector[11]
    
    # APC generation -
    push!(control_parameter_vector,2.0)        # 12 21 alpha_shutdown_APC = control_parameter_vector[12]
    push!(control_parameter_vector,2.0)        # 13 22 order_shutdown_APC = control_parameter_vector[13]

    push!(control_parameter_vector,1.0) #14 alpha_FV_activation
    push!(control_parameter_vector,1.0) #15 order_FV_activation
    push!(control_parameter_vector,1.0) #16 alpha_FX_activation
    push!(control_parameter_vector,1.0) #17 order_FX_activation
    push!(control_parameter_vector,1.0) #14 alpha_FX_inhibition
    push!(control_parameter_vector,1.0) #15 order_FX_inhibition
    PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"] = control_parameter_vector

   #platlet controls
	platelet_parameter_vector = Float64[]
	normal_platelet_count = 300 #*10^6 #/mL
	push!(platelet_parameter_vector, .005*5) #1 rate constant
	push!(platelet_parameter_vector, 1.6123) #2 power for control function
	push!(platelet_parameter_vector, 2.4279E-9) #3 adjustment in denominator
	push!(platelet_parameter_vector, .01) #4 Epsmax0
	push!(platelet_parameter_vector, 1.05*platelet_count/normal_platelet_count)  #5 aida
	push!(platelet_parameter_vector, .005) #koffplatelets
	PROBLEM_DICTIONARY["PLATELET_PARAMS"] = platelet_parameter_vector
   

   PROBLEM_DICTIONARY["ALEPH"] = initial_condition_vector[2]

	timing = Float64[]
	push!(timing, 0.0) #time_delay
	push!(timing, 1.0) #coeff
	PROBLEM_DICTIONARY["TIME_DELAY"]=timing
    
    # QFactor vector - from simulateFig3Butenas2002
    qualitative_factor_vector =Float64[]
   push!(qualitative_factor_vector,2.5)           # 1 TFPI
   push!(qualitative_factor_vector,20.0)          # 2 FV
   push!(qualitative_factor_vector,0.7)           # 3 FVIII
   push!(qualitative_factor_vector,90.0)           # 4 FIX
   push!(qualitative_factor_vector,170.0)         # 5 FX
   push!(qualitative_factor_vector,1.0)           # 6 Platelets
   push!(qualitative_factor_vector, 70)		#7 TAFI
   push!(qualitative_factor_vector,93)		#8 FXIII
    PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"] = qualitative_factor_vector

	fibrin_kinetic_parameters = Float64[]
	  # Kinetic parameters for fibrin generation
	  push!(fibrin_kinetic_parameters,0.01)            # 1 k_cat_Fibrinogen
	  push!(fibrin_kinetic_parameters,1000.0)         # 2 Km_Fibrinogen
	  push!(fibrin_kinetic_parameters,0.01)         # 3 k_fibrin_monomer_association
	  push!(fibrin_kinetic_parameters,0.01)         # 4 k_protofibril_association
	  push!(fibrin_kinetic_parameters,0.01)         # 5 k_protofibril_monomer_association

	  # Kinetic parameters for fibrinolysis
	  push!(fibrin_kinetic_parameters,2.0)           # 6 k_cat_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,10.0)          # 7 Km_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,0.1)           # 8 k_cat_plasminogen_Fibrin_uPA
	  push!(fibrin_kinetic_parameters,10.0)          # 9 Km_plasminogen_Fibrin_uPA
	  push!(fibrin_kinetic_parameters,1.0)           # 10 k_cat_Fibrin
	  push!(fibrin_kinetic_parameters,10.0)          # 11 Km_Fibrin
	  push!(fibrin_kinetic_parameters,0.01)          # 12 k_inhibit_plasmin
	  push!(fibrin_kinetic_parameters,0.01)          # 13 k_inhibit_PAI_1_APC
	  push!(fibrin_kinetic_parameters,0.01)          # 14 k_inhibit_tPA_PAI_1
	  push!(fibrin_kinetic_parameters,0.01)          # 15 k_inhibit_uPA_PAI_1

	  # Kinetic parameter for fibrin formation
	  push!(fibrin_kinetic_parameters,0.01)          # 16 k_fibrin_formation

	  # Fiber degradation
	  push!(fibrin_kinetic_parameters,100.0)          # 17 k_cat_fiber
	  push!(fibrin_kinetic_parameters,10.0)          # 18 Km_fiber

	  # Plasmin activation
	  push!(fibrin_kinetic_parameters,0.01)          # 19 Km_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,0.01)          # 20 K_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,1.00)          # 21 k_cat_fibrinogen_deg
	  push!(fibrin_kinetic_parameters,10.0)          # 22 Km_fibrinogen_deg
	PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]=fibrin_kinetic_parameters

	fibrin_control_parameters = Float64[]
	# Control constants for TAFI inhibition -
	  push!(fibrin_control_parameters,0.1)          # 1 alpha_fib_inh_fXIII
	  push!(fibrin_control_parameters,2.0)          # 2 order_fib_inh_fXIII
	  push!(fibrin_control_parameters,0.1)          # 3 alpha_fib_inh_TAFI
	  push!(fibrin_control_parameters,2.0)          # 4 order_fib_inh_TAFI

	  # Control constants for plasmin activation term
	  push!(fibrin_control_parameters,0.1)          # 5 alpha_tPA_inh_PAI_1
	  push!(fibrin_control_parameters,2.0)          # 6 order_tPA_inh_PAI_1
	  push!(fibrin_control_parameters,0.1)          # 7 alpha_uPA_inh_PAI_1
	  push!(fibrin_control_parameters,2.0)          # 8 order_uPA_inh_PAI_1

	PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]=fibrin_control_parameters

	PROBLEM_DICTIONARY["FVIII_CONTROL"] =  5.20895
	PROBLEM_DICTIONARY["parameter_name_mapping_array"]=["k_trigger","K_trigger" ,"k_amplification" ,"K_FII_amplification" ,"k_APC_formation" ,"K_PC_formation" ,"k_inhibition" ,"K_FIIa_inhibition" ,"k_inhibition_ATIII" ,"k_FV_X_activation", "K_FV_X_actiation" ,"k_FX_activation","k_FX_activation","k_complex","k_amp_prothombinase","K_FII_amp_prothombinase","k_amp_active_factors", "k_amp_active_factors","alpha_trigger_activation","order_trigger_activation","alpha_trigger_inhibition_APC","order_trigger_inhibition_APC","alpha_trigger_inhibition_TFPI", "order_trigger_inhibition_TFPI","alpha_amplification_FIIa","order_amplification_FIIa","alpha_amplification_APC","order_amplification_APC","alpha_amplification_TFPI","order_amplification_TFPI","alpha_shutdown_APC","order_shutdown_APC","alpha_FV_activation","order_FV_activation","alpha_FX_activation","order_FX_activation","alpha_FX_inhibition","order_FX_inhibition","platelet_rate_constant",
"platelet_power","platelet_denominator","Epsmax0","aida","koffplatelets","time_delay","coeff", "platelets"
,"k_cat_Fibrinogen", "Km_Fibrinogen", "k_fibrin_monomer_association", "k_protofibril_association","k_protofibril_monomer_association", "k_cat_plasminogen_Fibrin_tPA", "Km_plasminogen_Fibrin_tPA",
"k_cat_plasminogen_Fibrin_uPA", "Km_plasminogen_Fibrin_uPA", "k_cat_Fibrin", "Km_Fibrin", "k_inhibit_plasmin","k_inhibit_PAI_1_APC", "k_inhibit_tPA_PAI_1", "k_inhibit_uPA_PAI_1",
"k_fibrin_formation", "k_cat_fiber", "Km_fiber", "Ka_plasminogen_Fibrin_tPA","K_plasminogen_Fibrin_tPA","k_cat_fibrinogen_deg", "Km_fibrinogen_deg",
"alpha_fib_inh_fXIII", "order_fib_inh_fXIII", "alpha_fib_inh_TAFI", "order_fib_inh_TAFI", "order_tPA_inh_PAI_1", "alpha_tPA_inh_PAI_1", "alpha_uPA_inh_PAI_1",
"order_uPA_inh_PAI_1"] 
	PROBLEM_DICTIONARY["number_of_states"]=22
    
    return PROBLEM_DICTIONARY
end

 function buildCoagulationModelDictionary(kinetic_parameter_vector, control_parameter_vector, platelet_parameter_vector,timing, platelet_count)
    
    # Initialize -
    PROBLEM_DICTIONARY = Dict()
    normal_platelet_count = 300 #*10^6 #/mL
    # Initial condition -
	#from simulateFig3Butenas2002
    initial_condition_vector = Float64[]
    push!(initial_condition_vector,1400)       # 0 FII
    push!(initial_condition_vector,0.001)        # 1 FIIa
    push!(initial_condition_vector,60.0)        # 2 PC
    push!(initial_condition_vector,0.0)        # 3 APC
    push!(initial_condition_vector,3400)       # 4 ATIII
    push!(initial_condition_vector,12)         # 5 TM
    push!(initial_condition_vector,5.0E-3)          # 6 TRIGGER
    push!(initial_condition_vector, 0.01) #Fraction of platelets activated
    push!(initial_condition_vector,22.0)          #  FV+FX
    push!(initial_condition_vector,0.0)          #  FVa+FXa
    push!(initial_condition_vector,0.000)         #  prothombinase complex
    push!(initial_condition_vector, 0.0)	#Fibrin
    push!(initial_condition_vector,0.0)		#plasmin
    push!(initial_condition_vector,7400)	#Fibrinogen
    push!(initial_condition_vector, 2000)	#plasminogen
    push!(initial_condition_vector,0.0)		#tPA
    push!(initial_condition_vector,0.0)		#uPA
    push!(initial_condition_vector,0.0)		#fibrin monomer
    push!(initial_condition_vector,0.0)		#protofibril
    push!(initial_condition_vector,1180)	#antiplasmin
    push!(initial_condition_vector,.56)		#PAI_1
    push!(initial_condition_vector,0.0)		#Fiber
    PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"] = initial_condition_vector

    PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"] = kinetic_parameter_vector
    
    PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"] = control_parameter_vector

   #platlet controls
	PROBLEM_DICTIONARY["PLATELET_PARAMS"] = platelet_parameter_vector
	PROBLEM_DICTIONARY["PLATELET_PARAMS"][5] = PROBLEM_DICTIONARY["PLATELET_PARAMS"][5]*platelet_count; #correct for differing amounts of platelets
    
    # Experimental output scaling -
    # Fig 5A scaling (3.20,1.0)
    # Fig 3 scaling (-1.0,1.0)
    # Fig 2 scaling (1.8)
    # Fig 1 allen scaling (4.0,1.0)
    scaling_parameter_vector =Float64[]
    push!(scaling_parameter_vector,4.0)        # 0 Time scale 
    push!(scaling_parameter_vector,1.0)        # 1 Abundance scale
    PROBLEM_DICTIONARY["SCALING_PARAMETER_VECTOR"] = scaling_parameter_vector

   PROBLEM_DICTIONARY["ALEPH"] = initial_condition_vector[2]+.001
	PROBLEM_DICTIONARY["TIME_DELAY"]=timing
     # QFactor vector - from simulateFig3Butenas2002
    qualitative_factor_vector =Float64[]
   push!(qualitative_factor_vector,2.5)           # 1 TFPI
   push!(qualitative_factor_vector,20.0)          # 2 FV
   push!(qualitative_factor_vector,0.7)           # 3 FVIII
   push!(qualitative_factor_vector,90.0)           # 4 FIX
   push!(qualitative_factor_vector,170.0)         # 5 FX
   push!(qualitative_factor_vector,1.0)           # 6 Platelets
   push!(qualitative_factor_vector, 70)		#7 TAFI
   push!(qualitative_factor_vector,93)		#8 FXIII
    PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"] = qualitative_factor_vector

	fibrin_kinetic_parameters = Float64[]
	  # Kinetic parameters for fibrin generation
	  push!(fibrin_kinetic_parameters,0.01)            # 1 k_cat_Fibrinogen
	  push!(fibrin_kinetic_parameters,1000.0)         # 2 Km_Fibrinogen
	  push!(fibrin_kinetic_parameters,0.01)         # 3 k_fibrin_monomer_association
	  push!(fibrin_kinetic_parameters,0.01)         # 4 k_protofibril_association
	  push!(fibrin_kinetic_parameters,0.01)         # 5 k_protofibril_monomer_association

	  # Kinetic parameters for fibrinolysis
	  push!(fibrin_kinetic_parameters,2.0)           # 6 k_cat_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,10.0)          # 7 Km_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,0.1)           # 8 k_cat_plasminogen_Fibrin_uPA
	  push!(fibrin_kinetic_parameters,10.0)          # 9 Km_plasminogen_Fibrin_uPA
	  push!(fibrin_kinetic_parameters,1.0)           # 10 k_cat_Fibrin
	  push!(fibrin_kinetic_parameters,10.0)          # 11 Km_Fibrin
	  push!(fibrin_kinetic_parameters,0.01)          # 12 k_inhibit_plasmin
	  push!(fibrin_kinetic_parameters,0.01)          # 13 k_inhibit_PAI_1_APC
	  push!(fibrin_kinetic_parameters,0.01)          # 14 k_inhibit_tPA_PAI_1
	  push!(fibrin_kinetic_parameters,0.01)          # 15 k_inhibit_uPA_PAI_1

	  # Kinetic parameter for fibrin formation
	  push!(fibrin_kinetic_parameters,0.01)          # 16 k_fibrin_formation

	  # Fiber degradation
	  push!(fibrin_kinetic_parameters,100.0)          # 17 k_cat_fiber
	  push!(fibrin_kinetic_parameters,10.0)          # 18 Km_fiber

	  # Plasmin activation
	  push!(fibrin_kinetic_parameters,0.01)          # 19 Km_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,0.01)          # 20 K_plasminogen_Fibrin_tPA
	  push!(fibrin_kinetic_parameters,1.00)          # 21 k_cat_fibrinogen_deg
	  push!(fibrin_kinetic_parameters,10.0)          # 22 Km_fibrinogen_deg
	PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]=fibrin_kinetic_parameters

	fibrin_control_parameters = Float64[]
	# Control constants for TAFI inhibition -
	  push!(fibrin_control_parameters,0.1)          # 1 alpha_fib_inh_fXIII
	  push!(fibrin_control_parameters,2.0)          # 2 order_fib_inh_fXIII
	  push!(fibrin_control_parameters,0.1)          # 3 alpha_fib_inh_TAFI
	  push!(fibrin_control_parameters,2.0)          # 4 order_fib_inh_TAFI

	  # Control constants for plasmin activation term
	  push!(fibrin_control_parameters,0.1)          # 5 alpha_tPA_inh_PAI_1
	  push!(fibrin_control_parameters,2.0)          # 6 order_tPA_inh_PAI_1
	  push!(fibrin_control_parameters,0.1)          # 7 alpha_uPA_inh_PAI_1
	  push!(fibrin_control_parameters,2.0)          # 8 order_uPA_inh_PAI_1

	PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]=fibrin_control_parameters

	PROBLEM_DICTIONARY["FVIII_CONTROL"] = 5.20895
	PROBLEM_DICTIONARY["parameter_name_mapping_array"]=["k_trigger","K_trigger" ,"k_amplification" ,"K_FII_amplification" ,"k_APC_formation" ,"K_PC_formation" ,"k_inhibition" ,"K_FIIa_inhibition" ,"k_inhibition_ATIII" ,"k_FV_X_activation", "K_FV_X_actiation" ,"k_FX_activation","k_FX_activation","k_complex","k_amp_prothombinase","K_FII_amp_prothombinase","k_amp_active_factors", "k_amp_active_factors","alpha_trigger_activation","order_trigger_activation","alpha_trigger_inhibition_APC","order_trigger_inhibition_APC","alpha_trigger_inhibition_TFPI", "order_trigger_inhibition_TFPI","alpha_amplification_FIIa","order_amplification_FIIa","alpha_amplification_APC","order_amplification_APC","alpha_amplification_TFPI","order_amplification_TFPI","alpha_shutdown_APC","order_shutdown_APC","alpha_FV_activation","order_FV_activation","alpha_FX_activation","order_FX_activation","alpha_FX_inhibition","order_FX_inhibition","platelet_rate_constant",
"platelet_power","platelet_denominator","Epsmax0","aida","koffplatelets","time_delay","coeff", "platelets"
,"k_cat_Fibrinogen", "Km_Fibrinogen", "k_fibrin_monomer_association", "k_protofibril_association","k_protofibril_monomer_association", "k_cat_plasminogen_Fibrin_tPA", "Km_plasminogen_Fibrin_tPA",
"k_cat_plasminogen_Fibrin_uPA", "Km_plasminogen_Fibrin_uPA", "k_cat_Fibrin", "Km_Fibrin", "k_inhibit_plasmin","k_inhibit_PAI_1_APC", "k_inhibit_tPA_PAI_1", "k_inhibit_uPA_PAI_1",
"k_fibrin_formation", "k_cat_fiber", "Km_fiber", "Ka_plasminogen_Fibrin_tPA","K_plasminogen_Fibrin_tPA","k_cat_fibrinogen_deg", "Km_fibrinogen_deg",
"alpha_fib_inh_fXIII", "order_fib_inh_fXIII", "alpha_fib_inh_TAFI", "order_fib_inh_TAFI", "order_tPA_inh_PAI_1", "alpha_tPA_inh_PAI_1", "alpha_uPA_inh_PAI_1",
"order_uPA_inh_PAI_1"] 

	PROBLEM_DICTIONARY["number_of_states"]=22
    return PROBLEM_DICTIONARY
end

 function buildCoagulationModelDictionary(kinetic_parameter_vector, control_parameter_vector, platelet_parameter_vector,timing, platelet_count,fibrin_kinetic_parameters,fibrin_control_parameters)
    
    # Initialize -
    PROBLEM_DICTIONARY = Dict()
    normal_platelet_count = 300 #*10^6 #/mL
    # Initial condition -
	#from simulateFig3Butenas2002
    initial_condition_vector = Float64[]
    push!(initial_condition_vector,1400)       # 1 FII
    push!(initial_condition_vector,0.001)        # 2 FIIa
    push!(initial_condition_vector,60.0)        # 3 PC
    push!(initial_condition_vector,0.0)        # 4 APC
    push!(initial_condition_vector,3400)       # 5 ATIII
    push!(initial_condition_vector,12)         # 6 TM
    push!(initial_condition_vector,5.0E-3)          # 7 TRIGGER
    push!(initial_condition_vector, 0.01) #8 Fraction of platelets activated
    push!(initial_condition_vector,22.0)          #9  FV+FX
    push!(initial_condition_vector,0.0)          #10  FVa+FXa
    push!(initial_condition_vector,0.000)         #11  prothombinase complex
    push!(initial_condition_vector, 0.0)	#12 Fibrin
    push!(initial_condition_vector,0.0)		#13 plasmin
    push!(initial_condition_vector,7400)	#14 Fibrinogen
    push!(initial_condition_vector, 2000)	#15 plasminogen
    push!(initial_condition_vector,0.0)		#16 tPA
    push!(initial_condition_vector,0.0)		#uPA
    push!(initial_condition_vector,0.0)		#fibrin monomer
    push!(initial_condition_vector,0.0)		#protofibril
    push!(initial_condition_vector,1180)	#antiplasmin
    push!(initial_condition_vector,.56)		#PAI_1
    push!(initial_condition_vector,0.0)		#Fiber
    PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"] = initial_condition_vector

    PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"] = kinetic_parameter_vector
    
    PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"] = control_parameter_vector
    PROBLEM_DICTIONARY["platelet_count"]=platelet_count

   #platlet controls
	PROBLEM_DICTIONARY["PLATELET_PARAMS"] = platelet_parameter_vector
	PROBLEM_DICTIONARY["PLATELET_PARAMS"][5] = PROBLEM_DICTIONARY["PLATELET_PARAMS"][5]* PROBLEM_DICTIONARY["platelet_count"]; #correct for differing amounts of platelets
    
    # Experimental output scaling -
    # Fig 5A scaling (3.20,1.0)
    # Fig 3 scaling (-1.0,1.0)
    # Fig 2 scaling (1.8)
    # Fig 1 allen scaling (4.0,1.0)
    scaling_parameter_vector =Float64[]
    push!(scaling_parameter_vector,4.0)        # 0 Time scale 
    push!(scaling_parameter_vector,1.0)        # 1 Abundance scale
    PROBLEM_DICTIONARY["SCALING_PARAMETER_VECTOR"] = scaling_parameter_vector

   PROBLEM_DICTIONARY["ALEPH"] = initial_condition_vector[2]+.001
	PROBLEM_DICTIONARY["TIME_DELAY"]=timing
     # QFactor vector - from simulateFig3Butenas2002
    qualitative_factor_vector =Float64[]
   push!(qualitative_factor_vector,2.5)           # 1 TFPI
   push!(qualitative_factor_vector,20.0)          # 2 FV
   push!(qualitative_factor_vector,0.7)           # 3 FVIII
   push!(qualitative_factor_vector,90.0)           # 4 FIX
   push!(qualitative_factor_vector,170.0)         # 5 FX
   push!(qualitative_factor_vector,1.0)           # 6 Platelets
   push!(qualitative_factor_vector, 70)		#7 TAFI
   push!(qualitative_factor_vector,93)		#8 FXIII
    PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"] = qualitative_factor_vector
	PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]=fibrin_kinetic_parameters
	PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]=fibrin_control_parameters

	PROBLEM_DICTIONARY["FVIII_CONTROL"] = 5.20895
	PROBLEM_DICTIONARY["parameter_name_mapping_array"]=["k_trigger","K_trigger" ,"k_amplification" ,"K_FII_amplification" ,"k_APC_formation" ,"K_PC_formation" ,"k_inhibition" ,"K_FIIa_inhibition" ,"k_inhibition_ATIII" ,"k_FV_X_activation", "K_FV_X_actiation" ,"k_FX_activation","k_FX_activation","k_complex","k_amp_prothombinase","K_FII_amp_prothombinase","k_amp_active_factors", "k_amp_active_factors","alpha_trigger_activation","order_trigger_activation","alpha_trigger_inhibition_APC","order_trigger_inhibition_APC","alpha_trigger_inhibition_TFPI", "order_trigger_inhibition_TFPI","alpha_amplification_FIIa","order_amplification_FIIa","alpha_amplification_APC","order_amplification_APC","alpha_amplification_TFPI","order_amplification_TFPI","alpha_shutdown_APC","order_shutdown_APC","alpha_FV_activation","order_FV_activation","alpha_FX_activation","order_FX_activation","alpha_FX_inhibition","order_FX_inhibition","platelet_rate_constant",
"platelet_power","platelet_denominator","Epsmax0","aida","koffplatelets","time_delay","coeff", "platelets"
,"k_cat_Fibrinogen", "Km_Fibrinogen", "k_fibrin_monomer_association", "k_protofibril_association","k_protofibril_monomer_association", "k_cat_plasminogen_Fibrin_tPA", "Km_plasminogen_Fibrin_tPA",
"k_cat_plasminogen_Fibrin_uPA", "Km_plasminogen_Fibrin_uPA", "k_cat_Fibrin", "Km_Fibrin", "k_inhibit_plasmin","k_inhibit_PAI_1_APC", "k_inhibit_tPA_PAI_1", "k_inhibit_uPA_PAI_1",
"k_fibrin_formation", "k_cat_fiber", "Km_fiber", "Ka_plasminogen_Fibrin_tPA","K_plasminogen_Fibrin_tPA","k_cat_fibrinogen_deg", "Km_fibrinogen_deg",
"alpha_fib_inh_fXIII", "order_fib_inh_fXIII", "alpha_fib_inh_TAFI", "order_fib_inh_TAFI", "order_tPA_inh_PAI_1", "alpha_tPA_inh_PAI_1", "alpha_uPA_inh_PAI_1",
"order_uPA_inh_PAI_1"] 

	PROBLEM_DICTIONARY["number_of_states"]=22
    return PROBLEM_DICTIONARY
end
