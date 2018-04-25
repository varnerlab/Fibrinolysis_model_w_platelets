include("CoagulationModelFactory.jl")
include("utilities.jl")

using DifferentialEquations
using Plots

function BalanceEquations(dxdt_total,x,p,t)
	num_rates = 21
	# Correct nagative x's = throws errors in control even if small - 
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
	kinetic_parameter_vector = p[1:18]#PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"]
	control_parameter_vector = p[19:38]#PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"]
	qualitative_factor_level_vector =p[39:46] #PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"]
	platelet_parameter_vector = p[47:52]#PROBLEM_DICTIONARY["PLATELET_PARAMS"]
	timing =p[53:54] #PROBLEM_DICTIONARY["TIME_DELAY"]
	fibrin_kinetic_parameter_vector =p[55:76] #PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]
	fibrin_control_vector =p[77:84] #PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]

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
	FVIII_control = p[85]#PROBLEM_DICTIONARY["FVIII_CONTROL"]

	#platelet control
	#update aleph so that it holds the maximum value of FIIa
	if(FIIa>p[86])
		p[86]=FIIa
	end

	aleph = p[86]
	faleph = aleph^platelet_pwr/(aleph^platelet_pwr + platelet_denom^platelet_pwr)
	EpsMax = EpsMax0+(1-EpsMax0)*faleph
	#@show aleph, falpeh, EpsMax

	#timing
	time_delay = timing[1]
	time_coeff = timing[2]

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

    control_vector = Array{typeof(initiation_trigger_term),2}(1,num_rates)
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
	#@show control_vector

    # Calculate the kinetics -
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
	#@show FIIa, FII, PROTHOMBINASE_PLATELETS, FV_FXA

    rate_vector =  Array{typeof(initiation_trigger_term),2}(1,num_rates)
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
  	rate_vector[21] = (Plasmin*k_cat_fibrinogen_deg*Fibrinogen)/(Km_fibrinogen_deg+Fibrinogen)

	for j in collect(1:length(rate_vector))
		if(isnan(rate_vector[j]) || !isfinite(rate_vector[j]))
		#	println("Found NAN at rate_vector $j")
			rate_vector[j] = 0.0
		end
	end

	#remove nans
	for j in collect(1:length(rate_vector))
		if(isnan(control_vector[j])|| !isfinite(control_vector[j]))
		#	println("Found NAN at control_vector $j")
			control_vector[j] = 1.0
		end
	end

	 # modified rate vector -
	modified_rate_vector = (rate_vector).*(control_vector);
	
	dxdt_total[1] = -1*modified_rate_vector[2] -modified_rate_vector[7]-modified_rate_vector[6]	# 1 FII
	dxdt_total[2] = modified_rate_vector[2] - modified_rate_vector[4]+modified_rate_vector[7]+modified_rate_vector[6]	 # 2 FIIa
	dxdt_total[3] = -1*modified_rate_vector[3] # 3 PC
	dxdt_total[4] = 1*modified_rate_vector[3]  # 4 APC
	dxdt_total[5] = -modified_rate_vector[4]# 5 ATIII
	dxdt_total[6] = 0.0 # 6 TM (acts as enzyme, so no change)
	dxdt_total[7] = -0.0*TRIGGER	# 7 TRIGGER
	dxdt_total[8] = kplatelts*(EpsMax-Eps)-koffplatelets*Eps#-rate_vector[5] #frac active platelets
	#dxdt_total[8] = 0.0 #if no platelets, they can't be activated
	dxdt_total[9] = -1*modified_rate_vector[1] #FV_FX
	dxdt_total[10] = modified_rate_vector[1] -modified_rate_vector[5]#FV_FXa
	dxdt_total[11] = modified_rate_vector[5] #Prothromibase_platelets
	dxdt_total[12] = modified_rate_vector[12]-modified_rate_vector[15]                     # 8 Fibrin
  	dxdt_total[13] = modified_rate_vector[13]+modified_rate_vector[14]-modified_rate_vector[16]     # 9 Plasmin
	dxdt_total[14] = -modified_rate_vector[8] - modified_rate_vector[21]                 # 10 Fibrinogen
	dxdt_total[15] = -modified_rate_vector[13] - modified_rate_vector[14]                # 11 Plasminogen
	dxdt_total[16] = -modified_rate_vector[18]                                 # 12 tPA
	dxdt_total[17] = -modified_rate_vector[19]                                 # 13 uPA
	dxdt_total[18] = 1.0(modified_rate_vector[8]-modified_rate_vector[9])                 # 14 Fibrin monomer
	dxdt_total[19] = modified_rate_vector[9]+modified_rate_vector[11]-modified_rate_vector[10]       # 15 Protofibril
	dxdt_total[20] = -modified_rate_vector[16]                                 # 16 Antiplasmin
	dxdt_total[21] = -modified_rate_vector[17]-modified_rate_vector[18]-modified_rate_vector[19]   # 17 PAI_1
	dxdt_total[22] = 1.0*(modified_rate_vector[10]-modified_rate_vector[12]-modified_rate_vector[20])# 18 Fibers

	#idx = find(x->(x<0),x);
	#x[idx] = 0.0;
	#dxdt_total[idx]= 0.0


	tau = time_coeff*(1-FIIa/aleph)
	time_scale =1-1*exp(-tau*(t-time_delay))
	#@show 1-1*exp(-tau*(t-time_delay))
	if(t<time_delay)
		time_scale = 0.0
		#@show t, time_scale, time_delaye
	end
	#idx = find(dxdt_total->(abs(dxdt_total)<1E-9),dxdt_total); #lets make anything less than 1E-9 zero
	#@show idx
	#dxdt_total[idx] = 0.0;

	dxdt_total = dxdt_total.*time_scale
end

function runSim()
	tstart = 0.0
	tend = 35.0
	PROBLEM_DICTIONARY = buildCoagulationModelDictionary(346.0)
	initial_condition_vector = PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"]
	#pull things out of dictionary
	kinetic_parameter_vector = PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"]
	control_parameter_vector = PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"]
	qualitative_factor_level_vector = PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"]
	platelet_parameter_vector = PROBLEM_DICTIONARY["PLATELET_PARAMS"]
	timing = PROBLEM_DICTIONARY["TIME_DELAY"]
	fibrin_kinetic_parameter_vector = PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]
	fibrin_control_vector = PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]
	FVIII_control = PROBLEM_DICTIONARY["FVIII_CONTROL"]
	aleph=PROBLEM_DICTIONARY["ALEPH"]

	allparams = vcat(kinetic_parameter_vector,control_parameter_vector,qualitative_factor_level_vector,platelet_parameter_vector,timing,fibrin_kinetic_parameter_vector,fibrin_control_vector, FVIII_control,aleph)
	#pf = ParameterizedFunction(pf_func, allparams)
	#prob = ODELocalSensitivityProblem(pf, initial_condition_vector[1:11] , (tstart, tend))
	#sol = solve(prob,alghints=[:stiff],dt = .02, abstol=1E-6, reltol=1E-6)
	prob=ODEProblem(BalanceEquations, initial_condition_vector, (tstart,tend), allparams)
	sol = solve(prob, Tsit5(),reltol=1e-8,abstol=1e-8)
	Plots.plot(sol)
	#tpa is 12
	A = convertToROTEM(sol.t, sol, initial_condition_vector[12])
	PyPlot.plot(sol.t, A)
	CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(A,sol.t)
	return sol
end

function runSim(params)
	tstart = 0.0
	tend = 35.0
	PROBLEM_DICTIONARY = buildCompleteDictFromOneVector(params)
	initial_condition_vector =PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"]
	#pull things out of dictionary
	kinetic_parameter_vector = PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"]
	control_parameter_vector = PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"]
	qualitative_factor_level_vector = PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"] #this is actually set, not adjustable, but need to feed it in
	platelet_parameter_vector = PROBLEM_DICTIONARY["PLATELET_PARAMS"]
	timing = PROBLEM_DICTIONARY["TIME_DELAY"]
	fibrin_kinetic_parameter_vector = PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]
	fibrin_control_vector = PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]
	FVIII_control = PROBLEM_DICTIONARY["FVIII_CONTROL"]
	aleph=PROBLEM_DICTIONARY["ALEPH"]

	allparams = vcat(kinetic_parameter_vector,control_parameter_vector,qualitative_factor_level_vector,platelet_parameter_vector,timing,fibrin_kinetic_parameter_vector,fibrin_control_vector, FVIII_control,aleph)
	#pf = ParameterizedFunction(pf_func, allparams)
	#prob = ODELocalSensitivityProblem(pf, initial_condition_vector[1:11] , (tstart, tend))
	#sol = solve(prob,alghints=[:stiff],dt = .02, abstol=1E-6, reltol=1E-6)
	prob=ODEProblem(BalanceEquations, initial_condition_vector, (tstart,tend), allparams)
	sol = solve(prob,  CVODE_BDF(method=:Functional),reltol=1e-8,abstol=1e-8, dtmax=1.0,saveat = 0.01)
	Plots.plot(sol)
	#tpa is 12
	A = convertToROTEM(sol.t, sol, initial_condition_vector[12])
	PyPlot.plot(sol.t, A)
	@show CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(A,sol.t)
	return sol
end

