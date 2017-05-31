include("BalanceEquations.jl")
include("CoagulationModelFactory.jl")
include("utilities.jl")
#using Sundials
using ODE
using NLopt
using POETs


#load data once
#experimentaldata = readdlm("../data/ButenasFig1B60nMFVIIa.csv", ',')
#experimentaldata = readdlm("../data/Luan2010Fig5A.csv", ',')
pathToData = "../data/fromOrfeo_Thrombin_BL_PRP.txt"
data = readdlm(pathToData)
time = data[:,1]
avg_run = mean(data[:,2:3],2);
experimentaldata = hcat(time/60, avg_run)

#pathsToData = ["../data/ButenasFig1B60nMFVIIa.csv","../data/Buentas1999Fig450PercentProthrombin.txt", "../data/Buentas1999Fig4100PercentProthrombin.txt", "../data/Buentas1999Fig4150PercentProthrombin.txt"]
poss_tPA = [0,2]
ids = ["3", "4", "5", "6", "7", "8", "9", "10"]
allexperimentaldata = Array[]
all_platelets = Float64[]
for j in collect(1:size(poss_tPA,1))
	for k in collect(1:size(ids,1))
		platelets,currdata = setROTEMIC(poss_tPA[j], ids[k])
		push!(allexperimentaldata, currdata)
		push!(all_platelets, platelets)
	end
end
selected_idxs = [3,4,5,6,11,12,13,14]

function objectiveForNLOpt(params::Vector, grad::Vector)
	tic()
	dict = buildDictFromOneVector(params)
	TSTART = 0.0
	Ts = .02
	TSTOP = 35.0
	TSIM = collect(TSTART:Ts:TSTOP)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	fbalances(t,y)= BalanceEquations(t,y,dict) 
	#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6,minstep=1E-9)
	t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-4, reltol = 1E-4,minstep=1E-9)
	FIIa = [a[2] for a in X]
	hasdynamics=checkForDynamics(FIIa, t)
	if(hasdynamics)
		MSE, interpolatedExperimentalData=calculateMSE(t, FIIa, allexperimentaldata[curridx])
		#MSE, interpolatedExperimentalData=calculateMSE(t, FIIa, experimentaldata)
	else
		MSE =1E7 #if it doesn't generate dynamics, make this parameter set very unfavorable
	end
	#hold("on")
	#plot(t, FIIa, alpha = .5)
	#write params to file
	f = open("parameterEstimation/NM_05_04_2017.txt", "a+")
	write(f, string(params, ",", MSE, "\n"))
	close(f)
	toc()
	@show MSE
	return MSE
end

function objectiveForPOETS(parameter_array)
	tic()
	obj_array = SharedArray(Float64,8,1)
	#obj_array=10^7*ones(8,1)
	TSTART = 0.0
	Ts = .02
	count = 1
	@show parameter_array
	@sync @parallel for j in selected_idxs
		#@show myid(), j
		temp_params = parameter_array
		temp_params[47] = all_platelets[j] #set platelets to experimental value
		dict = buildCompleteDictFromOneVector(temp_params)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		if(j<10) #no tPA, and set experimental run time
			initial_condition_vector[16]=0.0
			TSTOP = 180.0
			tPA = 0.0
		else
			initial_condition_vector[16]=2.0
			TSTOP = 60.0
			tPA = 2.0
		end
		TSIM = collect(TSTART:Ts:TSTOP)
		fbalances(t,y)= BalanceEquations(t,y,dict) 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		FIIa = [a[2] for a in X]
		A = convertToROTEM(t,X,tPA)
		hasdynamics=checkForDynamics(FIIa, t)
		if(hasdynamics)
			print("has dynamics")
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end
		@show myid(), count,MSE
		obj_array[findin(selected_idxs,j),1]=MSE
		count = count+1
		@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	toc()
	return obj_array
end

function attemptOptimizationPOETS()
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 8
	initial_parameter_estimate = vec(readdlm("parameterEstimation/fittingPatient6_tPA2_05_24_17.txt", ','))
	#inital_parameter_estimate= readdlm("parameterEstimation/paramsToRestart03_30_2017.txt")
	outputfile = "parameterEstimation/POETS_info_24_05_2017maxstep1_OriginalShapeFunction.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*1.05, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/1.05, initial_parameter_estimate[47:end]/1000)
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		initial_parameter_estimate =initial_parameter_estimate+initial_parameter_estimate*rand()*.1

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveForPOETS,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)

		# Package -
		@show (EC, PC, RA)
		ec_array = [ec_array EC]
		pc_array = [pc_array PC]
		f = open(outputfile, "a")
		write(f, string(EC, ",", PC, ",", RA, "\n"))
		close(f)
	end

	return (ec_array,pc_array)
end


function attemptOptimizationSingleAndMulti()
	number_of_subdivisions = 10
	number_of_parameters = 46
	number_of_objectives = 2
	platelet_counts = [420,364]
	initial_parameter_estimate = readdlm("parameterEstimation/BestAfterNM_05_04_2017.txt")
	outputfile = "parameterEstimation/POETS_info_05_04_2017_SingleAndMulti.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	initial_parameter_estimate =vec(initial_parameter_estimate) #get initial parameter estimate
	global up_arr = initial_parameter_estimate*1000
	global lb_arr = initial_parameter_estimate/1000
	global curridx = 1
	for index in collect(1:number_of_subdivisions)

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveForPOETS,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)
		#get the index of the parameter set with smallest mean error (on both sets)
		best_index=indmax(mean(EC,1))
		best_params = PC[:,best_index]

		# Package -
		ec_array = [ec_array EC]
		pc_array = [pc_array PC]
		@show (EC, PC, RA)
		f = open(outputfile, "a")
		write(f, string(EC, ",", PC, ",", RA, "\n"))
		close(f)
		#estimate params for single objective, alternating which platelet count to use
		curridx =mod(index,2)+1 
		inital_parameter_estimate=attemptOptimizationNLOpt(best_params, platelet_counts[curridx])  
	end

	return (ec_array,pc_array)
end

function attemptOptimizationNLOpt()
	numvars = 47
	platelet_count = 420
	opt = Opt(:LN_NELDERMEAD,numvars)
	lowerbounds =fill(0, 1, numvars)
	upperbounds = fill(1E7, 1, numvars)
	upperbounds[3] = 70.0 #bound k_amplication to be small
	upperbounds[17] = 70.0 #bound k_amp_active_factors to be small
	upperbounds[47]=platelet_count #make it so platelet count can't move
	lowerbounds[47]=platelet_count
	upper_bounds!(opt, vec(upperbounds))
	lower_bounds!(opt, vec(lowerbounds))
	min_objective!(opt, objectiveForNLOpt)
	#ftol_abs!(opt, 1E-4)
	#min_objective!(opt, objectiveForNLOptFourDataSets)
#		
#		# Kinetic parameters -
#	 kinetic_parameter_vector = Float64[]
#		push!(kinetic_parameter_vector,7200/2)			 #	k_trigger
#		push!(kinetic_parameter_vector,100)			 # 1 K_trigger
#		push!(kinetic_parameter_vector,1.0)				# 2 k_amplification
#		push!(kinetic_parameter_vector,12)			 # 3 K_FII_amplification
#		push!(kinetic_parameter_vector,.5)				# 4 k_APC_formation
#		push!(kinetic_parameter_vector,3/100.0)				 # 5 K_PC_formation
#		push!(kinetic_parameter_vector,0.2*10)				# 6 k_inhibition
#		push!(kinetic_parameter_vector,120)			 # 7 K_FIIa_inhibition
#		push!(kinetic_parameter_vector,6E-4)		 # 8 k_inhibition_ATIII
#		#push!(kinetic_parameter_vector,0.001)			# 9 K_inhibition_ATIII
#		#push!(kinetic_parameter_vector,100.0)			# 10 K_inhibition_FIIa
#		push!(kinetic_parameter_vector, 2E7*60*10.0^-4) #9 k_FV_activation, from reaction 16 in Diamond 2010 paper
#		push!(kinetic_parameter_vector, 1E8*80*10.0^-6/100) #10 K_FV_activation 
#		push!(kinetic_parameter_vector, 600.0) #11 k_FX_activation from reaction 6 in Diamond 2010 paper
#		push!(kinetic_parameter_vector, .28) #12 K_FX_activation
#		push!(kinetic_parameter_vector, 24) #13 k_complex
#		push!(kinetic_parameter_vector, 63.5*60*2.5 )#14 k_amp_prothombinase from reaction 18 in Diamond 2010
#		push!(kinetic_parameter_vector,10) #13 K_FII_amp_prothombinase
#		push!(kinetic_parameter_vector, 6.0) #k_amp_active_factors
#		push!(kinetic_parameter_vector, 1.0) #K_amp_active_factor

#		
#		# Control parameters -
#		control_parameter_vector =Float64[]
#		# Trigger -
#		push!(control_parameter_vector,140.0/10)			# 0 9 alpha_trigger_activation = control_parameter_vector[0]
#		push!(control_parameter_vector,2.0)				# 1 10 order_trigger_activation = control_parameter_vector[1]
#		push!(control_parameter_vector,1.0)				# 2 11 alpha_trigger_inhibition_APC = control_parameter_vector[2]
#		push!(control_parameter_vector,1.0)				# 3 12 order_trigger_inhibition_APC = control_parameter_vector[3]
#		push!(control_parameter_vector,0.1)				# 4 13 alpha_trigger_inhibition_TFPI = control_parameter_vector[4]
#		push!(control_parameter_vector,2.0)				# 5 14 order_trigger_inhibition_TFPI = control_parameter_vector[5]
#		
#		# Amplification -
#		push!(control_parameter_vector,0.1)				# 6 15 alpha_amplification_FIIa = control_parameter_vector[6]
#		push!(control_parameter_vector,2.0)				# 7 16 order_amplification_FIIa = control_parameter_vector[7]
#		push!(control_parameter_vector,0.4)				# 8 17 alpha_amplification_APC = control_parameter_vector[8]
#		push!(control_parameter_vector,2.0)				# 9 18 order_amplification_APC = control_parameter_vector[9]
#		push!(control_parameter_vector,0.01)			 # 10 19 alpha_amplification_TFPI = control_parameter_vector[10]
#		push!(control_parameter_vector,2.0)				# 11 20 order_amplification_TFPI = control_parameter_vector[11]
#		
#		# APC generation -
#		push!(control_parameter_vector,2.0)				# 12 21 alpha_shutdown_APC = control_parameter_vector[12]
#		push!(control_parameter_vector,2.0)				# 13 22 order_shutdown_APC = control_parameter_vector[13]

#		push!(control_parameter_vector,1.0) #14 alpha_FV_activation
#		push!(control_parameter_vector,1.0/10) #15 order_FV_activation
#		push!(control_parameter_vector,1.0) #16 alpha_FX_activation
#		push!(control_parameter_vector,1.0) #17 order_FX_activation
#		push!(control_parameter_vector,4.0) #14 alpha_FX_inhibition
#		push!(control_parameter_vector,1.0) #15 order_FX_inhibition
#	 #platlet controls
#	platelet_parameter_vector = Float64[]
#	normal_platelet_count = 300 #*10^6 #/mL
#	push!(platelet_parameter_vector, .005*5) #1 rate constant
#	push!(platelet_parameter_vector, 1.6123) #2 power for control function
#	push!(platelet_parameter_vector, 2.4279E-9) #3 adjustment in denominator
#	push!(platelet_parameter_vector, .01) #4 Epsmax0
#	push!(platelet_parameter_vector, 1.05*platelet_count/normal_platelet_count)	#5 aida
#	push!(platelet_parameter_vector, .005) #koffplatelets

#	timing = Float64[]
#	push!(timing, 1.0) #time_delay
#	push!(timing, 1.0) #coeff
#	 
	#inital_parameter_estimate = vcat(kinetic_parameter_vector, control_parameter_vector, platelet_parameter_vector, timing, platelet_count)	
	inital_parameter_estimate = readdlm("parameterEstimation/startingPoint05_04_2017.txt")	
	@show inital_parameter_estimate
	(minf, minx, ret) = NLopt.optimize(opt, vec(inital_parameter_estimate))
	println("got $minf at $minx after $count iterations (returned $ret)")
end

function attemptOptimizationNLOpt(parameter_estimate, platelet_count)
	numvars = 47
	outputfile = "parameterEstimation/DuringSingleObj_04_04_2017.txt"
	opt = Opt(:LN_NELDERMEAD,numvars)
	lowerbounds =vec(lb_arr[1:46])
	upperbounds = vec(up_arr[1:46])
	#upperbounds[47]=platelet_count #make it so platelet count can't move
	#lowerbounds[47]=platelet_count
	push!(upperbounds, platelet_count)
	push!(lowerbounds, platelet_count)
	min_objective!(opt, objectiveForNLOpt)
	inital_parameter_estimate = parameter_estimate[1:46]
	push!(inital_parameter_estimate, platelet_count)
	
	@show size(inital_parameter_estimate)
	@show lowerbounds
	@show upperbounds
	@show inital_parameter_estimate
	upper_bounds!(opt, vec(upperbounds))
	lower_bounds!(opt, vec(lowerbounds))
	(minf, minx, ret) = NLopt.optimize(opt, vec(inital_parameter_estimate))
	println("got $minf at $minx after $count iterations (returned $ret)")
	f = open(outputfile, "a")
	write(f, string(minf, ",", minx, ",", ret, "\n"))
	close(f)
	return minx
end


