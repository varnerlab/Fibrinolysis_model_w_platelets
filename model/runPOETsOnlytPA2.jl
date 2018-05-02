include("Balances.jl")
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
selected_idxs = [11,12,13,14]

function objectiveForPOETS(parameter_array)
	tic()
	obj_array = SharedArray(Float64,4,1)
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
		#@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	toc()
	return obj_array
end

function objectiveForPOETSPlatletContribution(parameter_array)
	tic()
	obj_array = SharedArray(Float64,4,1)
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
		fbalances(t,y)= Balances(t,y,dict) 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		FIIa = [a[2] for a in X]
		fibrinogen = [a[14] for a in X]
		A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[j])
		hasdynamics=checkForDynamics(FIIa, t)
		if(hasdynamics)
			print("has dynamics")
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end
		#check to make sure we used up fibrinogen, penalize if we haven't
		if(fibrinogen[end]<100)
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE = 10^7
		end

		@show myid(), count,MSE
		obj_array[findin(selected_idxs,j),1]=MSE
		count = count+1
		#@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	toc()
	return obj_array
end

function objectiveForPOETSPlatletContributionDiffROTEM(parameter_array)
	tic()
	obj_array = SharedArray(Float64,4,1)
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
		fbalances(t,y)= Balances(t,y,dict) 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		FIIa = [a[2] for a in X]
		A = convertToROTEMPlateletContributionScaledF(t,X,tPA,all_platelets[j],initial_condition_vector[14])
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
		#@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	toc()
	return obj_array
end

function attemptOptimizationPOETSOnlytPA2()
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 4
	initial_parameter_estimate = vec(readdlm("parameterEstimation/fittingPatient6_tPA2_05_24_17.txt", ','))
	#inital_parameter_estimate= readdlm("parameterEstimation/paramsToRestart03_30_2017.txt")
	outputfile = "parameterEstimation/POETS_info_25_05_2017maxstep1_OriginalShapeFunction_OnlytPA2.txt"
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

function attemptOptimizationPOETSOnlytPA2PlateletContribution()
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 4
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/handfitting_23_3_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/goodParamsForPatient7_28_3_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/meanParamsStartPoint_05_04_18.txt"))
	#initial_parameter_estimate=vec(readdlm("../parameterEstimation/useUpFibrinogen_04_11_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/startingPoint_19_4_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_25_04_18.txt"))
	initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_29_04_18.txt"))
	outputfile = "../parameterEstimation/POETS_info_29_04_18_PlateletContributionToROTEM.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*1.005, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/1.005, initial_parameter_estimate[47:end]/1000)
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		initial_parameter_estimate =initial_parameter_estimate+initial_parameter_estimate*rand()*.1

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveForPOETSPlatletContribution,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)

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


function attemptOptimizationPOETSOnlytPA2PlateletContributionDiffROTEM()
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 4
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/handfitting_23_3_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/goodParamsForPatient7_28_3_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/startingPoint_09_04_18_useupallfibrinogen.txt"))
	outputfile = "../parameterEstimation/POETS_info_09_04_18_PlateletContributionToROTEMDiffROTEM.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*1.05, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/1.05, initial_parameter_estimate[47:end]/1000)
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		initial_parameter_estimate =initial_parameter_estimate+initial_parameter_estimate*rand()*.1

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveForPOETSPlatletContribution,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)

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
