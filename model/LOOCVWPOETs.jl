using Distributed
@everywhere include("Balances.jl")
@everywhere include("Kinetics.jl")
@everywhere include("Control.jl")
@everywhere include("CoagulationModelFactory.jl")
@everywhere include("utilities.jl")
@everywhere include("LOOCVutils.jl")
#using Sundials
using DifferentialEquations
using NLopt
using POETs
using SharedArrays
using DelimitedFiles
using Statistics


#load data once
#experimentaldata = readdlm("../data/ButenasFig1B60nMFVIIa.csv", ',')
#experimentaldata = readdlm("../data/Luan2010Fig5A.csv", ',')
pathToData = "../data/fromOrfeo_Thrombin_BL_PRP.txt"
data = readdlm(pathToData)
time = data[:,1]
avg_run = mean(data[:,2:3],dims=2);
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

function objectiveForPOETSPlatletContributionOnlyTPA2_leaveOnOutFromTraining(parameter_array)
	obj_array = SharedArray{Float64}(3,1)
	#obj_array=10^7*ones(8,1)
	TSTART = 0.0
	Ts = .02
	count = 1
	@show parameter_array
	@sync @distributed for j in selected_idxs
	#for j in selected_idxs
		@show myid(), j
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
		#fbalances(t,y)= Balances(t,y,dict) 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		#t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		fbalances(y,p,t)= Balances(t,y,dict) 
		#fbalances(t,y)= Balances(t,y,dict) 
		#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
		prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
		@time sol = solve(prob)
		t =sol.t
		X = sol
		FIIa = X[2,:]
		fibrinogen = X[14,:]
		A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[j])
		AnyFlats=checkForFlatness(t,A)
		hasdynamics=checkForDynamics(FIIa, t)
		#make sure it has dynamics, used up fibrinogen and doesn't have any flat patches'
		#@show hasdynamics, fibrinogen[end], AnyFlats
		if(hasdynamics && fibrinogen[end]<370 && AnyFlats==false)
			print("has dynamics")
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end
		#check to make sure we used up fibrinogen, penalize if we haven't
		@show myid(), count,MSE
		#obj_array[findin(selected_idxs,j),1]=MSE #old, julia .6
		obj_array[findall(x->x==j,selected_idxs)[1],1]=MSE
		count = count+1
		#@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	return obj_array
end


function objectiveForPOETSPlatletContribution(parameter_array)
	obj_array = SharedArray{Float64}(7,1)
	#obj_array=10^7*ones(8,1)
	TSTART = 0.0
	Ts = .02
	count = 1
	@show parameter_array
	@sync @distributed for j in selected_idxs
	#for j in selected_idxs
		@show myid(), j
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
		AnyFlats=checkForFlatness(t,A)
		hasdynamics=checkForDynamics(FIIa, t)
		#make sure it has dynamics, used up fibrinogen and doesn't have any flat patches'
		#@show hasdynamics, fibrinogen[end], AnyFlats
		if(hasdynamics && fibrinogen[end]<370 && AnyFlats==false)
			print("has dynamics")
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end
		#check to make sure we used up fibrinogen, penalize if we haven't
		@show myid(), count,MSE
		#obj_array[findin(selected_idxs,j),1]=MSE #old, julia .6
		obj_array[findall(x->x==j,selected_idxs)[1],1]=MSE
		count = count+1
		#@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	return obj_array
end

function objectiveForPOETSPlatletContributionBothCases(parameter_array)
	tic()
	obj_array = SharedArray{Float64}(7,1)
	#obj_array=10^7*ones(8,1)
	TSTART = 0.0
	Ts = .02
	count = 1
	@show parameter_array
	@distributed for j in selected_idxs
		#@show myid(), j
		temp_params = parameter_array
		temp_params[47] = all_platelets[j] #set platelets to experimental value
		dict = buildCompleteDictFromOneVector(temp_params)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		if(j<10) #no tPA, and set experimental run time
			initial_condition_vector[16]=0.0
			TSTOP = 180.0
			tPA = 0.0
			initial_condition_vector[16]=tPA
		else
			initial_condition_vector[16]=2.0
			TSTOP = 60.0
			tPA = 2.0
		end
		#run tPA = 2 case
		TSIM = collect(TSTART:Ts:TSTOP)
		fbalances(t,y)= Balances(t,y,dict) 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		FIIa = [a[2] for a in X]
		fibrinogen = [a[14] for a in X]
		A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[j])
		AnyFlats=checkForFlatness(t,A)
		hasdynamics=checkForDynamics(FIIa, t)
		#make sure it has dynamics, used up fibrinogen and doesn't have any flat patches'
		#@show hasdynamics, fibrinogen[end], AnyFlats
		if(hasdynamics && fibrinogen[end]<370 && AnyFlats==false)
			print("has dynamics")
			MSE2, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE2 =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end
		#run tPA=0 case
		TSTOP = 180.0
		TSIM = collect(TSTART:Ts:TSTOP)
		fbalances(t,y)= Balances(t,y,dict)
		tPA = 0.0
		initial_condition_vector[16]=tPA 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		FIIa = [a[2] for a in X]
		fibrinogen = [a[14] for a in X]
		A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[j])
		AnyFlats=checkForFlatness(t,A)
		hasdynamics=checkForDynamics(FIIa, t)
		#make sure it has dynamics, used up fibrinogen and doesn't have any flat patches'
		#@show hasdynamics, fibrinogen[end], AnyFlats
		if(hasdynamics && fibrinogen[end]<370 && AnyFlats==false)
			print("has dynamics")
			MSE0, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE0 =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end



		#check to make sure we used up fibrinogen, penalize if we haven't
		@show myid(), count,MSE0, MSE2
		obj_array[findin(selected_idxs,j),1]=MSE2+MSE0
		count = count+1
		#@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	toc()
	return obj_array
end

function objectiveForPOETSPlatletContributionBothCases14D(parameter_array)
	tic()
	obj_array = SharedArray{Float64}(14,1)
	#obj_array=10^7*ones(8,1)
	TSTART = 0.0
	Ts = .02
	count = 1
	@show parameter_array
	@distributed for j in selected_idxs
		#@show myid(), j
		temp_params = parameter_array
		temp_params[47] = all_platelets[j] #set platelets to experimental value
		dict = buildCompleteDictFromOneVector(temp_params)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		if(j<10) #no tPA, and set experimental run time
			TSTOP = 180.0
			tPA = 0.0
			initial_condition_vector[16]=tPA
		else
			tPA = 2.0
			initial_condition_vector[16]=2.0
			TSTOP = 60.0
			tPA = 2.0
		end
		#run tPA = 2 case
		TSIM = collect(TSTART:Ts:TSTOP)
		fbalances(t,y)= Balances(t,y,dict) 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		FIIa = [a[2] for a in X]
		fibrinogen = [a[14] for a in X]
		A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[j])
		AnyFlats=checkForFlatness(t,A)
		hasdynamics=checkForDynamics(FIIa, t)
		#make sure it has dynamics, used up fibrinogen and doesn't have any flat patches'
		println("tPA = 2")
		@show hasdynamics, fibrinogen[end], AnyFlats
		if(hasdynamics && fibrinogen[end]<370 && AnyFlats==false)
			#print("has dynamics")
			MSE2, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE2 =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end
		#run tPA=0 case
		TSTOP = 180.0
		TSIM = collect(TSTART:Ts:TSTOP)
		fbalances(t,y)= Balances(t,y,dict)
		tPA = 0.0
		initial_condition_vector[16]=tPA 
		#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
		FIIa = [a[2] for a in X]
		fibrinogen = [a[14] for a in X]
		A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[j])
		AnyFlats=checkForFlatness(t,A)
		hasdynamics=checkForDynamics(FIIa, t)
		#make sure it has dynamics, used up fibrinogen and doesn't have any flat patches'
		println("tPA = 0")
		@show hasdynamics, fibrinogen[end], AnyFlats
		#for tPA =0, we expect some flatness, so let's not penalize parameter sets that have flatness
		if(hasdynamics && fibrinogen[end]<370)
			#print("has dynamics")
			MSE0, interpData = calculateMSE(t,A, allexperimentaldata[j])
		else
			MSE0 =10^7 #if it doesn't generate dynamics, make this parameter set very unfavorable
		end



		#check to make sure we used up fibrinogen, penalize if we haven't
		@show myid(), count,MSE0, MSE2
		obj_array[findin(selected_idxs,j),1]=MSE0
		obj_array[findin(selected_idxs,j)+7,1]=MSE2
		count = count+1
		#@show obj_array
	end
	@show obj_array
	#@show size(parameter_array)
	toc()
	return obj_array
end

function attemptOptimizationPOETSOnlytPA2PlateletContribution(outputfilename, sel_idxs)
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 7
	global selected_idxs = sel_idxs
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/handfitting_23_3_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/goodParamsForPatient7_28_3_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/meanParamsStartPoint_05_04_18.txt"))
	#initial_parameter_estimate=vec(readdlm("../parameterEstimation/useUpFibrinogen_04_11_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/startingPoint_19_4_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_25_04_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_29_04_18.txt"))
	initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_02_05_18.txt"))
	outputfile = outputfilename
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

function attemptOptimizationPOETSOnlytPA2PlateletContributionRestartable(outputfilename, sel_idxs)
	number_of_subdivisions = 5
	number_of_parameters = 77
	number_of_objectives = 7
	global selected_idxs = sel_idxs
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/handfitting_23_3_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/goodParamsForPatient7_28_3_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/meanParamsStartPoint_05_04_18.txt"))
	#initial_parameter_estimate=vec(readdlm("../parameterEstimation/useUpFibrinogen_04_11_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/startingPoint_19_4_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_25_04_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_29_04_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_02_05_18.txt"))
	initial_parameter_estimate, round = checkForPreviousAndLoad(outputfilename)
	@show round, size(initial_parameter_estimate)

	#append round to output file name
	outputfilename = string(outputfilename, "Round_", round, ".txt")
	@show outputfilename
	outputfile = outputfilename
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*1.005, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/1.005, initial_parameter_estimate[47:end]/1000)
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		initial_parameter_estimate =initial_parameter_estimate#+initial_parameter_estimate*rand()*.1

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveForPOETSPlatletContribution,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=5,show_trace=true)

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

function attemptOptimizationPOETSOnlytPA2PlateletContributionRestartable_LeaveOneOutFromTraining(outputfilename, sel_idxs)
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 3
	max_iters = 5
	global selected_idxs = sel_idxs
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/handfitting_23_3_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/goodParamsForPatient7_28_3_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/meanParamsStartPoint_05_04_18.txt"))
	#initial_parameter_estimate=vec(readdlm("../parameterEstimation/useUpFibrinogen_04_11_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/startingPoint_19_4_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_25_04_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_29_04_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_02_05_18.txt"))
	initial_parameter_estimate, round = checkForPreviousAndLoad(outputfilename)
	@show round, size(initial_parameter_estimate)

	#append round to output file name
	outputfilename = string(outputfilename, "Round_", round,"Max_iters",number_of_subdivisions, ".txt")
	@show outputfilename
	outputfile = outputfilename
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*1.005, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/1.005, initial_parameter_estimate[47:end]/1000)
	#running into undefined variable errors-will this fix it?
	EC = []
	PC = []
	for index in collect(1:number_of_subdivisions)
		output_eachiter =string(outputfilename, "Round_", index, "outOf", number_of_subdivisions, ".txt")
		# Grab a starting point -
		if index ==1
			initial_parameter_estimate =initial_parameter_estimate#+initial_parameter_estimate*rand()*.1
		else
			total_error = sum(EC[:,1:end],dims=1)
			min_index = argmin(total_error)[2]
			@show size(PC)
			@show min_index
			initial_parameter_estimate= PC[:,min_index]
		end

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveForPOETSPlatletContributionOnlyTPA2_leaveOnOutFromTraining,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=max_iters,show_trace=true)

		# Package -
		@show (EC, PC, RA)
		ec_array = [ec_array EC]
		pc_array = [pc_array PC]
		f = open(outputfile, "a")
		write(f, string(EC, ",", PC, ",", RA, "\n"))
		close(f)
		f = open(output_eachiter, "a")
		write(f, string(EC, ",", PC, ",", RA, "\n"))
		close(f)
	end

	return (ec_array,pc_array)
end


function attemptOptimizationPOETSBothCasesPlateletContributionRestartable(outputfilename, sel_idxs)
	number_of_subdivisions = 10
	number_of_parameters = 77
	#number_of_objectives = 7
	number_of_objectives = 14
	global selected_idxs = sel_idxs
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/handfitting_23_3_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/goodParamsForPatient7_28_3_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/meanParamsStartPoint_05_04_18.txt"))
	#initial_parameter_estimate=vec(readdlm("../parameterEstimation/useUpFibrinogen_04_11_18.txt"))
	#initial_parameter_estimate = vec(readdlm("../parameterEstimation/startingPoint_19_4_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_25_04_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_29_04_18.txt"))
	#initial_parameter_estimate= vec(readdlm("../parameterEstimation/startingPoint_02_05_18.txt"))
	initial_parameter_estimate, round = checkForPreviousAndLoad(outputfilename)
	@show round, size(initial_parameter_estimate)

	#append round to output file name
	outputfilename = string(outputfilename, "Round_", round, ".txt")
	@show outputfilename
	outputfile = outputfilename
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*1.005, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/1.005, initial_parameter_estimate[47:end]/1000)
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		initial_parameter_estimate =initial_parameter_estimate#+initial_parameter_estimate*rand()*.1

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveForPOETSPlatletContributionBothCases14D,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=5,show_trace=true)

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

#allids = collect(9:16)
#leaveoutidx = parse(Int64, ARGS[1])
#leaveoutid = allids[leaveoutidx]
#selids = [allids[1:leaveoutidx-1] ;allids[ leaveoutidx+1:end]]
#outputfilename = string("../LOOCV/POETS_info_23_01_19_PlateletContributionToROTEMFlatness1ToBeTestedOn",strip(string(leaveoutid)))
#@show outputfilename
#attemptOptimizationPOETSOnlytPA2PlateletContributionRestartable(outputfilename, selids)
#attemptOptimizationPOETSBothCasesPlateletContributionRestartable(outputfilename, selids)

training_ids = collect(11:14)
leaveoutidx=1
leaveoutid=training_ids[leaveoutidx]
selids=[training_ids[1:leaveoutidx-1];training_ids[leaveoutidx+1:end] ]
outputfilename = string("../LOOCV/POETS_info_14_02_19_PlateletContributionToROTEMFlatness1ToBeTestedOn",strip(string(leaveoutid)))
attemptOptimizationPOETSOnlytPA2PlateletContributionRestartable_LeaveOneOutFromTraining(outputfilename, selids)
