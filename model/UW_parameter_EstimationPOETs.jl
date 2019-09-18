using POETs
using CSV

include("runModel.jl")
include("inverseProbUtils.jl")
include("LOOCVutils.jl")

#load data once
pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
all_data = CSV.read(pathToData)

function plotAllResults()
	params = readdlm("../parameterEstimation/Best24_UW_Rescaled_16_09_19.txt")
	ids = [855, 856, 880]
	hours = [0,36]
	temp_params = params[1,:]
	dict = buildCompleteDictFromOneVector(temp_params)
	nominal_ICs= dict["INITIAL_CONDITION_VECTOR"]
	for h in hours
		for id in ids
			plotResult(id, h, nominal_ICs, all_data, params)
		end
	end
end

function plotResult(id::Int, hour::Int, nominal_ICs, all_data, params)
	close("all")
	initial_condition_vector = updateICsForPatient(id, hour, all_data, nominal_ICs)
	sel_target  = constructTarget(hour, id, all_data)
	figure(figsize =[15,15])
	num_param_sets = size(params,1)
	TSTART =0.0
	TSTOP = 60.0
	curr_platelets=300.0
	@show id, hour
	@show sel_target
	for j in 1:num_param_sets
		curr_params = params[j,:]
		dict = buildCompleteDictFromOneVector(curr_params)
		fbalances(y,p,t)= Balances(t,y,dict) 
		prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
		sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = .1,maxiters = 1e6)
		t =sol.t
		X = sol
		#A = convertToROTEM(t,X,tPA)
		tPA = initial_condition_vector[16]
		A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
		plot(t,A, "grey")
		metrics = calculateCommonMetrics(A,t)
		CT = metrics[1]
		CFT = metrics[2]
		alpha = metrics[3]
		MCF = metrics[4]
		LI30 = calculateLysisAtTime(A,t,30.0)
		sel_metrics = [CT,CFT,alpha,MCF,LI30]
		#@show broadcast(abs, sel_metrics-sel_target)
	end
	xlabel("Time (minutes)", fontsize = 20)
	ylabel("Clot Amplitude (mm)", fontsize = 20)
	titlestr =string("Patient ", id, " Hour ", hour)
	title(titlestr)
	#let's plot ideal metrics on top, in red
	plot(sel_target[1]/60, 2, "rx", markersize = 20) #CT
	plot(sel_target[1]/60+sel_target[2]/60, 20, "rx", markersize = 20) #CFT
	plot([sel_target[1]/60, sel_target[1]/60+sel_target[2]/60], [0,20], color = "red", linestyle = "dotted") #alpha
	plot([20,40], [sel_target[4], sel_target[4]], color = "red", linestyle ="dashed") #MCF
	#@show 1-sel_target[5]/100
	plot(30, sel_target[4]*(1-sel_target[5]/100), "rx", markersize=20) #LI30

	savefig(string("../figures/UW_ParamEst/",titlestr,".pdf"))

end



function updateICsForPatient(id::Int, hour::Int, all_data, initial_parameter_estimate)
	curr_dat = all_data[.&(all_data.id .== id, all_data.hour.==hour), :]
	nominal_ICs = initial_parameter_estimate
	offset =1
	names = ["FII", "PC", "ATIII", "Fibrinogen", "plasminogen"]

	#all in nM
	nominal_II = nominal_ICs[1]
	nominal_PC =nominal_ICs[3]
	nominal_ATIII = nominal_ICs[5]
	nominal_Fibrinogen = nominal_ICs[14]
	nominal_plasminogen = nominal_ICs[15]

	actual_FII = nominal_II*curr_dat.II[1]/100
	actual_PC = nominal_PC*curr_dat.ProtC[1]/100
	actual_ATIII = nominal_ATIII*curr_dat.AT[1]/100
	actual_Fibrinogen = convertFibrinogenTonM(curr_dat.Fbgn[1])
	actual_plasminogen = nominal_plasminogen*curr_dat.Plgn[1]/100

	#update intiial condition vector
	initial_parameter_estimate[1]=actual_FII
	initial_parameter_estimate[3]=actual_PC
	initial_parameter_estimate[5]=actual_ATIII
	initial_parameter_estimate[14]=actual_Fibrinogen
	initial_parameter_estimate[15]=actual_plasminogen
	return initial_parameter_estimate
	
end

function constructTarget(hour::Int, patientID::Int, all_data)
	curr_dat = all_data[.&(all_data.id .== patientID, all_data.hour.==hour), :]
	target_CT =  curr_dat.teg_rapidteg_r[1]*60
	target_CFT = curr_dat.teg_rapidteg_[1]*60 #k is missing because of find and replace used to clean up data
	target_alpha = curr_dat.teg_rapidteg_angle[1]
	target_MCF = curr_dat.teg_rapidteg_ma[1]
	target_LI30 = curr_dat.teg_rapidteg_ly30[1]
	
	target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	return target
end

function objectiveUWHour0(parameter_array)
	obj_array =ones(6,1)
	ids = [855, 856,880]
	TSTART = 0.0
	TSTOP =60.0
	#since we don't actually know platelets, give a normal platelet count'
	curr_platelets =300
	#use the median as characteritic scaling
	scale_CT = 6.7*60
	scale_CFT = 2*60
	scale_alpha = 62.3
	scale_MCF = 60.6
	scale_MaxLysis = .12
	scale_AUC = 250.0
	hours = [0,36]

	sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]
	#weight all metrics equally for now
	weights = ones(size(sel_scales))
	#heavily weight CT and CFT
	weights[1]=1000.0
	weights[2]=1000.0
	count = 1
	for (idx,id) in enumerate(ids)
		for hour in hours
			#@show hour, id, count
			temp_params = parameter_array
			sel_target  = constructTarget(hour, id, all_data)
			dict = buildCompleteDictFromOneVector(temp_params)
			initial_parameter_estimate = dict["INITIAL_CONDITION_VECTOR"]
			initial_condition_vector = updateICsForPatient(id, hour,all_data, initial_parameter_estimate)
			#now use our actually measured ICs to update initial_parameter_estimate
			fbalances(y,p,t)= Balances(t,y,dict) 
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = .1,maxiters = 1e6)
			t =sol.t
			X = sol
			#A = convertToROTEM(t,X,tPA)
			tPA = initial_condition_vector[16]
			A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
			#plot(t,A)
			metrics = calculateCommonMetrics(A,t)
			CT = metrics[1]
			CFT = metrics[2]
			alpha = metrics[3]
			MCF = metrics[4]
			LI30 = calculateLysisAtTime(A,t,30.0)
			sel_metrics = [CT,CFT,alpha,MCF,LI30]
			#@show sel_metrics
			#if any of the metrics are negative, something went wrong, penalize our parameters
			if(true in (metrics .<0)) #if any of our metrics are negative
				calc_obj = 1E8
			else
				sum = 0.0
				#use MSE to determine how far away we are from hitting our target
				for j in 1:maximum(size(sel_metrics))
					sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2*weights[j]
				end
				#@show params
				calc_obj = sqrt(sum)
			end
			obj_array[count]=calc_obj
			count = count+1
		end
	end
	@show obj_array
	return obj_array

end

function attemptOptimization_UW()
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 6
	outputfile = "../parameterEstimation/POETS_info_18_09_19_PlateletContributionToROTEM_UW_Scaled.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	initial_parameter_estimate=vec(mean(allp, dims=1))

	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*4, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/4, initial_parameter_estimate[47:end]/1000)
	EC = []
	PC = []
	RA = []
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		if(index==1)
			initial_parameter_estimate =initial_parameter_estimate+initial_parameter_estimate*rand()*.1
		else
			#use results from previous round to start
			sum_errs = sum(EC, dims=1)
			min_idx =  findmin(sum_errs)[2][2] #weird indexing because it returns a cartesian idx
			bestParams = PC[:, min_idx]
			initial_parameter_estimate = bestParams

		end

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveUWHour0,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)

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
