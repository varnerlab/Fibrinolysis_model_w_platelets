using POETs

include("runModel.jl")
include("inverseProbUtils.jl")

#load data once
pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
all_data = CSV.read(pathToData)

function updateICsForPatient(id, all_data, initial_parameter_estimate)
	hour = 0
	curr_dat = all_data[.&(all_data.id .== id, all_data.hour.==hour), :]
	nominal_ICs = initial_parameter_estimate
	offset =1
	names = ["FII", "PC", "ATIII", "Fibrinogen", "plasminogen"]
	est_FII = estimate[offset]
	est_PC = estimate[offset+2]
	est_ATIII= estimate[offset+4]
	est_Fibrinogen = estimate[offset+13]
	est_Plasminogen = estimate[offset+14]
	all_est = [est_FII, est_PC, est_ATIII, est_Fibrinogen, est_Plasminogen]

	#all in nM
	nominal_II = nominal_ICs[1]
	nominal_PC =nominal_ICs[3]
	nominal_ATIII = nominal_ICs[5]
	nominal_Fibrinogen = nominal_ICs[14]
	nominal_plasminogen = nominal_ICs[15]

	actual_FII = nominal_II*rel_dat.II[1]/100
	actual_PC = nominal_PC*rel_dat.ProtC[1]/100
	actual_ATIII = nominal_ATIII*rel_dat.AT[1]/100
	actual_Fibrinogen = convertFibrinogenTonM(rel_dat.Fbgn[1])
	actual_plasminogen = nominal_plasminogen*rel_dat.Plgn[1]/100

	#update intiial condition vector
	initial_parameter_estimate[1]=actual_FII
	initial_parameter_estimate[3]=actual_PC
	initial_parameter_estimate[5]=actual_ATIII
	initial_parameter_estimate[14]=actual_Fibrinogen
	initial_parameter_estimate[15]=actual_plasminogen
	return initial_parameter_estimate
	
end

function objectiveUWHour0(parameter_array)
	obj_array =SharedArray{Float64}(3,1)
	ids = [855, 856,880]
	for (idx,id) in ids
		temp_params = parameter_array
		dict = buildCompleteDictFromOneVector(temp_params)
		initial_parameter_estimate = dict["INITIAL_CONDITION_VECTOR"]
		initial_parameter_estimate = updateICsForPatient(id, all_data, initial_parameter_estimate)
		#now use our actually measured ICs to update initial_parameter_estimate
		fbalances(y,p,t)= Balances(t,y,dict) 
		prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
		sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = .1,maxiters = 1e6)
		t =sol.t
		X = sol
		#A = convertToROTEM(t,X,tPA)
		A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
		metrics = calculateCommonMetrics(A,t)
		CT = metrics[1]
		CFT = metrics[2]
		alpha = metrics[3]
		MCF = metrics[4]
		LI30 = calculateLysisAtTime(R,T,30.0)
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
		obj_array[idx]=calc_obj
	
	end
	return calc_obj

end

function attemptOptimization_UW()
	number_of_subdivisions = 10
	number_of_parameters = 77
	number_of_objectives = 3
	outputfile = "../parameterEstimation/POETS_info_13_09_19_PlateletContributionToROTEM_UW_Scaled.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	initial_parameter_estimate=mean(allp, dims=1)

	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*2, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/2, initial_parameter_estimate[47:end]/1000)
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		initial_parameter_estimate =initial_parameter_estimate+initial_parameter_estimate*rand()*.1

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
