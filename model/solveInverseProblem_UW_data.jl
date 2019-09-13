#to repeat solving the same inverse problem. To test how well our optimization technique and scaling do on actual data
using CSV
using Statistics
using Random
using NLopt
using Optim #for PSO
using Distributed

include("runModel.jl")
include("inverseProbUtils.jl")

function getPatientData(patient_number::Int, hour::Int)
	pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
	all_data = CSV.read(pathToData)
	rel_dat =all_data[.&(all_data.id .== patient_number, all_data.hour.==hour), :]
	return rel_dat
end


function compareEstimateToReal(nominal_ICs, estimate, rel_dat, output_fn)
	offset =9
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
	all_actual = [actual_FII, actual_PC, actual_ATIII, actual_Fibrinogen, actual_plasminogen]
	for j in 1:maximum(size(all_est))
		println(names[j], ": ", "Actual: ", all_actual[j], " Estimated: ", all_est[j])
	end
	f =open(output_fn, "w")
	write(f, "Name \t Actual \t Estimated \n")
	for j in 1:maximum(size(all_est))
		write(f, string(names[j], "\t", all_actual[j], "\t", all_est[j], "\n"))
	end
	close(f)
	
end

function compareMetricsToReal(foundICs, sel_metrics, output_fn_metrics)
	curr_exp = foundICs[1:8]
	curr_ICs = foundICs[9:end-1]
	curr_platelets = foundICs[end]
	T,R=runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	LI30 = calculateLysisAtTime(R,T,30.0)
	est_metrics = [CT,CFT,alpha,MCF,LI30]
	names = ["R", "K", "alpha", "MA", "LI30"]
	println("Metrics")
	for j in 1:maximum(size(est_metrics))
		println(names[j], ": ", "Actual: ", sel_metrics[j], " Estimated: ", est_metrics[j])
	end
	f =open(output_fn_metrics, "w")
	write(f, "Name \t Actual \t Estimated \n")
	for j in 1:maximum(size(est_metrics))
		write(f, string(names[j], "\t", sel_metrics[j], "\t", est_metrics[j], "\n"))
	end
	close(f)
end

function runAllPatientsAndHours()
	pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
	all_data = CSV.read(pathToData)
	for j in 1:size(all_data, 1)
		curr_id = all_data.id[j]
		curr_hr = all_data.hour[j]
		poseInverseProblem(curr_id, curr_hr)
	end
end

function runAllPatientsAndHoursPSO()
	pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
	all_data = CSV.read(pathToData)
	for j in 1:size(all_data, 1)
		curr_id = all_data.id[j]
		curr_hr = all_data.hour[j]
		poseInverseProblemPSO(curr_id, curr_hr)
	end
end

function runAllPatientsAndHoursThreeStages()
	pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
	all_data = CSV.read(pathToData)
	for j in 1:size(all_data, 1)
		curr_id = all_data.id[j]
		curr_hr = all_data.hour[j]
		poseInverseProblemThreeStages(curr_id, curr_hr)
	end
end

function poseInverseProblem(patient_number::Int, hour::Int)
	rel_dat = getPatientData(patient_number, hour)
	# Load data -
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	global kin_params=mean(allp, dims=1)
	@show size(kin_params)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	#add some tPA
	nominal_ICs[16]=.073
	#add tissue factor to 35 pM
	nominal_ICs[7]=.035
	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =2.01*all_nominal
	trigger_idx = size(nominal_experimental,1)+6
	ups[trigger_idx]=all_nominal[trigger_idx]*1.1
	#make it possible for things that are zero to move some
	ups[ups.==0]=fill(2.0, size(ups[ups.==0]))
	#replace zeros in nominal conditions with eps, otherwise NLopt gets upset
	all_nominal[all_nominal.==0]=fill(eps(), size(all_nominal[all_nominal.==0]))

	#create Optization problem using DIRECT
	opt = NLopt.Opt(:GN_DIRECT, size(all_nominal,1))
	#set bounds
	lower_bounds!(opt, vec(lbs))
	upper_bounds!(opt, vec(ups))

	#set objective
	min_objective!(opt, objective_UW)
	ftol_rel!(opt, 1E-6)

	#create our objective 
	#times are in minutes, so need to convert into seconds
	target_CT =  rel_dat.teg_rapidteg_r[1]*60
	target_CFT = rel_dat.teg_rapidteg_[1]*60 #k is missing because of find and replace used to clean up data
	target_alpha = rel_dat.teg_rapidteg_angle[1]
	target_MCF = rel_dat.teg_rapidteg_ma[1]
	target_LI30 = rel_dat.teg_rapidteg_ly30[1]
	#don't have AUC
	#use the median as characteritic scaling
	scale_CT = 6.7*60
	scale_CFT = 2*60
	scale_alpha = 62.3
	scale_MCF = 60.6
	scale_MaxLysis = .12
	scale_AUC = 250.0

	global sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	@show sel_target
	global weights = [1,1,1.0,1,1000.0,1] 
	@show opt
	#run optimization
	
	#use nominal ICs as initial guess
	(minf,minx,ret) = NLopt.optimize(opt, vec(vcat(nominal_experimental,nominal_ICs, platelets)))
	
	println("got $minf at $minx after $count iterations (returned $ret)")

	println("Staring Local Optimization")
	optlocal = NLopt.Opt(:LN_NELDERMEAD, size(all_nominal,1))
	#set bounds
	lower_bounds!(optlocal, vec(minx*.5))
	upper_bounds!(optlocal, vec(minx*1.5))
	min_objective!(optlocal,objective_UW)
	maxeval!(optlocal, 1000)
	
	@time (minf_local,minx_local,ret_local) = NLopt.optimize(opt, vec(minx))
	
	println("got $minf_local at $minx_local after $count iterations (returned $ret_local)")
	writedlm(string("../solveInverseProb/UW/foundIcs_11_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), minx_local)
	#store and print our estimates
	output_fn = string("../solveInverseProb/UW/Acc_11_09_19_comparisonPatient", patient_number, "Hour", hour, ".txt")
	compareEstimateToReal(nominal_ICs, minx_local, rel_dat, output_fn)
	output_fn_metrics = string("../solveInverseProb/UW/Metrics_11_09_19_comparisonPatient", patient_number, "Hour", hour, ".txt")
	compareMetricsToReal(minx_local, sel_target, output_fn_metrics)	
end

function poseInverseProblemThreeStages(patient_number::Int, hour::Int)
	rel_dat = getPatientData(patient_number, hour)
	# Load data -
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	global kin_params=mean(allp, dims=1)
	@show size(kin_params)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	#add some tPA
	nominal_ICs[16]=.073
	#add tissue factor to 35 pM
	nominal_ICs[7]=.035
	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =3.01*all_nominal

	#make it possible for things that are zero to move some
	ups[ups.==0]=fill(2.0, size(ups[ups.==0]))
	#Let's be a bit smarter with bounds
#	trigger_idx = size(nominal_experimental,1)+7
#	ups[trigger_idx]=all_nominal[trigger_idx]*1.1
#	fIIa_idx = size(nominal_experimental,1)+2
#	ups[fIIa_idx]=.2
#	APC_idx = size(nominal_experimental,1)+4
#	ups[APC_idx]=.2
	#replace zeros in nominal conditions with eps, otherwise NLopt gets upset
	all_nominal[all_nominal.==0]=fill(eps(), size(all_nominal[all_nominal.==0]))

	#create Optization problem using DIRECT
	opt = Opt(:GN_CRS2_LM, size(all_nominal,1))
	#set bounds
	lower_bounds!(opt, vec(lbs))
	upper_bounds!(opt, vec(ups))

	#set objective
	min_objective!(opt, objective_UW)
	#give 5 minutes to solve
	#maxtime!(opt, 1*60)
	ftol_rel!(opt, 1E-2)

	#create our objective 
	#times are in minutes, so need to convert into seconds
	target_CT =  rel_dat.teg_rapidteg_r[1]*60
	target_CFT = rel_dat.teg_rapidteg_[1]*60 #k is missing because of find and replace used to clean up data
	target_alpha = rel_dat.teg_rapidteg_angle[1]
	target_MCF = rel_dat.teg_rapidteg_ma[1]
	target_LI30 = rel_dat.teg_rapidteg_ly30[1]
	#don't have AUC
	#use the median as characteritic scaling
	scale_CT = 6.7*60
	scale_CFT = 2*60
	scale_alpha = 62.3
	scale_MCF = 60.6
	scale_MaxLysis = .12
	scale_AUC = 250.0

	global sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	@show sel_target
	global weights = [50.0,1,1.0,100.0,1] 
	@show opt
	#run optimization
	
	#use nominal ICs as initial guess
	(minf,minx,ret) = NLopt.optimize(opt, vec(vcat(nominal_experimental,nominal_ICs, platelets)))
	
	println("got $minf at $minx after $count iterations (returned $ret)")

	println("Staring Stage 2 Optimization")
	optlocal = Opt(:LN_NELDERMEAD, size(all_nominal,1))
	#set bounds
	lower_bounds!(optlocal, .25*minx)
	upper_bounds!(optlocal, 4*minx)
	min_objective!(optlocal,objective_UW)
	#maxeval!(optlocal, 500)
	#give 5 minutes to solve
	#maxtime!(optlocal, 1*60)
	ftol_rel!(optlocal, 1E-2)
	
	@time (minf_local,minx_local,ret_local) = NLopt.optimize(opt, vec(minx))
	println("got $minf_local at $minx_local after $count iterations (returned $ret_local)")

	println("Starting Stage 3 Optimization")
	optStage3 = Opt(:LN_COBYLA, size(all_nominal,1))
	#set bounds
	lower_bounds!(optStage3, .5*minx)
	upper_bounds!(optStage3, 2*minx)
	min_objective!(optStage3,objective_UW)
	#maxtime!(optStage3, 1*60)
	ftol_rel!(optStage3, 1E-2)
	
	@time (minf_Stage3,minx_Stage3,ret_Stage3) = NLopt.optimize(opt, vec(minx_local))
	
	
	println("got $minf_Stage3 at $minx_Stage3 after $count iterations (returned $ret_Stage3)")
	writedlm(string("../solveInverseProb/UW/ThreeStagefoundIcs_09_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), minx_Stage3)
	#store and print our estimates
	output_fn = string("../solveInverseProb/UW/ThreeStageAcc_09_09_19_comparisonPatient", patient_number, "Hour", hour, ".txt")
	output_fn_metrics = string("../solveInverseProb/UW/ThreeStageMetrics_09_09_19_comparisonPatient", patient_number, "Hour", hour, ".txt")
	compareEstimateToReal(nominal_ICs, minx_Stage3, rel_dat, output_fn)
	compareMetricsToReal(minx_Stage3, sel_target, output_fn_metrics)	
end


function poseInverseProblemPSO(patient_number::Int, hour::Int)
	rel_dat = getPatientData(patient_number, hour)
	# Load data -
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	global kin_params=mean(allp, dims=1)
	@show size(kin_params)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	#add some tPA
	nominal_ICs[16]=.073
	#add tissue factor to 35 pM
	nominal_ICs[7]=.035
	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =2.01*all_nominal
	#make it possible for things that are zero to move some
	ups[ups.==0]=fill(1.0, size(ups[ups.==0]))
	#replace zeros in nominal conditions with eps, otherwise NLopt gets upset
	all_nominal[all_nominal.==0]=fill(eps(), size(all_nominal[all_nominal.==0]))

	#create our objective 
	#times are in minutes, so need to convert into seconds
	target_CT =  rel_dat.teg_rapidteg_r[1]*60
	target_CFT = rel_dat.teg_rapidteg_[1]*60 #k is missing because of find and replace used to clean up data
	target_alpha = rel_dat.teg_rapidteg_angle[1]
	target_MCF = rel_dat.teg_rapidteg_ma[1]
	target_LI30 = rel_dat.teg_rapidteg_ly30[1]
	#don't have AUC
	#use the median as characteritic scaling
	scale_CT = 6.7*60
	scale_CFT = 2*60
	scale_alpha = 62.3
	scale_MCF = 60.6
	scale_MaxLysis = .12
	scale_AUC = 250.0

	#global sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]
	global sel_scales = ones(1,5)

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	@show sel_target
	global weights = [1,1,1.0,100.0,1] 
	#run optimization
	
	#use nominal ICs as initial guess
	initial_guess = vec(vcat(nominal_experimental,nominal_ICs, platelets))
	results= Optim.optimize(objective_UW, initial_guess, ParticleSwarm(lower = lbs, upper = ups, n_particles = 100), Optim.Options( show_trace = true, show_every = 1, time_limit = 15*60))	
	@show results
	writedlm(string("../solveInverseProb/UW/foundIcs_13_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), results.minimizer)
	writedlm(string("../solveInverseProb/UW/CompletePSOResults_13_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), replace(string(results), "\n"=> ""))
	#store and print our estimates
	output_fn = string("../solveInverseProb/UW/Acc_13_09_19_comparisonPSOPatient", patient_number, "Hour", hour, ".txt")
	compareEstimateToReal(nominal_ICs, results.minimizer, rel_dat, output_fn)
	output_fn_metrics = string("../solveInverseProb/UW/PSOeMetrics_13_09_19_comparisonPatient", patient_number, "Hour", hour, ".txt")
	compareMetricsToReal(results.minimizer, sel_target, output_fn_metrics)		
end

function poseInverseProblemPSOIterative(patient_number::Int, hour::Int)
	#do we do better if we try to fit things one item at a time?
	rel_dat = getPatientData(patient_number, hour)
	# Load data -
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	global kin_params=mean(allp, dims=1)
	@show size(kin_params)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	#add some tPA
	nominal_ICs[16]=.073
	#add tissue factor to 35 pM
	nominal_ICs[7]=.035

	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =4.01*all_nominal
	#make it possible for things that are zero to move some
	ups[ups.==0]=fill(5.0, size(ups[ups.==0]))
	#replace zeros in nominal conditions with eps, otherwise NLopt gets upset
	all_nominal[all_nominal.==0]=fill(eps(), size(all_nominal[all_nominal.==0]))

	#create our objective 
	#times are in minutes, so need to convert into seconds
	target_CT =  rel_dat.teg_rapidteg_r[1]*60
	target_CFT = rel_dat.teg_rapidteg_[1]*60 #k is missing because of find and replace used to clean up data
	target_alpha = rel_dat.teg_rapidteg_angle[1]
	target_MCF = rel_dat.teg_rapidteg_ma[1]
	target_LI30 = rel_dat.teg_rapidteg_ly30[1]
	#don't have AUC
	#use the median as characteritic scaling
	scale_CT = 6.7*60
	scale_CFT = 2*60
	scale_alpha = 62.3
	scale_MCF = 60.6
	scale_MaxLysis = .12
	scale_AUC = 250.0

	global sel_scales = [scale_MCF, scale_CFT, scale_alpha, scale_CT, scale_MaxLysis]

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_MCF,target_CFT, target_CT, target_alpha,target_LI30]
	@show sel_target
	global weights = [1,1,1.0,1,1,1] 
	#run optimization
	
	#use nominal ICs as initial guess
	initial_guess = vec(vcat(nominal_experimental,nominal_ICs, platelets))
	global metric_number = 1
	results = []
	output_fn = string("../solveInverseProb/UW/Acc_13_09_19_comparisonPSOPatient", patient_number, "Hour", hour, ".txt")
	output_fn_metrics = string("../solveInverseProb/UW/Metrics_13_09_19_comparisonPSOPatient", patient_number, "Hour", hour, ".txt")
	for j=1:maximum(size(sel_target))
		results= Optim.optimize(objective_UW_iterative, initial_guess, ParticleSwarm(lower = lbs, upper = ups, n_particles = 100), Optim.Options( show_trace = true, show_every = 1,iterations=30))
		initial_guess = results.minimizer
		metric_number = metric_number+1
		#update weights so we really care about what we've already optimized
		weights[j]=weights[j]*1000
		@show weights
		compareEstimateToReal(nominal_ICs, results.minimizer, rel_dat, output_fn)
		compareMetricsToReal(results.minimizer, sel_target, output_fn_metrics)
	end	
	@show results
	writedlm(string("../solveInverseProb/UW/foundIcs_13_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), results.minimizer)
	writedlm(string("../solveInverseProb/UW/CompletePSOResults_13_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), replace(string(results), "\n"=> ""))
	#store and print our estimates
	compareEstimateToReal(nominal_ICs, results.minimizer, rel_dat, output_fn)
	compareMetricsToReal(results.minimizer, sel_target, output_fn_metrics)	
end
