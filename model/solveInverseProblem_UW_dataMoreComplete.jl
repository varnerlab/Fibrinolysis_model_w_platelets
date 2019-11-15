#to repeat solving the same inverse problem. To test how well our optimization technique and scaling do on actual data
using CSV
using Statistics
using Random
using NLopt
using Optim #for PSO
using Distributed
using Interpolations

include("runModel.jl")
include("inverseProbUtils.jl")

function getPatientData(patient_number::Int, hour::Int)
	pathToData ="/home/rachel/Documents/washington_burns/PatientsWTEGAndMoreFactors.csv"
	all_data = CSV.read(pathToData)
	rel_dat =all_data[.&(all_data.id .== patient_number, all_data.hour.==hour), :]
	return rel_dat
end

function runAllPatientsAndHoursPSO()
	pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndMoreFactors.csv"
	all_data = CSV.read(pathToData)
	for j in 1:size(all_data, 1)
		curr_id = all_data.id[j]
		curr_hr = all_data.hour[j]
		poseInverseProblemPSO_MoreCompleteData(curr_id, curr_hr)
		#poseInverseProblemPSO_MoreCompleteData_UseContructedCurves(curr_id, curr_hr)
	end
end

function constructCurveFromMetrics(metrics)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	LI30 = metrics[5]
	est_trajectory = zeros(60*60) #use time scale of one point per second
	est_t = 0:1:60*60
	#from zero to CT
	start_increase_t = floor(CT/2)
	itp_1 = LinearInterpolation([start_increase_t-1, CT+1], [0, 2.0])
	
	for k in Int(start_increase_t):Int(CT)
		est_trajectory[k]=itp_1(est_t[k])
	end
	itp_2 = LinearInterpolation([CT-1, CT+CFT+1], [2,20.0])
	for k in Int(CT):Int(CT+CFT)
		est_trajectory[k]=itp_2(est_t[k])
	end
	tmax = 12*60 #assume we hit max around 12 minutes
	if(CT+CFT>tmax) #deal with the very slow person
		tmax = CT+CFT+300.00
	end
	@show CT+CFT-1, tmax
	itp_3 = LinearInterpolation([CT+CFT-1, tmax], [20.0, MCF])
	for k in Int(CT+CFT):Int(tmax)
		est_trajectory[k]=itp_3(est_t[k])
	end

	tLysis = 30.02*60
	if(tmax-1>tLysis)
		tLysis = tmax+10*60
	end
	itp_4 = LinearInterpolation([tmax-1, tLysis], [MCF,MCF*(1-LI30/100) ])
	for k in tmax:60*30
		est_trajectory[k]=itp_4(est_t[k])
	end

	est_trajectory[60*30:length(est_trajectory)]= fill(MCF*(1-LI30/100), length(est_trajectory[60*30:length(est_trajectory)]))
	
	#now actually want scale of 10 points per minute
	est_trajectory = est_trajectory[1:6:3600]
	return est_trajectory
end

function compareEstimateToReal(nominal_ICs, estimate, rel_dat, output_fn)
	offset =9
	names = ["FII", "PC", "ATIII", "Fibrinogen", "plasminogen", "FV_FX", "PAI-1", "TFPI", "TAFI"]
	est_FII = estimate[offset]
	est_PC = estimate[offset+2]
	est_ATIII= estimate[offset+4]
	est_Fibrinogen = estimate[offset+13]
	est_Plasminogen = estimate[offset+14]
	est_PAI1 = estimate[offset+21]
	est_TFPI=estimate[1]
	est_TAFI = estimate[7]
	all_est = [est_FII, est_PC, est_ATIII, est_Fibrinogen, est_Plasminogen, est_PAI1, est_TFPI, est_TAFI]

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
	actual_PAI1 = convertPAIToM(rel_dat.PAI1[1])
	actual_TFPI = convertTFPIToM(rel_dat.TFPIFree[1])
	actual_TAFI = convertTAFIToM(rel_dat.TAFIaai[1])

	all_actual = [actual_FII, actual_PC, actual_ATIII, actual_Fibrinogen, actual_plasminogen, actual_PAI1, actual_TFPI, actual_TAFI]
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

function poseInverseProblemPSO_MoreCompleteData(patient_number::Int, hour::Int)
	rel_dat = getPatientData(patient_number, hour)
	# Load data -
	#allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	#allp = readdlm("../parameterEstimation/Best24_UW_Rescaled_16_09_19.txt")
	ec,pc,ra = parsePOETsoutput( "../parameterEstimation/POETS_info_09_10_19_PlateletContributionToROTEM_UW_Scaled_MoreComplete.txt",4)
	numParamSets=3 #the number of paramter sets per objective
	allp=generateNbestPerObjective(numParamSets,ec,pc)
	global kin_params=mean(allp, dims=1)
	@show size(kin_params)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	TFPI = d["FACTOR_LEVEL_VECTOR"][1]
	TAFI = d["FACTOR_LEVEL_VECTOR"][7]
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

	global sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]
	#global sel_scales = ones(1,5)

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	@show sel_target
	global weights = [1,1,1.0,1.0,1] 
	#run optimization
	
	#use nominal ICs as initial guess
	initial_guess = vec(vcat(nominal_experimental,nominal_ICs, platelets))
	results= Optim.optimize(objective_UW_moreComplete, initial_guess, ParticleSwarm(lower = lbs, upper = ups, n_particles = 400), Optim.Options( show_trace = true, show_every = 1, time_limit = 20*60))	
	@show results
	writedlm(string("../solveInverseProb/UW/foundIcs_12_11_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), results.minimizer)
	writedlm(string("../solveInverseProb/UW/CompletePSOResults_12_11_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), replace(string(results), "\n"=> ""))
	#store and print our estimates
	output_fn = string("../solveInverseProb/UW/Acc_12_11_19_comparisonPSOPatient", patient_number, "Hour", hour, ".txt")
	compareEstimateToReal(nominal_ICs, results.minimizer, rel_dat, output_fn)
	output_fn_metrics = string("../solveInverseProb/UW/PSOeMetrics_12_11_19_comparisonPatient", patient_number, "Hour", hour, ".txt")
	compareMetricsToReal(results.minimizer, sel_target, output_fn_metrics)		
end

function poseInverseProblemPSO_MoreCompleteData_UseContructedCurves(patient_number::Int, hour::Int)
	rel_dat = getPatientData(patient_number, hour)
	# Load data -
	#allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	#allp = readdlm("../parameterEstimation/Best24_UW_Rescaled_16_09_19.txt")
	ec,pc,ra = parsePOETsoutput( "../parameterEstimation/POETS_info_09_10_19_PlateletContributionToROTEM_UW_Scaled_MoreComplete.txt",4)
	numParamSets=3 #the number of paramter sets per objective
	allp=generateNbestPerObjective(numParamSets,ec,pc)
	global kin_params=mean(allp, dims=1)
	@show size(kin_params)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	TFPI = d["FACTOR_LEVEL_VECTOR"][1]
	TAFI = d["FACTOR_LEVEL_VECTOR"][7]
	#add some tPA
	nominal_ICs[16]=.073
	#add tissue factor to 35 pM
	nominal_ICs[7]=.035
	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = .5*all_nominal#zeros(size(all_nominal))
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

	global sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]
	#global sel_scales = ones(1,5)

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	@show sel_target
	global weights = [1,1,1.0,1.0,1]
	
	#create a trajectory from these metrics 
	global R_to_match = constructCurveFromMetrics(sel_target)
	#run optimization
	
	#use nominal ICs as initial guess
	initial_guess = vec(vcat(nominal_experimental,nominal_ICs, platelets))
	results= Optim.optimize(objective_MSE_moreComplete, initial_guess, ParticleSwarm(lower = lbs, upper = ups, n_particles = 500), Optim.Options( show_trace = true, show_every = 1, time_limit = 15*60))	
	@show results
	writedlm(string("../solveInverseProb/UW/foundIcs_11_11_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), results.minimizer)
	writedlm(string("../solveInverseProb/UW/CompletePSOResults_11_11_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), replace(string(results), "\n"=> ""))
	#store and print our estimates
	output_fn = string("../solveInverseProb/UW/Acc_11_11_19_comparisonPSOPatient", patient_number, "Hour", hour, ".txt")
	compareEstimateToReal(nominal_ICs, results.minimizer, rel_dat, output_fn)
	output_fn_metrics = string("../solveInverseProb/UW/PSOeMetrics_11_11_19_comparisonPatient", patient_number, "Hour", hour, ".txt")
	compareMetricsToReal(results.minimizer, sel_target, output_fn_metrics)		
end
