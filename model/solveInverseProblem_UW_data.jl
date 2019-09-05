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

function convertFibrinogenTonM(fibrinogen_in_mg_dl::Real)
	mw_fibrinogen = 340.0 #kDa
	kilodaltons_in_a_mg =  6.0221366516752E17
	deci_liters_per_liter =10.0
	NA = 6.0221409e23 #avagadro's constant'
	conv_fibrinogen = fibrinogen_in_mg_dl*1/mw_fibrinogen*kilodaltons_in_a_mg*deci_liters_per_liter*1/NA*10^9
	return conv_fibrinogen
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

function runAllPatientsAndHours()
	pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
	all_data = CSV.read(pathToData)
	for j in 1:size(all_data, 1)
		curr_id = all_data.id[j]
		curr_hr = all_data.hour[j]
		poseInverseProblem(curr_id, curr_hr)
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
	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =2.01*all_nominal
	#make it possible for things that are zero to move some
	ups[ups.==0]=fill(2.0, size(ups[ups.==0]))
	#replace zeros in nominal conditions with eps, otherwise NLopt gets upset
	all_nominal[all_nominal.==0]=fill(eps(), size(all_nominal[all_nominal.==0]))

	#create Optization problem using DIRECT
	opt = Opt(:GN_DIRECT, size(all_nominal,1))
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
	global weights = [1,1,1.0,1,1,1] 
	@show opt
	#run optimization
	
	#use nominal ICs as initial guess
	(minf,minx,ret) = optimize(opt, vec(vcat(nominal_experimental,nominal_ICs, platelets)))
	
	println("got $minf at $minx after $count iterations (returned $ret)")

	println("Staring Local Optimization")
	optlocal = Opt(:LN_NELDERMEAD, size(all_nominal,1))
	#set bounds
	lower_bounds!(optlocal, vec(minx*.5))
	upper_bounds!(optlocal, vec(minx*1.5))
	min_objective!(optlocal,objective_UW)
	maxeval!(optlocal, 500)
	
	@time (minf_local,minx_local,ret_local) = optimize(opt, vec(minx))
	
	println("got $minf_local at $minx_local after $count iterations (returned $ret_local)")
	writedlm(string("../solveInverseProb/UW/foundIcs_04_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), minx_local)
	#store and print our estimates
	output_fn = string("../solveInverseProb/UW/Acc_comparisonPatient", patient_number, "Hour", hour, ".txt")
	compareEstimateToReal(nominal_ICs, minx_local, rel_dat, output_fn)	
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
	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =2.01*all_nominal
	#make it possible for things that are zero to move some
	ups[ups.==0]=fill(2.0, size(ups[ups.==0]))
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

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	@show sel_target
	global weights = [1,1,1.0,1,1,1] 
	#run optimization
	
	#use nominal ICs as initial guess
	initial_guess = vec(vcat(nominal_experimental,nominal_ICs, platelets))
	results= Optim.optimize(objective_UW, initial_guess, ParticleSwarm(lower = lbs, upper = ups, n_particles = 40), Optim.Options( show_trace = true, show_every = 1, time_limit = 15*60))	
	@show results
	writedlm(string("../solveInverseProb/UW/foundIcs_05_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), results.minimizer)
	writedlm(string("../solveInverseProb/UW/CompletePSOResults_05_09_19_", "solvingPatient", patient_number, "Hour", hour, ".txt"), string(results))
	#store and print our estimates
	output_fn = string("../solveInverseProb/UW/Acc_comparisonPSOPatient", patient_number, "Hour", hour, ".txt")
	compareEstimateToReal(nominal_ICs, results.minimizer, rel_dat, output_fn)	
end
