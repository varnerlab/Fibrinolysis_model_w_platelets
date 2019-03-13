#to solve the inverse problem-we know ROTEM metrics, but can we back out ICs?
#using NLopt
using Random #to set seed

include("runModel.jl")
include("inverseProbUtils.jl")


function testIfPhysical(params,genIC,genExp,genPlatelets)
	T,R=runModelWithParamsChangeICReturnA(params,genIC,genExp,genPlatelets)

	#describe the curve we've created'
	metrics = calculateCommonMetrics(R,T)
	target_CT = metrics[1]
	target_CFT = metrics[2]
	target_alpha = metrics[3]
	target_MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	target_MaximumLysis = calculateLysisAtTime(R,T,60.0)
	target_LI30 = calculateLysisAtTime(R,T,30.0)
	target_AUC = calculateAUC(R,T)


	temp_target=[target_CT, target_CFT, target_alpha, target_MCF,target_MaximumLysis, target_AUC, target_LI30]
	#@show temp_target
	if(true in (temp_target .<0))
		isPhysical = false
	else
		isPhysical = true
	end
	#based on ""Normal range values for thromboelastography in healthy adult volunteers""
	#want MCF between 50 and 70
	#MA = MCF
	#R = CT
	#K = CFT
	if(target_MCF<49.7 || target_MCF>72.7)
		isPhysical=false
	end
	if(target_CT<3.8*60  || target_CT>9.8*60)
		isPhysical=false
	end
	if(target_CFT<.7*60  || target_CFT>3.4*60)
		isPhysical=false
	end
#	if(target_alpha<47.8  || target_CFT>77.7)
#		isPhysical=false
#	end

	return isPhysical
end

function runNLopt(seed,iter)
	# Load data -
	global kin_params = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	#nominal_ICs=[1106.0, 0.0, 60.0, 0.0, 2686.0, 0.79, 3.95, 0.0, 0.0, 7900.0, 971.7, 1.58, 0.0, 0.0, 0.0, 932.2, 0.4424, 0.0]
	#nominal_experimental=[1.975,15.8,0.553,71.1,134.3,73.47,70.0]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]

	@show size(nominal_experimental), size(nominal_ICs), size(platelets)

	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =4*all_nominal
	#make it possible for things that are zero to move some
	#ups[ups.==0]=2.0
	ups[ups.==0]=fill(2.0, size(ups[ups.==0]))
	#ups[9]=.1 #limit the amount of FIIa we can start with
	#@show all_nominal[9]
	#bound that we must have some prothrombin
	#lbs[8]=100.0

	#replace zeros in nominal conditions with eps, otherwise NLopt gets upset
	all_nominal[all_nominal.==0]=fill(eps(), size(all_nominal[all_nominal.==0]))

	#create Optization problem using DIRECT
	opt = Opt(:GN_DIRECT, size(all_nominal,1))
	#set bounds
	lower_bounds!(opt, vec(lbs))
	upper_bounds!(opt, vec(ups))

	@show lbs.<ups
	#set objective
	min_objective!(opt, objective_six_metrics_weighted)
	ftol_rel!(opt, 1E-4)

	#normal ranges for NATEM from https://www.sciencedirect.com/science/article/pii/S0049384810004081?via%3Dihub
	#CT-461-917 s
	#CFT 105-409 S
	#MCF 26-60 mm
	#alpha 29-65 degrees
	#maximumLysis 8-24 percent
	#maximumVelocity 2-10 mm/min

	#set seed so this is repeatable
	#srand(seed)
	Random.seed!(seed)
	#create our data
	#conditions taken from SampleSingleParameterSet adjustments to IC and experimental conditions
	temp_IC = d["INITIAL_CONDITION_VECTOR"]
	temp_exp = d["FACTOR_LEVEL_VECTOR"]
	temp_platelets = d["PLATELET_PARAMS"][5] 
	stretchfactor = 4.0 #for setting upper and lower bounds of our generated ROTEM curve
	genIC =temp_IC/stretchfactor+rand(size(temp_IC,1)).*(stretchfactor*temp_IC-temp_IC*1/stretchfactor)
	genExp =temp_exp/stretchfactor+rand(size(temp_exp,1)).*(stretchfactor*temp_exp-temp_exp*1/stretchfactor)
#	@show temp_platelets/stretchfactor, rand(size(temp_platelets,1))[1]
	genPlatelets =temp_platelets/stretchfactor+rand(size(temp_platelets,1))[1].*(stretchfactor*temp_platelets-temp_platelets*1/stretchfactor)

	genIC = vec(genIC)
	genExp = vec(genExp)
	tPA = genIC[12]

	isPhysical = testIfPhysical(kin_params,genIC,genExp,genPlatelets)
	while(!isPhysical)
		#this parameters produce a nonsensical results, regenerate and repeat as necceassary 
		genIC =temp_IC/stretchfactor+rand(size(temp_IC,1)).*(stretchfactor*temp_IC-temp_IC*1/stretchfactor)
		genExp =temp_exp/stretchfactor+rand(size(temp_exp,1)).*(stretchfactor*temp_exp-temp_exp*1/stretchfactor)
		genPlatelets =temp_platelets/stretchfactor+rand(size(temp_platelets,1))[1].*(stretchfactor*temp_platelets-temp_platelets*1/stretchfactor)
		genIC = vec(genIC)
		genExp = vec(genExp)
		#let time delay range between 400 and 1500
		genMX = 400+rand()*(1500-400)
		tPA = genIC[12]
		@show genIC, genExp, genPlatelets, genMX
		isPhysical = testIfPhysical(kin_params,genIC,genExp,genPlatelets)
	end

	#Store our ICS to disk
	writedlm(string("../solveInverseProb/ics_to_match_11_03_19_iter",iter ,".txt"), hcat(genExp', genIC', genPlatelets), ',')

	#generate the curve we're fitting
	R,T =runModelWithParamsChangeICReturnA(kin_params,genIC,genExp,genPlatelets)

	#describe the curve we've created'
	metrics = calculateCommonMetrics(R,T)
	target_CT = metrics[1]
	target_CFT = metrics[2]
	target_alpha = metrics[3]
	target_MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	target_MaximumLysis = calculateLysisAtTime(R,T,120.0)
	target_AUC = calculateAUC(R,T)

	#use the midpoint as characteritic scaling
	scale_CT = 500.0
	scale_CFT = 250.0
	scale_alpha = 40.0
	scale_MCF = 40.0
	scale_MaxLysis = .84
	scale_AUC = 2500.0

	global sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis, scale_AUC]

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_MaximumLysis, target_AUC]
	@show sel_target
	global weights = [1,1,1.0,1,1,1] 
	@show opt
	#run optimization
	
	@show (all_nominal .< ups)
	@show (all_nominal .> lbs)
	@time (minf,minx,ret) = optimize(opt, vec(all_nominal))
	
	println("got $minf at $minx after $count iterations (returned $ret)")

	println("Staring Local Optimization")
	optlocal = opt = Opt(:LN_NELDERMEAD, size(all_nominal,1))
	#set bounds
	lower_bounds!(optlocal, vec(minx*.5))
	upper_bounds!(optlocal, vec(minx*1.5))
	min_objective!(optlocal, objective_six_metrics_weighted)
	
	@time (minf_local,minx_local,ret_local) = optimize(opt, vec(minx))
	
	println("got $minf_local at $minx_local after $count iterations (returned $ret_local)")
	writedlm(string("../solveInverseProb/foundIcs_11_03_19_", iter, ".txt"), minx_local)
end

#for j = 1:10
#	runNLopt(j,j)
#end






