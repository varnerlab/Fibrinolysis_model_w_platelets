include("runModel.jl")
using Optim
using Statistics
using Random #so we can set the seed
#to estimate initial conditions for a patient
#let's load the things we're going to want globally

#POETs_data = "../parameterEstimation/POETS_info_05_12_18_PlateletContributionToROTEMFlatness1SmallerConversion.txt"
POETs_data ="../parameterEstimation/POETS_info_02_01_19_PlateletContributionToROTEMFlatness1SmallerConversion.txt"

ec,pc,ra=parsePOETsoutput(POETs_data)
numParamSets = 2
bestparams=generateNbestPerObjective(numParamSets,ec,pc)

function objectiveICs(params)
	tPAs = [0,2]
	totalMSE = 0.0
	weight = 5.0
	for tPA in tPAs
		#run model with current ICs
		alldata, meanROTEM, stdROTEM, TSIM=testROTEMPredicitionGivenParamsPlatetContributionToROTEM(bestparams, patient_id ,tPA, "findingICs.txt", params)
		@show size(meanROTEM)
		#calculate error
		platelets,expdata = setROTEMIC(tPA, patient_id)
		currMSE,interpolatedExperimentalData = calculateMSE(TSIM, meanROTEM, expdata)
		@show currMSE
		if(tPA ==0)
			endMCF = expdata[end,2]
			@show endMCF
			endMCFsim = meanROTEM[end]
			diff = abs(endMCF-endMCFsim)
			@show diff
			totalMSE = totalMSE+diff*weight
		end
		totalMSE = totalMSE+currMSE	
	end
	@show totalMSE, params
	return totalMSE
end

function testROTEMPredicitionGivenParamsPlatetContributionToROTEM(allparams,patient_id,tPA,savestr, currICs)
	numparams = 77
	pathToThrombinData="../data/fromOrfeo_Thrombin_BL_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 100.0
	end
	TSIM = collect(TSTART:Ts:TSTOP)
	platelets,usefuldata = setROTEMIC(tPA, patient_id)
	platelet_count =platelets
	alldata = zeros(1,size(TSIM,1))
	if(size(allparams,1)==numparams) #deal with parameters being stored either vertically or horizontally
		itridx = 2
	else
		itridx = 1
	end
	
	for j in collect(1:size(allparams,itridx))
		if(itridx ==2)
			currparams = vec(allparams[:,j])
		else
			currparams = vec(allparams[j,:])
		end
		#@show currparams
		if(typeof(currparams)==Array{Array,1}) #deal with params being inside an extra layer of array
			currparams= currparams[1]
		end
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = currICs
		initial_condition_vector[16]=tPA #set tPA level
		#initial_condition_vector=setCompleteModelIC(initial_condition_vector,patient_id)
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(y,p,t)= Balances(t,y,dict) 
		#fbalances(t,y)= Balances(t,y,dict) 
		#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
		prob = DifferentialEquations.ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
		@time sol = DifferentialEquations.solve(prob, saveat=.02,maxiters=1E6)
		t =sol.t
		X = sol
		#@show size([a[2] for a in X])
		A = convertToROTEMPlateletContribution(t,X,tPA,platelet_count)
		@show size(A), size(TSIM)
		count = 0
		#deal with failure to solve, but give us a limit number of tries so we don't get stuck
		while(size(A)!=size(TSIM) && count <10)
			prob = DifferentialEquations.ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			@time sol = DifferentialEquations.solve(prob, saveat=.02,maxiters=1E8)
			t =sol.t
			X = sol
			A = convertToROTEMPlateletContribution(t,X,tPA,platelet_count)
			count = count+1
			if(count ==10)
				A = fill(10^6, size(TSIM))
			end
		end

		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,dims=1)
	stdROTEM = std(alldata,dims=1)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function runOptimization(currID)
	rseed =89239
	numFevals =10
	global patient_id = currID #change me as nesseccary
	#kin_params = mean(readdlm("../parameterEstimation/best8_02_05_18.txt"),dims=1)
	kin_params = mean(readdlm("../parameterEstimation/Best4PerObjectiveParameters_01_02_19PlateletContributionToROTEM.txt"), dims=1)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	lbs = nominal_ICs*.5
	ubs = nominal_ICs*1.5
	@show typeof(lbs)
	@show typeof(ubs)
	Random.seed!(rseed)
	#create problem and run optimization

	#res = optimize(objectiveICs, lbs, ubs, ParticleSwarm(n_particles=40), Optim.Options(iterations=numFevals))
	@time res = optimize(objectiveICs,nominal_ICs, lbs, ubs,NelderMead(), Optim.Options(iterations=numFevals))
	print(res)
	writedlm(string("../parameterEstimation/ICEstimation/MinimumICsForPatient", patient_id, "usingParametersFrom_02_01_19UsingNM_W", numFevals,"Evals.txt"), res.minimizer)
	writedlm(string("../parameterEstimation/ICEstimation/ResultICsForPatient", patient_id, "usingParametersFrom_02_01_19UsingNM_W", numFevals,"Evals.txt"),string(res))
end
