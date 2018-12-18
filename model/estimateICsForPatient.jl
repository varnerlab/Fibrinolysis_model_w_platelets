include("runModel.jl")
using Optim
#to estimate initial conditions for a patient
#let's load the things we're going to want globally

POETs_data = "../parameterEstimation/POETS_info_05_12_18_PlateletContributionToROTEMFlatness1SmallerConversion.txt"

ec,pc,ra=parsePOETsoutput(POETs_data)
numParamSets = 2
bestparams=generateNbestPerObjective(numParamSets,ec,pc)

function objectiveICs(params)
	tPAs = [0,2]
	totalMSE = 0.0
	for tPA in tPAs
		#run model with current ICs
		alldata, meanROTEM, stdROTEM, TSIM=testROTEMPredicitionGivenParamsPlatetContributionToROTEM(bestparams, patient_id ,tPA, "findingICs.txt", params)
		@show size(meanROTEM)
		#calculate error
		platelets,expdata = setROTEMIC(tPA, patient_id)
		currMSE,interpolatedExperimentalData = calculateMSE(TSIM, meanROTEM, expdata)
		@show currMSE
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
		fbalances(t,y)= Balances(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		#@show size([a[2] for a in X])
		A = convertToROTEMPlateletContribution(t,X,tPA,platelet_count)
		@show size(A), size(TSIM)
		while(size(A)!=size(TSIM))
			tic() 
			t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
			toc()	
			#@show size([a[2] for a in X])
			A = convertToROTEMPlateletContribution(t,X,tPA,platelet_count)
		end

		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function runOptimization(currID)
	seed =89239
	numFevals =15
	global patient_id = currID #change me as nesseccary
	kin_params = mean(readdlm("../parameterEstimation/best8_02_05_18.txt"),1)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	lbs = nominal_ICs*.5
	ubs = nominal_ICs*1.5
	@show typeof(lbs)
	@show typeof(ubs)
	srand(seed)
	#create problem and run optimization
	tic()
	#res = optimize(objectiveICs, lbs, ubs, ParticleSwarm(n_particles=40), Optim.Options(iterations=numFevals))
	res = optimize(objectiveICs,nominal_ICs, lbs, ubs,NelderMead(), Optim.Options(iterations=numFevals))
	toc()
	print(res)
	writedlm(string("../parameterEstimation/ICEstimation/MinimumICsForPatient", patient_id, "usingParametersFrom_05_12_18UsingNM_W", numFevals,"Evals.txt"), res.minimizer)
	writedlm(string("../parameterEstimation/ICEstimation/ResultICsForPatient", patient_id, "usingParametersFrom_05_12_18UsingNM_W", numFevals,"Evals.txt"),string(res))
end
