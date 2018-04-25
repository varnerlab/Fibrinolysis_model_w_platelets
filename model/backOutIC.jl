include("utilities.jl")
include("CoagulationModelFactory.jl")
include("BalanceEquations.jl")

using ODE
using PyPlot

function estimateICFromROTEM(target)
	patientID = 3
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
	numparams = 77
	allparams = readdlm(pathToParams, '\t')

	#let's use set 4
	j =4

	if(size(allparams,1)==numparams) #deal with parameters being stored either vertically or horizontally
		itridx = 2
	else
		itridx = 1
	end

	if(itridx ==2)
			currparams = vec(allparams[:,j])
		else
			currparams = vec(allparams[j,:])
		end
	dict = buildCompleteDictFromOneVector(currparams)

	alldata, meanROTEM, stdROTEM, TSIM=generateAvgROTEMCurve(patientID,allparams, [],0.0,numparams)
	curr_curve=calculateCommonMetrics(meanROTEM,TSIM)
	#plot(TSIM,meanROTEM)
	
	#proably just want to do grid search-but expensive-at .5,1,1.5 currparams-(22)^3
	#also, normal ranges-what are they
	return alldata, meanROTEM, stdROTEM, TSIM
	
	
end

function runSobolChangeIC()
	patientID = 3
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
	numparams = 77
	allparams = readdlm(pathToParams, '\t')

	#let's use set 4
	j =4

	if(size(allparams,1)==numparams) #deal with parameters being stored either vertically or horizontally
		itridx = 2
	else
		itridx = 1
	end

	if(itridx ==2)
			currparams = vec(allparams[:,j])
		else
			currparams = vec(allparams[j,:])
		end
	dict = buildCompleteDictFromOneVector(currparams)
	A,TSIM = generateROTEMCurve(patient_id,dict)
	curr_curve=calculateCommonMetrics(A,TSIM)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function generateROTEMCurve(patient_id,dict)
	pathToThrombinData="../data/fromOrfeo_Thrombin_HT_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 60.0
	end
	TSIM = collect(TSTART:Ts:TSTOP)
	platelets,usefuldata = setROTEMIC(tPA, patient_id)
	fbalances(t,y)= BalanceEquations(t,y,dict)
	tic() 
	t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
	toc()	
	A = convertToROTEM(t,X,tPA)
	return A,TSIM
	
end

function generateAvgROTEMCurve(patient_id,d)
	pathToThrombinData="../data/fromOrfeo_Thrombin_HT_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 60.0
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
		@show currparams
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		initial_condition_vector =alterIC(initial_condition_vector,IC_to_alter)
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(t,y)= BalanceEquations(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		A = convertToROTEM(t,X,tPA)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	return alldata, meanROTEM, stdROTEM, TSIM
end

