include("Balances.jl")
include("CoagulationModelFactory.jl")
include("utilities.jl")
include("plotData.jl")
include("runModel.jl")
#using Sundials
using ODE
using PyPlot
using PyCall
PyCall.PyDict(matplotlib["rcParams"])["font.sans-serif"] = ["Helvetica"]

function basicRunModel()
	#allparams = readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt", '\t')
	#params = allparams[4,:]
	params = vec(readdlm("../parameterEstimation/startingPoint_02_05_18.txt"))
	close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP = 90.0
	TSIM = collect(TSTART:Ts:TSTOP)
	tPA = 2.0
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	pathToData = "../data/fromOrfeo_Thrombin_BL_PRP.txt"
	
	data = readdlm(pathToData)
	time = data[:,1]
	avg_run = mean(data[:,2:3],2);
	usefuldata = hcat(time, avg_run)

	curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"6")
	fig = figure(figsize = (15,15))
	@show curr_platelets
	params[47]=curr_platelets
	dict = buildCompleteDictFromOneVector(params)
	@show dict["PLATELET_PARAMS"]
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector[16]=tPA
	#fbalances(t,y)= Balances(t,y,dict) 
	fbalances(t,y)= Balances(t,y,dict) 
	t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
	plotThrombinWData(t,X,pathToData)
	figure()
	plotFibrinSpecies(t,X)
	A = convertToROTEMPlateletContribution(t,X,tPA,curr_platelets)
	#A = convertToROTEMPlateletContributionScaledF(t,X,tPA,curr_platelets, initial_condition_vector[14])
	figure()
	plot(t, A)
	plot(usefulROTEMdata[:,1], usefulROTEMdata[:,2], "k.")
	#savefig("figures/AfterNM_24_03_2017.pdf")
	figure(figsize=[15,15])
	makeLoopPlots(t,X)
end

function basicRunModel(params)
	close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP = 60
	TSIM = collect(TSTART:Ts:TSTOP)
	tPA = 2.0
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#pathToData = "../data/fromOrfeo_Thrombin_BL_PRP.txt"
	pathToData = "../data/fromOrfeo_Thrombin_HT_PRP.txt"
	#for observing thrombin generation
	#HT = 364, BL 420
	HT_platelets = 364
	BL_platelets = 420	


	data = readdlm(pathToData)
	time = data[:,1]
	avg_run = mean(data[:,2:3],2);
	usefuldata = hcat(time, avg_run)

	curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"6")
	curr_platelets = HT_platelets
	fig = figure(figsize = (15,15))
	params[47]=curr_platelets
	dict = buildCompleteDictFromOneVector(params)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	#temp adjustments for HT_platelets
	initial_condition_vector[1] = initial_condition_vector[1]*.7
	initial_condition_vector[3] = initial_condition_vector[3]*1.3
	initial_condition_vector[5] = initial_condition_vector[5]*1.3
	initial_condition_vector[6] = initial_condition_vector[6]*1.3
	initial_condition_vector[9] = initial_condition_vector[9]*.7

	initial_condition_vector[16]=tPA
	#fbalances(t,y)= Balances(t,y,dict) 
	fbalances(t,y)= Balances(t,y,dict) 
	t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
	plotThrombinWData(t,X,pathToData)
	figure()
	plotFibrinSpecies(t,X)
	A = convertToROTEMPlateletContribution(t,X,tPA,curr_platelets)
	#A = convertToROTEMPlateletContributionScaledF(t,X,tPA,curr_platelets, initial_condition_vector[14])
	MSE, interpData = calculateMSE(t,A, usefulROTEMdata)
	@show MSE
	figure()
	plot(t, A)
	plot(usefulROTEMdata[:,1], usefulROTEMdata[:,2], "k.")
	#savefig("figures/AfterNM_24_03_2017.pdf")
	figure(figsize=[15,15])
	makeLoopPlots(t,X)
	AnyFlats=checkForFlatness(t,A)
	@show AnyFlats
	return t,A
end


