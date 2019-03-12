include("Balances.jl")
include("CoagulationModelFactory.jl")
include("utilities.jl")
include("Kinetics.jl")
include("Control.jl")
#include("plotData.jl")
#using Sundials
#using ODE
using DifferentialEquations
#using PyPlot
#using PyCall
#PyCall.PyDict(matplotlib["rcParams"])["font.sans-serif"] = ["Helvetica"]

function runModel(TSTART,Ts,TSTOP, platelet_count)
	#TSTART = 0.0
	#Ts = .02
	#TSTOP = 1.0
	PROBLEM_DICTIONARY = Dict()
	PROBLEM_DICTIONARY = buildCoagulationModelDictionary(platelet_count)
	TSIM = TSTART:Ts:TSTOP
	initial_condition_vector = PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"]
	reshaped_IC = vec(reshape(initial_condition_vector,22,1)) #may need to cast to vector for Sundials

	#calling solver
#	fbalances(t,y,ydot)=Balances(t,y,ydot,PROBLEM_DICTIONARY)
#	X = Sundials.cvode(fbalances,reshaped_IC,TSIM, abstol =1E-4, reltol=1E-4);
	TSIM = collect(TSTART:Ts:TSTOP)
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
#	println("got here")
	return (t,X);
end

function solveAdjBalances(TSTART,TSTOP,Ts,parameter_index, PROBLEM_DICTIONARY)
	TSIM = TSTART:Ts:TSTOP
	initial_condition_array = PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector=vcat(initial_condition_array, vec(zeros(22,1)))
	PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"] = initial_condition_vector
	@show size(initial_condition_vector)
	epsilon = 1e-6
	fbalances(t,y) = AdjBalances(t,y,parameter_index, PROBLEM_DICTIONARY)
	(t,y) =ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
	#Map
	  number_of_timesteps = length(t)
	  number_of_states = length(initial_condition_vector)
	  X = zeros(number_of_timesteps,number_of_states)
	  for state_index = 1:number_of_states
	    tmp = map(y->y[state_index],y)
	    for time_index = 1:number_of_timesteps
	      X[time_index,state_index] = tmp[time_index]
	    end
	  end

	  # # Check for smalls -
	  idx_n = find(abs(X).<epsilon)
	  X[idx_n] = 0.0

	  # return time and state -
  return (t,X);
	
end

#function runModelSundials(TSTART,Ts,TSTOP)

#	PROBLEM_DICTIONARY = Dict()
#	PROBLEM_DICTIONARY = buildCoagulationModelDictionary()
#	TSIM = collect(TSTART:Ts:TSTOP)
#	initial_condition_vector = PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"]
#	reshaped_IC = vec(reshape(initial_condition_vector,11,1)) #may need to cast to vector for Sundials

#	#calling solver
#	fbalances(t,y,ydot)=Balances(t,y,ydot,PROBLEM_DICTIONARY)
#	mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON) 
#	Sundials.@checkflag Sundials.CVodeSetMaxStep(mem, .01)
#	#X = Sundials.cvode(fbalances,reshaped_IC,TSIM, integrator=:Adams,reltol=1E-8, abstol=1E-8)#, abstol =1E-4, reltol=1E-4);
#	flag = Sundials.CVode(mem, TSIM,fbalances, TSIM, Sundials.CV_NORMAL)
##	println("got here")
#	return (X);
#end

function makePlots(t,x)
	FII = x[:,1]
	FIIA = x[:,2]
	PC = x[:,3]
	APC = x[:,4]
	ATIII = x[:,5]
	TM = x[:,6]
	TRIGGER = x[:,7]

	fig = figure(figsize = (10,10))
#	y_formatter = PyPlot.ticker.ScalarFormatter(useOffset=False)
#	ax = fig.gca()
#	println(ax)
#	ax.yaxis.set_major_formatter(y_formatter)
	plt[:subplot](2,4,1)
	plot(t, FII)
	title("FII")
	plt[:subplot](2,4,2)
	plot(t, FIIA)
	title("FIIa")
	plt[:subplot](2,4,3)
	plot(t, PC)
	title("PC")
	plt[:subplot](2,4,4)
	plot(t, APC)
	title("APC")
	plt[:subplot](2,4,5)
	plot(t, ATIII)
	title("ATIII")
	plt[:subplot](2,4,6)
	plot(t, TM)
	title("TM")
	plt[:subplot](2,4,7)
	plot(t, TRIGGER)
	title("TRIGGER")

	
end

function makeLoopPlots(t,x)
	names = ["FII", "FIIa", "PC", "APC", "ATIII", "TM", "TRIGGER", "Fraction Activated Platelets", "FV_FX", "FV_FXa", "Prothombinase-Platelets",
	"Fibrin", "Plasmin", "Fibrinogen", "Plasminogen", "tPA", "uPA", "fibrin monomer", "protofibril", "antiplasmin", "PAI_1", "Fiber"]
	#fig = figure(figsize = (15,15))
#	y_formatter = PyPlot.ticker.ScalarFormatter(useOffset=false)
#	ax = fig.gca()
#	println(ax)
#	ax.yaxis.set_major_formatter(y_formatter)
	#@show size(t)
	for j in collect(1:size(names,1))
		plt[:subplot](5,5,j)
		#@show size([a[j] for a in x])
		#plot(t, [a[j] for a in x], "k")
		plot(t, x[j,:], "k")
		title(names[j])
	end
	#savefig("figures/Dec19_BeforeOpt.pdf")
end

function plotFibrinSpecies(t,x)
	selectedidxs = [12,14,18,19,22]
	legarr = ["Fibrin", "fibrinogen", "fibrin monomer", "protofibril", "Fiber", "Sum of Clot Forming Species"]
	#sumofspecies = zeros(size(a[18] for a in x))
	sumofspecies = zeros(size(x[18,:]))
	for j in selectedidxs
			#@show j
			#semilogy(t, [a[j] for a in x])
			semilogy(t, x[j,:])
			if(j !=14)
				#sumofspecies= sumofspecies +[a[j] for a in x]
				sumofspecies =sumofspecies+x[j,:]
			end
	end
	semilogy(t, sumofspecies)
	legend(legarr, loc="best")
	
end

function makePlotsfromODE4s(t,x)
	FII = Float64[]
	FIIA = Float64[]
	PC = Float64[]
	APC = Float64[]
	ATIII = Float64[]
	TM = Float64[]
	TRIGGER= Float64[]
	frac_acviated = Float64[]
	for j = 1:length(x)
		currx = x[j]
		for k = 1:length(currx)
			curritem = currx[k]
			if(k == 1)
				push!(FII, curritem)
			elseif(k == 2)
				push!(FIIA,curritem)
			elseif(k == 3)
				push!(PC, curritem)
			elseif(k == 4)
				push!(APC, curritem)
			elseif(k == 5)
				push!(ATIII, curritem)
			elseif(k == 6)
				push!(TM, curritem)
			elseif(k ==7)
				push!(TRIGGER,curritem)
			elseif(k ==8)
				push!(frac_acviated, curritem)
			end
		end
	end
	fig = figure(figsize = (10,10))
	plt[:subplot](2,4,1)
	plot(t, FII)
	title("FII")
	plt[:subplot](2,4,2)
	plot(t, FIIA)
	title("FIIa")
	plt[:subplot](2,4,3)
	plot(t, PC)
	title("PC")
	plt[:subplot](2,4,4)
	plot(t, APC)
	title("APC")
	plt[:subplot](2,4,5)
	plot(t, ATIII)
	title("ATIII")
	plt[:subplot](2,4,6)
	plot(t, TM)
	title("TM")
	plt[:subplot](2,4,7)
	plot(t, TRIGGER)
	title("TRIGGER")
	plt[:subplot](2,4,8)
	plot(t, frac_acviated)
	title("Fraction of Platelets Activated")

end

function plotFluxes(pathToData,t)
	data = readdlm(pathToData)
	fig = figure(figsize = (15,15))
	for j in collect(1:size(data,2))
		plt[:subplot](size(data,2),1,j)
		@show size(t)
		@show size(data[:,j])
		currdata = data[:,j]
		currdata[currdata.>=1E6] = 0.0
		plot(t, currdata, ".k")
		title(string("reaction ", j))
	end
end


function main()
	pathToData = "../data/fromOrfeo_Thrombin_HT_PRP.txt"
	#HT-346, BL-420
	close("all")
	#rm("ratevector.txt")
	#rm("modifiedratevector.txt")
	#rm("times.txt")
	(t,x) = runModel(0.0, .01, 35.0, 346)
	@show size(t)
	#remove tiny elements that are causing plotting problems
	#x[x.<=1E-20] = 0.0
	#times = readdlm("times.txt")
	#plotFluxes("ratevector.txt",times)
	#plotFluxes("modifiedratevector.txt",times)

	#println(x)
	makeLoopPlots(t,x)
	plotThrombinWData(t,x,pathToData)
	#savefig("figures/BeforeNLOptFeb8.pdf")
	data = readdlm(pathToData)
	time = data[:,1]/60 #convert to minutes
	avg_run = mean(data[:,2:3],2);
	usefuldata = hcat(time, avg_run)
	MSE, interpolatedExperimentalData=calculateMSE(t, [a[2] for a in x], usefuldata)
	#figure()
	#plot(t, interpolatedExperimentalData)
	estimatedAUC = calculateAUC(t, [a[2] for a in x])
	experimentalAUC = calculateAUC(t, interpolatedExperimentalData)
	println("MSE: %f, AUC Difference %f", MSE, abs(estimatedAUC-experimentalAUC) )
	return MSE, abs(estimatedAUC-experimentalAUC)
end


function plotThrombinWData(t,x,pathToData)
	#close("all")
	#figure()
	data = readdlm(pathToData)
	time = data[:,1]
	avg_run = mean(data[:,2:3],dims=2);
	plotcolor = "k"
	#fig = figure(figsize = (15,15))
	#plot(t, [a[2] for a in x], "-", color = plotcolor)
	plot(t, x[2,:],"-", color = plotcolor)
	plot(time/60, avg_run, ".", color = plotcolor)
	ylabel("Thrombin Concentration, nM")
	xlabel("Time, in minutes")
	#savefig("figures/UsingNMParameters.pdf")
end

function plotAverageThrobinWData(t,meanThrombin,stdThrombin,pathToData,MSE, savestr)
	expdata = readdlm(pathToData,'\t')
	fig = figure(figsize = (15,15))
	ax = gca()
	ax[:tick_params]("both",labelsize=24) 
	plot(expdata[:,1]./60, (expdata[:,4]), ".k", markersize=20)
	ylabel("Thrombin Concentration, nM", fontsize=28)
	xlabel("Time, in minutes",fontsize =28)
	plot(t, transpose(meanThrombin), "k")
	axis([0, 35, 0, 250])
	@show size(meanThrombin)
	@show size(stdThrombin)
	@show size(t)
	upper = transpose(meanThrombin+1.96*stdThrombin)
	lower = transpose(meanThrombin-1.96*stdThrombin)
	@show size(vec(upper))
	@show size(vec(lower))
	fill_between((t), vec(upper), vec(lower), color = ".5")
#	annotate(string("MSE=", MSE),
#	xy=[.85;.85],
#	xycoords="figure fraction",
#	xytext=[0,0],
#	textcoords="offset points",
#	ha="right",
#	va="top")
	savefig(savestr)
end

function plotAverageThrobinWData(t,meanThrombin,stdThrombin,expdata, savestr)
	fig = figure(figsize = (15,15))
	println("here")
	ylabel("Thrombin Concentration, nM", fontsize = 20)
	xlabel("Time, in minutes", fontsize = 20)
	plot(t, transpose(meanThrombin), "k")
	axis([0, 35, 0, 200])
	@show size(meanThrombin)
	@show size(stdThrombin)
	@show size(t)
	upper = transpose(meanThrombin+stdThrombin)
	lower = transpose(meanThrombin-stdThrombin)
	@show size(vec(upper))
	@show size(vec(lower))
	fill_between((t), vec(upper), vec(lower), color = ".5", alpha =.5)
	plot(expdata[:,1], expdata[:,2], ".k")
	savefig(savestr)
end

function runModelWithMultipleParams(pathToParams,pathToData,savestr)
	close("all")
	allparams = readdlm(pathToParams, '\t')
	TSTART = 0.0
	Ts = .02
	TSTOP =180.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Luan2010Fig5F.csv"
	#pathToData = "../data/fromOrfeo_Thrombin_BL_PRP.txt"
	data = readdlm(pathToData)
	time = data[:,1]
	avg_run = mean(data[:,2:end],2);
	usefuldata = hcat(time/60, avg_run)
	fig1 = figure(figsize = (15,15))
	fig2 = figure(figsize = (15,15))
	fig3 = figure(figsize = (15,15))
	platelet_count =418
	tPA = 0
	alldata = zeros(1,size(TSIM,1))
	@show size(allparams)
	if(size(allparams,1)==46) #deal with parameters being stored either vertically or horizontally
		itridx = 2
	else
		itridx = 1
	end
	
	for j in collect(1:size(allparams,itridx))
		if(itridx ==2)
			currparams = allparams[:,j]
		else
			currparams = allparams[j,:]
		end
		@show currparams
		push!(currparams, platelet_count)
		dict = buildDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(t,y)= Balances(t,y,dict) 
		t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-4, reltol = 1E-4, minstep = 1E-9, points=:specified)
		figure(1)
		plotThrombinWData(t,X,pathToData)
		figure(2)
		makeLoopPlots(t,X)
		#@show alldata
		#@show size([a[2] for a in X])
		A = convertToROTEM(t,X,tPA)
		figure(3)
		plot(t, A)
		alldata=vcat(alldata,transpose([a[2] for a in X]))
	end
	figure(3)
	#plotROTEM_given_tPA(tPA)
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	#hasdynamics = checkForDynamics(alldata)
	#idx = find(hasdynamics->(hasdynamics==1),hasdynamics); #find the indices for the parameter sets that actually produce dynamics
	#alldata = alldata[idx,:]
	meanThrombin = mean(alldata,1)
	stdThrombin = std(alldata,1)
	plotAverageThrobinWData(TSIM, meanThrombin, stdThrombin, usefuldata,savestr)
	return alldata
end

function runModelWithMultipleParams(pathToParams,pathToData,index,savestr)
	close("all")
	allparams = readdlm(pathToParams, ',')
	TSTART = 0.0
	Ts = .02
	TSTOP = 90.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Luan2010Fig5F.csv"
	fig = figure(figsize = (15,15))
	alldata = zeros(1,size(TSIM,1))
	
	for j in collect(1:size(allparams,1))
		currparams = allparams[j,:]
		dict = buildDictFromOneVector(currparams)
		dict = createCorrectDict(dict, index)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		fbalances(t,y)= Balances(t,y,dict) 
		t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		@show size(TSIM)
		@show size(t)
		@show size(X)
		@show size(alldata)
		@show size(transpose([a[2] for a in X]))
		alldata=vcat(alldata,transpose([a[2] for a in X]))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanThrombin = mean(alldata,1)
	stdThrombin = std(alldata,1)
	MSE, interpolatedExperimentalData=calculateMSE(TSIM, transpose(meanThrombin), readdlm(pathToData, ','))
	@show MSE
	plotAverageThrobinWData(TSIM, meanThrombin, stdThrombin, pathToData, MSE,savestr)
	return alldata
end

function runModelWithParams(params)
	tic()
	close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP = 60
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
	params[47]=curr_platelets
	dict = buildCompleteDictFromOneVector(params)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector[16]=tPA
	#fbalances(t,y)= Balances(t,y,dict) 
	fbalances(t,y)= Balances(t,y,dict)
	toc()
	tic() 
	t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
	toc()
#	plotThrombinWData(t,X,pathToData)
#	figure()
#	plotFibrinSpecies(t,X)
	A = convertToROTEM(t,X,tPA)
#	figure()
#	plot(t, A)
#	plot(usefulROTEMdata[:,1], usefulROTEMdata[:,2], "k.")
#	#savefig("figures/AfterNM_24_03_2017.pdf")
#	figure(figsize=[15,15])
#	makeLoopPlots(t,X)
	MSE, interpolatedExperimentalData=calculateMSE(t, [a[2] for a in X], usefuldata)
	ROTEM_MSE, interpROTEM = calculateMSE(t, A, usefulROTEMdata)
	return MSE, ROTEM_MSE
end


function runModelWithParamsReturnAUC(params,tPA)
	close("all")
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP=90.0
	else
		TSTOP=180
	end
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	dict = buildCompleteDictFromOneVector(params)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector[16]=tPA
	fbalances(y,p,t)= Balances(t,y,dict) 
	#fbalances(t,y)= Balances(t,y,dict) 
	#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	A = convertToROTEM(t,X,tPA)
	AUC=calculateAUC(t, A)
	return AUC
end

function runModelWithParamsReturnA(params,tPA)
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP=90.0
	else
		TSTOP=180.0
	end
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	dict = buildCompleteDictFromOneVector(params)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector[16]=tPA
	platelet_count = params[47]
	fbalances(y,p,t)= Balances(t,y,dict) 
	#fbalances(t,y)= Balances(t,y,dict) 
	#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	A = convertToROTEMPlateletContribution(t,X, tPA,platelet_count)
	#AUC=calculateAUC(t, A)
	return t,A
end

function runModelWithParamsSetICReturnAUC(params)
	close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP=120.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	modelparams = params[1:77]
	dict = buildCompleteDictFromOneVector(modelparams)
	initial_condition_vector = params[78:end]
	tPA = initial_condition_vector[16]
	fbalances(t,y)= Balances(t,y,dict) 
	t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)
	A = convertToROTEM(t,X,tPA)
	AUC=calculateAUC(t, A)
	return AUC
end

function runModelWithParamsChangeICReturnAUC(params, currIC)
	TSTART = 0.0
	Ts = .02
	TSTOP=60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	modelparams = params[1:77]
	dict = buildCompleteDictFromOneVector(modelparams)
	initial_condition_vector = currIC
	tPA = initial_condition_vector[16]
	TSIM = collect(TSTART:Ts:TSTOP)
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	A = convertToROTEM(t,X,tPA)
	AUC=calculateAUC(t, A)
end

function runModelWithParamsChangeICReturnA(params, currIC)
	TSTART = 0.0
	Ts = .02
	TSTOP=60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	modelparams = params[1:77]
	curr_platelets=params[47]
	dict = buildCompleteDictFromOneVector(modelparams)
	initial_condition_vector = currIC
	tPA = initial_condition_vector[16]
	TSIM = collect(TSTART:Ts:TSTOP)
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	#A = convertToROTEM(t,X,tPA)
	A =  convertToROTEMPlateletContribution(t,X, tPA,curr_platelets)
	return t,A
end

function runModelWithParamsChangeICReturnA(params,genIC,genExp,genPlatelets)
	TSTART = 0.0
	Ts = .02
	TSTOP=60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	modelparams = params[1:77]
	curr_platelets=genPlatelets
	initial_condition_vector=genIC
	dict = buildCompleteDictFromOneVector(modelparams)
	tPA = genIC[16]
	dict["FACTOR_LEVEL_VECTOR"]=genExp
	TSIM = collect(TSTART:Ts:TSTOP)
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	#A = convertToROTEM(t,X,tPA)
	A =  convertToROTEMPlateletContribution(t,X, tPA,curr_platelets)
	return t,A
end

function runModelWithParamsChangeICReturnA(params)
	TSTART = 0.0
	Ts = .02
	TSTOP=60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	curr_platelets=params[47]
	modelparams = params
	dict = buildCompleteDictFromOneVector(modelparams)
	initial_condition_vector = currIC
	tPA = initial_condition_vector[16]
	TSIM = collect(TSTART:Ts:TSTOP)
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	A =  convertToROTEMPlateletContribution(t,x, tPA,curr_platelets)
	return t,A
end

function runModelWithParamsSetICReturnROTEM(params)
	close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP=60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	modelparams = params[1:77]
	dict = buildCompleteDictFromOneVector(modelparams)
	initial_condition_vector = initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	#tPA = initial_condition_vector[16]
	tPA = 2.0
	initial_condition_vector[16]=tPA
	TSIM = collect(TSTART:Ts:TSTOP)
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	A = convertToROTEM(t,X,tPA)
	return t,A
end

function runModelWithParamsSetICReturnAUCFibrinIC(params, fibrinIC)
	close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP=120.0
	TSIM = collect(TSTART:Ts:TSTOP)
	#curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"5")
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	#use default platelets	
	#params[47]=curr_platelets
	modelparams = params[1:77]
	dict = buildCompleteDictFromOneVector(modelparams)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	tPA = initial_condition_vector[16]
	TSIM = collect(TSTART:Ts:TSTOP)
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob)
	t =sol.t
	X = sol
	A = convertToROTEM(t,X,tPA)
	AUC=calculateAUC(t, A)
	return AUC
end






function runModelWithParamsPeturbIC(params, num_runs)
	close("all")
	#rm("dataforvarner.txt")
	savestr = string("../figures/With_", num_runs, "different_IC.png")
	TSTART = 0.0
	Ts = .02
	TSTOP = 90.0
	TSIM = collect(TSTART:Ts:TSTOP)
	pathToData = "../data/fromOrfeo_Thrombin_BL_PRP.txt"
	#pathToData = "../data/Luan2010Fig5A.csv"
	fig = figure(figsize = (15,15))
	alldata = zeros(1,size(TSIM,1))
	seeds = [24,101,1000,3,4,5,11,14,17,23423,13124,123235,1232,132234,33,45,345,456,12434,100,101,102,105,109,111,1111,22,2222,33,3333,44,4444,55,5555,66,6666,66,77,777,777]
	
	for j in collect(1:num_runs)
		currparams = params
		dict = buildDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector = peturbIC(initial_condition_vector, j)
		fbalances(t,y)= Balances(t,y,dict) 
		t,X = ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6,points=:specified)
		@show size(TSIM)
		@show size(t)
		@show size(X)
		@show size(alldata)
		@show size(transpose([a[2] for a in X]))
		alldata=vcat(alldata,transpose([a[2] for a in X]))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanThrombin = mean(alldata,1)
	stdThrombin = std(alldata,1)
	MSE, interpolatedExperimentalData=calculateMSE(TSIM, transpose(meanThrombin), readdlm(pathToData, '\t'))
	#@show MSE
	plotAverageThrobinWData(TSIM, meanThrombin, stdThrombin, pathToData, MSE,savestr)
	#f = open("dataforvarner.txt", "a+")
	#writedlm(f, transpose(TSIM), ',')
	#writedlm(f, (meanThrombin), ',')
	#writedlm(f, (stdThrombin), ',')
	#close(f)
	return alldata
end

function runModelWithParamsSetF8(params, FVIIIcontrol, index)
	close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP = 90.0
	TSIM = collect(TSTART:Ts:TSTOP)
	letters = ["A", "B", "C", "D", "E", "F"]
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	pathToData = string("../data/Luan2010Fig5",letters[index], ".csv")
	fig = figure(figsize = (15,15))
	
	dict = buildDictFromOneVector(params)
	dict = createCorrectDict(dict, index)
	dict["FVIII_CONTROL"]= FVIIIcontrol
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector = setIC(initial_condition_vector, index)
	@show initial_condition_vector
	@show dict["FVIII_CONTROL"]
	fbalances(t,y)= Balances(t,y,dict) 
	t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6)
	plotThrombinWData(t,X,pathToData)
	savefig(string("figures/AttemtingF8FittingSet",index,"_02_16_2017.pdf"))
	#makeLoopPlots(t,X)
	MSE, interpolatedExperimentalData=calculateMSE(t, [a[2] for a in X], readdlm(pathToData, ','))
	return MSE
end

function runModelWithParamsSetF8OnePlot(fig,params, FVIIIcontrol, index)
	#close("all")
	TSTART = 0.0
	Ts = .02
	TSTOP = 90.0
	TSIM = collect(TSTART:Ts:TSTOP)
	letters = ["A", "B", "C", "D", "E", "F"]
	#pathToData = "../data/ButenasFig1B60nMFVIIa.csv"
	#pathToData = "../data/Buentas1999Fig4100PercentProthrombin.txt"
	pathToData = string("../data/Luan2010Fig5",letters[index], ".csv")
	
	dict = buildDictFromOneVector(params)
	dict = createCorrectDict(dict, index)
	dict["FVIII_CONTROL"]= FVIIIcontrol
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector = setIC(initial_condition_vector, index)
	@show initial_condition_vector
	@show dict["FVIII_CONTROL"]
	fbalances(t,y)= Balances(t,y,dict) 
	t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6)
	plotThrombinWData(t,X,pathToData, string(index/(6+.1)))
	#savefig(string("figures/AttemtingF8FittingSet",index,"_02_16_2017.pdf"))
	#makeLoopPlots(t,X)
	MSE, interpolatedExperimentalData=calculateMSE(t, [a[2] for a in X], readdlm(pathToData, ','))
	return fig
end

function makeThrominFigureForPoster()
	allparams = readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt", '\t')
	runModelWithParamsPeturbIC(allparams[4,:],30)

end
