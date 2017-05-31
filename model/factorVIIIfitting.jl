include("BalanceEquations.jl")
include("CoagulationModelFactory.jl")
include("utilities.jl")
include("runModel.jl")
#using Sundials
using ODE
using NLopt

# some global parameters -
BIG = 1e10
SMALL = 1e-9

#globally load data so that we don't have to reload on each iteration
pathsToData = ["../data/Luan2010Fig5A.csv", "../data/Luan2010Fig5B.csv","../data/Luan2010Fig5C.csv", "../data/Luan2010Fig5D.csv", "../data/Luan2010Fig5E.csv","../data/Luan2010Fig5F.csv"]
allexperimentaldata = Array[]
for j in collect(1:size(pathsToData,1))
	currdata = readdlm(pathsToData[j], ',')
	push!(allexperimentaldata, currdata)
end


function objectiveForNLOpt(params::Vector, grad::Vector)
	#tic()
	TSTART = 0.0
	Ts = .02
	TSTOP = 60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	dict["FVIII_CONTROL"] =params[1]
	@show params
	fbalances(t,y)= BalanceEquations(t,y,dict)
	#println("got here") 
	t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6)
	FIIa = [a[2] for a in X]
	MSE, interpolatedExperimentalData=calculateMSE(t, FIIa, experimentaldata)
	#hold("on")
	#plot(t, FIIa, alpha = .5)
	#write params to file
	f = open(savestr, "a+")
	write(f, string(params, ",", MSE, "\n"))
	close(f)
	#toc()
	@show MSE
	return MSE
end

function attemptOptimizationNLOptF8(exp_index,dict)
	numvars = 1
	opt = Opt(:LN_NELDERMEAD, numvars)
	lower_bounds!(opt, vec(fill(1E-9,1,numvars)))
	upper_bounds!(opt, vec(fill(1E2,1,numvars)))
	min_objective!(opt, objectiveForNLOpt)
	dict = createCorrectDict(dict, exp_index)
	(minf, minx, ret) = NLopt.optimize(opt, vec([50]))
	println("got $minf at $minx after $count iterations (returned $ret)")
	return minf, minx, ret
	
	
end

function runOpt(exp_index)
	global experimentaldata = allexperimentaldata[exp_index]
	global savestr = string("parameterEstimation/fitF8/NM15022017EstFVIIIparamDataSet",exp_index, ".txt")
	params = readdlm("parameterEstimation/AfterCalculatingF8FunctionPostNM.txt", ',')
	global dict =buildDictFromOneVector(params)
	dict = createCorrectDict(dict, exp_index)
	minf, minx, ret=attemptOptimizationNLOptF8(exp_index, dict)
	return minf, minx, ret
end

function runAllOpt()
	results = Array{Float64}(5,2)
	for j in collect(2:6)
		minf, minx, ret=runOpt(j)
		results[j-1,1] = minf
		results[j-1,2] = minx[1]
	end
	f8 = [1.07,.39,.07,.01, 0.0, 1.0]
	f8control = results[:,2]
	push!(f8control, 5.2089)
	figure()
	semilogy(f8, f8control, "k.")
	xlabel("FVIII Level")
	ylabel("FVIII Control")
	savefig("figures/FVIIIControl.pdf")	
	return results
end

function makeFigures()
	params = readdlm("parameterEstimation/AfterCalculatingF8FunctionPostNM.txt", ',')
	indices = [2,3,4,5,6,1]
	f8 = [1.07,.39,.07,.01, 0.0, 1.0]
	j = 1
	for ind in indices
		MSE=runModelWithParamsSetF8(params, calculateF8control(f8[j]), ind)
		@show MSE
		j = j+1
	end
end

function makeOneFigure()
	close("all")
	params = readdlm("parameterEstimation/AfterCalculatingF8FunctionPostNM.txt", ',')
	indices = [2,3,4,5,6,1]
	f8 = [1.07,.39,.07,.01, 0.0, 1.0]
	j = 1
	fig = figure(figsize = (15,15))
	hold(true)
	for ind in indices
		fig=runModelWithParamsSetF8OnePlot(fig, params, calculateF8control(f8[j]), ind)
		#@show MSE
		j = j+1
	end
	#legend(["FVIII=107%", "FVIII=100%", "FVIII=39%", "FVIII=7%", "FVIII=1%", "FVIII=0%"])
	axis([0, 60, 0, 400])
	savefig("figures/MasterFigure_08_03_2017.pdf")
end

function calculateF8control(f8)
	control = 4.6863*f8^2+.5357*f8-.01319
	return abs(control)
end
