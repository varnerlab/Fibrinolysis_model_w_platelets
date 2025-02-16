#for comparing generated initial conditions to found initial conditions

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

function runModelCompareConditions(givenICparams, foundICparams)
	# Load data -
	global kin_params = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	d = buildCompleteDictFromOneVector(kin_params)
	close("all")
	#run for givenIC
	curr_exp = givenICparams[1:8]
	curr_ICs = givenICparams[9:end-1]
	curr_platelets = givenICparams[end]
	tPA = givenICparams[12]
	println("Running Simulation with given parameters")
	Tgiven,Rgiven =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
	#run for foundIC
	curr_exp = foundICparams[1:8]
	curr_ICs = foundICparams[9:end-1]
	curr_platelets = foundICparams[end]
	tPA = foundICparams[12]
	println("Running Simulation with estimated parameters")
	Tfound,Rfound =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)

	figure(figsize = [15,15])
	@show Rgiven
	plot(Tgiven, Rgiven, "k-")
	plot(Tfound, Rfound, "k--")
	legend(["Given Initial Conditions", "Estimated Initial Conditions"])
	savefig("../figures/ComparingICs.pdf")
end

function runModelAllICs(numICs)
	close("all")
	figure(figsize = [15,15])
	global kin_params = readdlm("../parameterEstimation/best8_02_05_18.txt")[1,:]
	d = buildCompleteDictFromOneVector(kin_params)
	for j = 1:numICs
		givenICparams = readdlm(string("../solveInverseProb/ics_to_match_29_07_18_hyperfibrinolysis_iter", j,".txt"), ',')
		curr_exp = givenICparams[1:8]
		curr_ICs = givenICparams[9:end-1]
		curr_platelets = givenICparams[end]
		tPA = givenICparams[12]
		println("Running Simulation with given parameters")
		Tgiven,Rgiven =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
		plot(Tgiven, Rgiven, "k-")
	end
end

#foundICs = readdlm("../solveInverseProb/foundIcs_24_07_18_Hyercoag1.txt")
#givenICs = readdlm("../solveInverseProb/ics_to_match_24_07_18_Hyercoag_iter1.txt", ',')

#runModelCompareConditions(givenICs, foundICs)
