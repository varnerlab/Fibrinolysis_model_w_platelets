include("BalanceEquations.jl")
include("CoagulationModelFactory.jl")
include("utilities.jl")
#using Sundials
using ODE

# some global parameters -
BIG = 1e10
SMALL = 1e-6

pathToData = "../data/ButenasFig1B60nMFVIIa.csv"

# Globally load data so we don't have load on each iteration -
MEASURED_ARRAY = readdlm(pathToData, ',')

function objective_function(parameter_array)
	#write to file current parameters
	for j in collect(1:size(parameter_array,1)) #enforce lower bound
		if(parameter_array[j]<SMALL)
			parameter_array[j] = SMALL
		end
	end

#	f = open("parameterEstimation/parameters.txt", "a+")
#	writedlm(f, string(parameter_array))
#	close(f)

	obj_array = BIG*ones(2,1)
	dict = buildDictFromOneVector(parameter_array)
	TSTART = 0.0
	Ts = .02
	TSTOP = 20.0
	TSIM = collect(TSTART:Ts:TSTOP)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	reshaped_IC = vec(reshape(initial_condition_vector,11,1))
	fbalances(t,y)= BalanceEquations(t,y,dict) 
	t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-8, reltol = 1E-8)
	FIIa = [a[2] for a in X]
	MSE, interpolatedExperimentalData=calculateMSE(t, FIIa, MEASURED_ARRAY)
	estimatedAUC = calculateAUC(t, FIIa)
	experimentalAUC = calculateAUC(t, interpolatedExperimentalData)
	obj_array[1,1] = MSE
	obj_array[2,1] = abs(estimatedAUC-experimentalAUC)
	return obj_array
end


