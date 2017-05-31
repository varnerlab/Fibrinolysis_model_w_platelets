include("BalanceEquations.jl")
include("CoagulationModelFactory.jl")
include("utilities.jl")
#using Sundials
using ODE

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

function objective_function(parameter_array, selected_sets)
	#write to file current parameters
	for j in collect(1:size(parameter_array,1)) #enforce lower bound
		if(parameter_array[j]<SMALL)
			parameter_array[j] = SMALL
		end
	end

#	f = open("parameterEstimation/parameters.txt", "a+")
#	writedlm(f, string(parameter_array))
#	close(f)

	obj_array = BIG*ones(5,1)
	@show parameter_array
	TSTART = 0.0
	Ts = .02
	TSTOP = 60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	counter = 1
	for j in selected_sets
		dict = buildDictFromOneVector(parameter_array)
		dict = createCorrectDict(dict, j)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		fbalances(t,y)= BalanceEquations(t,y,dict) 
		tic()
		t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-4, reltol = 1E-6)
		toc()
		FIIa = [a[2] for a in X]
		MSE, interpolatedExperimentalData=calculateMSE(t, FIIa, allexperimentaldata[j])
		obj_array[counter,1] = MSE
		counter = counter+1
	end
	@show obj_array
	return obj_array
end


