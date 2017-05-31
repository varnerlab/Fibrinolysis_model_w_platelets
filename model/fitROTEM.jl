include("utilities.jl")
include("BalanceEquations.jl")
include("CoagulationModelFactory.jl")
include("utilities.jl")
#using Sundials
using ODE
using NLopt

#load data once
tPA = 2.0
curr_platelets,usefulROTEMdata = setROTEMIC(tPA,"6")

function objectiveForROTEM(params::Vector, grad::Vector)
	tic()
	TSTART = 0.0
	Ts = .02
	TSTOP = 60
	TSIM = collect(TSTART:Ts:TSTOP)
	dict = buildCompleteDictFromOneVector(params)
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector[16]=tPA
	fbalances(t,y)= BalanceEquations(t,y,dict) 
	t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = .10)
	A = convertToROTEM(t,X,tPA)
	ROTEM_MSE, interpROTEM = calculateMSE(t, A, usefulROTEMdata)
	f = open("parameterEstimation/NM_ROTEM_24_05_2017.txt", "a+")
	write(f, string(params, ",", ROTEM_MSE, "\n"))
	close(f)
	toc()
	@show ROTEM_MSE
	return ROTEM_MSE
end

function attemptOptimizationNLOptROTEM(initial_parameter_estimate)
	numvars = 77
	outputfile = "parameterEstimation/fittingROTEM_05_24_17.txt"
	global up_arr = vcat(initial_parameter_estimate[1:46]*1.05, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/1.05, initial_parameter_estimate[47:end]/1000)
	#set platelets
	up_arr[47]=curr_platelets
	lb_arr[47]=curr_platelets
	#define opt problem
	opt = Opt(:LN_NELDERMEAD,numvars)

	#set bounds for opt problem
	upper_bounds!(opt, vec(up_arr))
	lower_bounds!(opt, vec(lb_arr))
	#set objective
	min_objective!(opt, objectiveForROTEM)
	(minf, minx, ret) = NLopt.optimize(opt, vec(initial_parameter_estimate))
	println("got $minf at $minx after $count iterations (returned $ret)")
	f = open(outputfile, "a")
	write(f, string(minf, ",", minx, ",", ret, "\n"))
	close(f)
	return minx
end
