using NLopt
using ODE
#using Optim
#include("DEBalanceEqns.jl")
include("utilities.jl")
include("Balances.jl")
include("CoagulationModelFactory.jl")

function runSim(params,ICs)
	tstart = 0.0
	tend = 120.0
	TS = .02
	TSIM = collect(tstart:TS:tend)
	PROBLEM_DICTIONARY = buildCompleteDictFromOneVector(params)
	initial_condition_vector =ICs
	#pull things out of dictionary
	kinetic_parameter_vector = PROBLEM_DICTIONARY["KINETIC_PARAMETER_VECTOR"]
	control_parameter_vector = PROBLEM_DICTIONARY["CONTROL_PARAMETER_VECTOR"]
	qualitative_factor_level_vector = PROBLEM_DICTIONARY["FACTOR_LEVEL_VECTOR"] #this is actually set, not adjustable, but need to feed it in
	platelet_parameter_vector = PROBLEM_DICTIONARY["PLATELET_PARAMS"]
	timing = PROBLEM_DICTIONARY["TIME_DELAY"]
	fibrin_kinetic_parameter_vector = PROBLEM_DICTIONARY["FIBRIN_KINETIC_PARAMETER_VECTOR"]
	fibrin_control_vector = PROBLEM_DICTIONARY["FIBRIN_CONTROL_PARAMETER_VECTOR"]
	FVIII_control = PROBLEM_DICTIONARY["FVIII_CONTROL"]
	aleph=PROBLEM_DICTIONARY["ALEPH"]

	allparams = vcat(kinetic_parameter_vector,control_parameter_vector,qualitative_factor_level_vector,platelet_parameter_vector,timing,fibrin_kinetic_parameter_vector,fibrin_control_vector, FVIII_control,aleph)
	fbalances(t,y)=Balances(t,y,PROBLEM_DICTIONARY)
	t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0)

	#prob=ODEProblem(BalanceEquations, ICs, (tstart,tend), allparams)
	#sol = solve(prob, CVODE_BDF(method=:Functional),reltol=1e-8,abstol=1e-8, dtmax=1.0,saveat = 0.01)
	#Plots.plot(sol)
	#tpa is 12
	tPA = initial_condition_vector[16]
	A = convertToROTEM(TSIM, X, tPA)
	#A = convertToROTEM(sol.t, sol, initial_condition_vector[16])
	
	#PyPlot.plot(sol.t, A)
	#@show CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(A,sol.t)
#	if(sol.retcode!=:Success)
#		sol = zeros(sol.t,22)
#	end
	return A, X,TSIM
end

function objective_f(initial_conditions::Vector, grad::Vector)
	@show initial_conditions
	#in this case, the params are the initial conditions we're adjusting
	#in this case, the params are the initial conditions we're adjusting
	A,X,TSIM = runSim(params,initial_conditions)
	R,CFT,alpha,MA,A10,A20,LI30,LI60=calculateCommonMetrics(A,TSIM)
	diffAlpha = alpha-target_alpha
	diffR = R-target_R
	diffMA=MA-target_MA
	sum = diffAlpha^2+diffR^2+ diffMA^2
	@show sqrt(sum), diffAlpha, diffR, diffMA
	return sqrt(sum)
	
end


allparams = readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt", '\t')
params = allparams[7,:]
#from https://www.sciencedirect.com/science/article/pii/S0049384805002434
#to be adjusted
target_alpha = 65.0
target_R = 400.0
target_MA = 64.0
PROBLEM_DICTIONARY = buildCompleteDictFromOneVector(params)
default_ICs = PROBLEM_DICTIONARY["INITIAL_CONDITION_VECTOR"]#set our variable size
default_ICs[[4,10,11,12,13,16,17,18,19,22]] = eps() #so they're technially larger than 0
lowerbounds = zeros(size(default_ICs)) #a patient could lack all species
#lowerbounds[1] = default_ICs[1]*.25 #but they must have some amount of probthrombin
upperbounds = 4*default_ICs
#make it possible for things that are zero to move
upperbounds[[4,10,11,12,13,16,17,18,19,22]] = 5.0
#setting bounds for opt problem

opt = Opt(:LN_NELDERMEAD, maximum(size(default_ICs))) 
#@show lower_bounds.<upper_bounds
lower_bounds!(opt, vec(lowerbounds))
upper_bounds!(opt, vec(upperbounds))
min_objective!(opt, objective_f)

@show opt
#run optimization
(minf,minx,ret) = optimize(opt, vec(default_ICs))
println("got $minf at $minx after $count iterations (returned $ret)")

#optimize(objective_f, default_ICs, lowerbounds, upperbounds)
