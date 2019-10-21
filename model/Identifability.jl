#for parameter identifability analysis
include("runModel.jl")

#using methods from Parameter identifiability analysis and visualization in large-scale kinetic models of biosystems

function calculateCentralDifferences()
	#load our estimated parameters
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	kin_params=mean(allp, dims=1)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]

	TSTART = 0.0
	Ts = .1 #set by saveat in runModelWithParamsChangeICReturnA
	TSTOP=60.0
	TSIM = collect(TSTART:Ts:TSTOP)
	total_num_params = maximum(size(nominal_ICs))+maximum(size(nominal_experimental))+1 #+1 for platelet count
	S = zeros(total_num_params, length(TSIM))

	#for now, let's do all of them
	param_idx = 1
	for j in 1:maximum(size(nominal_ICs))
		@show param_idx
		if(nominal_ICs[j]==0)
			#wiggle things set to zero
			p_plus = 2.0
			p_minus = 0.0
		else
			p_plus = nominal_ICs[j]*1.5
			p_minus = nominal_ICs[j]*.5
		end
		@show p_plus, p_minus

		h = p_plus - p_minus
		ICs_plus = deepcopy(nominal_ICs)
		ICs_plus[j]=p_plus
		ICs_minus = deepcopy(nominal_ICs)
		ICs_minus[j]=p_minus
		t_plus, A_plus = runModelWithParamsChangeICReturnA(kin_params, ICs_plus, nominal_experimental, platelets)
		t_minus, A_minus = runModelWithParamsChangeICReturnA(kin_params, ICs_minus, nominal_experimental, platelets)
		diff = (A_plus-A_minus)/2*h
		@show maximum(diff)
		S[param_idx, :]=diff
		param_idx = param_idx+1
	end

	for j in 1:maximum(size(nominal_experimental))
		@show param_idx
		if(nominal_experimental[j]==0)
			#wiggle things set to zero
			p_plus = 2.0
			p_minus = 0.0
		else
			p_plus = nominal_ICs[j]*1.5
			p_minus = nominal_ICs[j]*.5
		end
		p_plus = nominal_experimental[j]*1.5
		p_minus = nominal_experimental[j]*.5
		h = p_plus - p_minus
		exp_plus = deepcopy(nominal_experimental)
		exp_plus[j]=p_plus
		exp_minus = deepcopy(nominal_experimental)
		exp_minus[j]=p_minus
		t_plus, A_plus = runModelWithParamsChangeICReturnA(kin_params, nominal_ICs, exp_plus, platelets)
		t_minus, A_minus = runModelWithParamsChangeICReturnA(kin_params, nominal_ICs, exp_minus, platelets)
		diff = (A_plus-A_minus)/2*h
		@show maximum(diff)
		S[param_idx, :]=diff
		param_idx = param_idx+1
	end

	#and finally, platelets!
	p_plus = platelets*1.5
	p_minus = platelets*.5
	h = p_plus -p_minus
	t_plus, A_plus = runModelWithParamsChangeICReturnA(kin_params, nominal_ICs, nominal_experimental, p_plus)
	t_minus, A_minus = runModelWithParamsChangeICReturnA(kin_params, nominal_ICs, nominal_experimental, p_minus)
	diff = (A_plus-A_minus)/2*h
	S[param_idx, :]=diff
	param_idx = param_idx+1
	return S	
end

function calculateRNS_Sensitivity()
	S = calculateCentralDifferences()
	ND = size(S)[2]
	num_params = size(S)[1]
	RMSs = zeros(num_params)
	for j in 1:num_params
		sum = 0
		for k in 1:ND
			sum = sum+S[j,k]^2
		end
		RMSs[j]=sqrt(1/ND*sum) #Eqn 15 
	end
	return RMSs
end
