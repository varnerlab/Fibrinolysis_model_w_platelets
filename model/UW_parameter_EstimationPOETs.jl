using POETs
using CSV

include("runModel.jl")
include("inverseProbUtils.jl")
include("LOOCVutils.jl")

#load data once
pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndFactors.csv"
all_data = CSV.read(pathToData)

pathToMoreCompleteData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndMoreFactors.csv"
more_data = CSV.read(pathToMoreCompleteData)

function plotAllResults()
	#params = readdlm("../parameterEstimation/Best24_UW_Rescaled_16_09_19.txt")
	#ec,pc,ra =parsePOETsoutput("../parameterEstimation/POETS_info_18_09_19_PlateletContributionToROTEM_UW_Scaled.txt", 6)
	#ec,pc,ra = parsePOETsoutput("../parameterEstimation/POETS_info_18_09_19_PlateletContributionToROTEM_UW_ScaledRound2.txt",6)
	#ec, pc, ra = parsePOETsoutput("../parameterEstimation/POETS_info_19_09_19_PlateletContributionToROTEM_UW_ScaledRound3.txt",6)
	#ec, pc, ra = parsePOETsoutput("../parameterEstimation/POETS_info_20_09_19_PlateletContributionToROTEM_UW_ScaledRound4.txt",6)
	#ec, pc, ra = parsePOETsoutput("../parameterEstimation/POETS_info_20_09_19_PlateletContributionToROTEM_UW_ScaledRound5.txt",6)
	#ec, pc, ra=parsePOETsoutput( "../parameterEstimation/POETS_info_23_09_19_PlateletContributionToROTEM_UW_Scaled.txt", 6)
	#ec, pc, ra = parsePOETsoutput("../parameterEstimation/POETS_info_24_09_19_PlateletContributionToROTEM_UW_Scaled.txt",6)
	#ec,pc,ra=parsePOETsoutput("../parameterEstimation/POETS_info_25_09_19_PlateletContributionToROTEM_UW_Scaled.txt",6)
	#ec,pc,ra = parsePOETsoutput( "../parameterEstimation/POETS_info_01_10_19_PlateletContributionToROTEM_UW_Scaled.txt", 6)
	#ec,pc,ra = parsePOETsoutput("../parameterEstimation/POETS_info_02_10_19_PlateletContributionToROTEM_UW_Scaled.txt",6)
	#ec,pc,ra = parsePOETsoutput("../parameterEstimation/POETS_info_03_10_19_PlateletContributionToROTEM_UW_Scaled.txt",6)
	#ec,pc,ra = parsePOETsoutput( "../parameterEstimation/POETS_info_03_10_19_PlateletContributionToROTEM_UW_Scaled.txt",6)
	ec,pc,ra = parsePOETsoutput( "../parameterEstimation/POETS_info_07_10_19_PlateletContributionToROTEM_UW_Scaled.txt",6)
	numParamSets=2 #the number of paramter sets per objective
	params=generateNbestPerObjective(numParamSets,ec,pc)
	#params = generateBestNparameters(12, ec, pc)
	ids = [855, 856, 880]
	hours = [0,36]
	temp_params = params[1,:]
	dict = buildCompleteDictFromOneVector(temp_params)
	nominal_ICs= dict["INITIAL_CONDITION_VECTOR"]
	for h in hours
		for id in ids
			plotResult(id, h, nominal_ICs, all_data, params)
		end
	end
end

function plotAllCompleteResults()
	#to plot the results we have more protein data for
	ids = [855, 856]
	hours = [36,168]

	ec,pc,ra = parsePOETsoutput( "../parameterEstimation/POETS_info_09_10_19_PlateletContributionToROTEM_UW_Scaled_MoreComplete.txt",4)
	numParamSets=3 #the number of paramter sets per objective
	params=generateNbestPerObjective(numParamSets,ec,pc)
	temp_params =params[1,:]
	dict = buildCompleteDictFromOneVector(temp_params)
	nominal_ICs= dict["INITIAL_CONDITION_VECTOR"]
	for h in hours
		for id in ids
			plotResultComplete(id, h, nominal_ICs, more_data, params)
		end
	end
end

function plotResultComplete(id::Int, hour::Int, nominal_ICs, all_data, params)
	close("all")
	sel_target  = constructTarget(hour, id, all_data)
	figure(figsize =[15,15])
	num_param_sets = size(params,1)
	TSTART =0.0
	TSTOP = 60.0
	step =.1
	curr_platelets=300.0
	@show id, hour
	@show sel_target
	all_trajectories = zeros(num_param_sets, maximum(size(collect(TSTART:step:TSTOP))))
	for j in 1:num_param_sets
		curr_params = params[j,:]
		dict = buildCompleteDictFromOneVector(curr_params)
		initial_parameter_estimate = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector,dict = updateICsForCompletePatient(id, hour,more_data, initial_parameter_estimate,dict)
		#update our time delay
		est_delay =sampleTimeDelay_UW_R()
		dict["TIME_DELAY"][1]=est_delay
		fbalances(y,p,t)= Balances(t,y,dict) 
		prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
		sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = step,maxiters = 1e6)
		t =sol.t
		X = sol
		#A = convertToROTEM(t,X,tPA)
		tPA = initial_condition_vector[16]
		A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
		all_trajectories[j,:]=A
		plot(t,A, "grey", alpha = .5)
		metrics = calculateCommonMetrics(A,t)
		CT = metrics[1]
		CFT = metrics[2]
		alpha = metrics[3]
		MCF = metrics[4]
		LI30 = calculateLysisAtTime(A,t,30.0)
		sel_metrics = [CT,CFT,alpha,MCF,LI30]
		#@show broadcast(abs, sel_metrics-sel_target)
	end
	#plot the results of using the mean paramters
	meanp = mean(params, dims=1)
	curr_params = meanp
	dict = buildCompleteDictFromOneVector(curr_params)
	initial_parameter_estimate = dict["INITIAL_CONDITION_VECTOR"]
	initial_condition_vector,dict = updateICsForCompletePatient(id, hour,more_data, initial_parameter_estimate,dict)
	#update our time delay
	est_delay =sampleTimeDelay_UW_R()
	dict["TIME_DELAY"][1]=est_delay
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = step,maxiters = 1e6)
	t =sol.t
	X = sol
	#A = convertToROTEM(t,X,tPA)
	tPA = initial_condition_vector[16]
	A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
	plot(t,A, linestyle = "-.")

	meanTrajectory = vec(mean(all_trajectories, dims=1))
	stdTrajectory = vec(std(all_trajectories, dims=1))
	TSIM = collect(TSTART:step:TSTOP)
	plot(TSIM, meanTrajectory, "k")
	upper = vec(meanTrajectory +stdTrajectory)
	lower = vec(meanTrajectory -stdTrajectory)
	fill_between((TSIM), vec(upper), vec(lower), color = "lightskyblue", alpha =.5)

	xlabel("Time (minutes)", fontsize = 20)
	ylabel("Clot Amplitude (mm)", fontsize = 20)
	titlestr =string("Patient ", id, " Hour ", hour)
	title(titlestr)
	#let's plot ideal metrics on top, in red
	plot(sel_target[1]/60, 2, "rx", markersize = 20) #CT
	plot(sel_target[1]/60+sel_target[2]/60, 20, "rx", markersize = 20) #CFT
	plot([sel_target[1]/60, sel_target[1]/60+sel_target[2]/60], [0,20], color = "red", linestyle = "dotted") #alpha
	plot([20,40], [sel_target[4], sel_target[4]], color = "red", linestyle ="dashed") #MCF
	#@show 1-sel_target[5]/100
	plot(30, sel_target[4]*(1-sel_target[5]/100), "rx", markersize=20) #LI30

	savefig(string("../figures/UW_ParamEst/Complete",titlestr,".pdf"))

end

function plotResult(id::Int, hour::Int, nominal_ICs, all_data, params)
	close("all")
	initial_condition_vector = updateICsForPatient(id, hour, all_data, nominal_ICs)
	sel_target  = constructTarget(hour, id, all_data)
	figure(figsize =[15,15])
	num_param_sets = size(params,1)
	TSTART =0.0
	TSTOP = 60.0
	step =.1
	curr_platelets=300.0
	@show id, hour
	@show sel_target
	all_trajectories = zeros(num_param_sets, maximum(size(collect(TSTART:step:TSTOP))))
	for j in 1:num_param_sets
		curr_params = params[j,:]
		dict = buildCompleteDictFromOneVector(curr_params)
		#update our time delay
		est_delay =sampleTimeDelay_UW_R()
		dict["TIME_DELAY"][1]=est_delay
		fbalances(y,p,t)= Balances(t,y,dict) 
		prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
		sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = step,maxiters = 1e6)
		t =sol.t
		X = sol
		#A = convertToROTEM(t,X,tPA)
		tPA = initial_condition_vector[16]
		A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
		all_trajectories[j,:]=A
		plot(t,A, "grey", alpha = .5)
		metrics = calculateCommonMetrics(A,t)
		CT = metrics[1]
		CFT = metrics[2]
		alpha = metrics[3]
		MCF = metrics[4]
		LI30 = calculateLysisAtTime(A,t,30.0)
		sel_metrics = [CT,CFT,alpha,MCF,LI30]
		#@show broadcast(abs, sel_metrics-sel_target)
	end
	meanTrajectory = vec(mean(all_trajectories, dims=1))
	stdTrajectory = vec(std(all_trajectories, dims=1))
	TSIM = collect(TSTART:step:TSTOP)
	plot(TSIM, meanTrajectory, "k")
	upper = vec(meanTrajectory +stdTrajectory)
	lower = vec(meanTrajectory -stdTrajectory)
	fill_between((TSIM), vec(upper), vec(lower), color = "lightskyblue", alpha =.5)

	xlabel("Time (minutes)", fontsize = 20)
	ylabel("Clot Amplitude (mm)", fontsize = 20)
	titlestr =string("Patient ", id, " Hour ", hour)
	title(titlestr)
	#let's plot ideal metrics on top, in red
	plot(sel_target[1]/60, 2, "rx", markersize = 20) #CT
	plot(sel_target[1]/60+sel_target[2]/60, 20, "rx", markersize = 20) #CFT
	plot([sel_target[1]/60, sel_target[1]/60+sel_target[2]/60], [0,20], color = "red", linestyle = "dotted") #alpha
	plot([20,40], [sel_target[4], sel_target[4]], color = "red", linestyle ="dashed") #MCF
	#@show 1-sel_target[5]/100
	plot(30, sel_target[4]*(1-sel_target[5]/100), "rx", markersize=20) #LI30

	savefig(string("../figures/UW_ParamEst/HeavierWeightedCTandCFT",titlestr,".pdf"))

end

function plotAllResults1Fig(params)
	close("all")
	temp_params = params[1,:]
	dict = buildCompleteDictFromOneVector(temp_params)
	nominal_ICs= dict["INITIAL_CONDITION_VECTOR"]
	#fig =figure(figsize =[15,15])
	fig, all_ax = subplots(nrows =3, ncols = 2, sharex = true, figsize = [15,15])
	num_param_sets = size(params,1)
	ids = [855, 856, 880]
	hours = [0,36]
	TSTART =0.0
	TSTOP = 60.0
	step =.1
	curr_platelets=300.0
	count = 1
	for (i, hour) in enumerate(hours)
		for (k,id) in enumerate(ids)
			all_trajectories = zeros(num_param_sets, maximum(size(collect(TSTART:step:TSTOP))))
			initial_condition_vector = updateICsForPatient(id, hour, all_data, nominal_ICs)
			@show hour, id, initial_condition_vector
			sel_target  = constructTarget(hour, id, all_data)
			#ax = fig.add_subplot(3,2,count, sharex = true)
			ax = all_ax[k,i]
			for j in 1:num_param_sets
				curr_params = params[j,:]
				dict = buildCompleteDictFromOneVector(curr_params)
				#update our time delay
				est_delay =sampleTimeDelay_UW_R()
				dict["TIME_DELAY"][1]=est_delay
				fbalances(y,p,t)= Balances(t,y,dict) 
				prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
				sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = step,maxiters = 1e6)
				t =sol.t
				X = sol
				#A = convertToROTEM(t,X,tPA)
				tPA = initial_condition_vector[16]
				A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
				all_trajectories[j,:]=A
				ax.plot(t,A, "grey", alpha = .5)
				metrics = calculateCommonMetrics(A,t)
				CT = metrics[1]
				CFT = metrics[2]
				alpha = metrics[3]
				MCF = metrics[4]
				LI30 = calculateLysisAtTime(A,t,30.0)
				sel_metrics = [CT,CFT,alpha,MCF,LI30]
				#@show broadcast(abs, sel_metrics-sel_target)
			end
			meanTrajectory = vec(mean(all_trajectories, dims=1))
			stdTrajectory = vec(std(all_trajectories, dims=1))
			TSIM = collect(TSTART:step:TSTOP)
			ax.plot(TSIM, meanTrajectory, "k")
			upper = vec(meanTrajectory +stdTrajectory)
			lower = vec(meanTrajectory -stdTrajectory)
			ax.fill_between((TSIM), vec(upper), vec(lower), color = "lightskyblue", alpha =.5)

			#ax.xlabel("Time (minutes)", fontsize = 20)
			#ax.ylabel("Clot Amplitude (mm)", fontsize = 20)
			titlestr =string("Patient ", id, " Hour ", hour)
			ax.set_title(titlestr)
			#let's plot ideal metrics on top, in red
			ax.plot(sel_target[1]/60, 2, "rx", markersize = 20) #CT
			ax.plot(sel_target[1]/60+sel_target[2]/60, 20, "rx", markersize = 20) #CFT
			ax.plot([sel_target[1]/60, sel_target[1]/60+sel_target[2]/60], [0,20], color = "red", linestyle = "dotted") #alpha
			ax.plot([20,40], [sel_target[4], sel_target[4]], color = "red", linestyle ="dashed") #MCF
			#@show 1-sel_target[5]/100
			ax.plot(30, sel_target[4]*(1-sel_target[5]/100), "rx", markersize=20) #LI30
			count  =count+1
		end
	end
	savefig(string("../figures/UW_ParamEst/CompositeFig.pdf"))
end



function updateICsForPatient(id::Int, hour::Int, all_data, initial_parameter_estimate)
	curr_dat = all_data[.&(all_data.id .== id, all_data.hour.==hour), :]
	nominal_ICs = initial_parameter_estimate
	offset =1
	names = ["FII", "PC", "ATIII", "Fibrinogen", "plasminogen", "FV_FX"]

	#all in nM
	nominal_II = nominal_ICs[1]
	nominal_PC =nominal_ICs[3]
	nominal_ATIII = nominal_ICs[5]
	nominal_FV_FX = nominal_ICs[9]
	nominal_Fibrinogen = nominal_ICs[14]
	nominal_plasminogen = nominal_ICs[15]

	actual_FII = nominal_II*curr_dat.II[1]/100
	actual_PC = nominal_PC*curr_dat.ProtC[1]/100
	actual_ATIII = nominal_ATIII*curr_dat.AT[1]/100
	actual_Fibrinogen = convertFibrinogenTonM(curr_dat.Fbgn[1])
	actual_plasminogen = nominal_plasminogen*curr_dat.Plgn[1]/100
	#take min of the two
	smaller =minimum([curr_dat.X[1], curr_dat.V[1]])
	actual_FV_FX = nominal_FV_FX*smaller/100

	#update intiial condition vector
	initial_parameter_estimate[1]=actual_FII
	initial_parameter_estimate[3]=actual_PC
	initial_parameter_estimate[5]=actual_ATIII
	initial_parameter_estimate[9]=actual_FV_FX
	initial_parameter_estimate[14]=actual_Fibrinogen
	initial_parameter_estimate[15]=actual_plasminogen
	return initial_parameter_estimate
	
end

function updateICsForCompletePatient(id::Int, hour::Int, all_data, initial_parameter_estimate,dict)
	curr_dat = all_data[.&(all_data.id .== id, all_data.hour.==hour), :]
	nominal_ICs = initial_parameter_estimate
	offset =1
	names = ["FII", "PC", "ATIII", "Fibrinogen", "plasminogen", "FV_FX", "PAI-1", "TFPI", "TAFI"]

	#all in nM
	nominal_II = nominal_ICs[1]
	nominal_PC =nominal_ICs[3]
	nominal_ATIII = nominal_ICs[5]
	nominal_FV_FX = nominal_ICs[9]
	nominal_Fibrinogen = nominal_ICs[14]
	nominal_plasminogen = nominal_ICs[15]
	nominal_PAI1=nominal_ICs[21]

	actual_FII = nominal_II*curr_dat.II[1]/100
	actual_PC = nominal_PC*curr_dat.ProtC[1]/100
	actual_ATIII = nominal_ATIII*curr_dat.AT[1]/100
	actual_Fibrinogen = convertFibrinogenTonM(curr_dat.Fbgn[1])
	actual_plasminogen = nominal_plasminogen*curr_dat.Plgn[1]/100
	#take min of the two
	smaller =minimum([curr_dat.X[1], curr_dat.V[1]])
	actual_FV_FX = nominal_FV_FX*smaller/100
	actual_PAI1 = convertPAIToM(curr_dat.PAI1[1])

	#update intiial condition vector
	initial_parameter_estimate[1]=actual_FII
	initial_parameter_estimate[3]=actual_PC
	initial_parameter_estimate[5]=actual_ATIII
	initial_parameter_estimate[9]=actual_FV_FX
	initial_parameter_estimate[14]=actual_Fibrinogen
	initial_parameter_estimate[15]=actual_plasminogen
	initial_parameter_estimate[21]=actual_PAI1
	
	dict["FACTOR_LEVEL_VECTOR"][1]=convertTFPIToM(curr_dat.TFPIFree[1])
	dict["FACTOR_LEVEL_VECTOR"][7]=convertTAFIToM(curr_dat.TAFIaai[1])

	return initial_parameter_estimate, dict
end

function constructTarget(hour::Int, patientID::Int, all_data)
	curr_dat = all_data[.&(all_data.id .== patientID, all_data.hour.==hour), :]
	target_CT =  curr_dat.teg_rapidteg_r[1]*60
	target_CFT = curr_dat.teg_rapidteg_[1]*60 #k is missing because of find and replace used to clean up data
	target_alpha = curr_dat.teg_rapidteg_angle[1]
	target_MCF = curr_dat.teg_rapidteg_ma[1]
	target_LI30 = curr_dat.teg_rapidteg_ly30[1]
	
	target = [target_CT, target_CFT, target_alpha, target_MCF,target_LI30]
	return target
end

function objectiveUWMoreComplete(parameter_array)
	obj_array =ones(4,1)
	ids = [855, 856]
	TSTART = 0.0
	TSTOP =60.0
	#since we don't actually know platelets, give a normal platelet count'
	curr_platelets =300
	#use the median as characteritic scaling
	scale_CT = 6.7*60
	scale_CFT = 2*60
	scale_alpha = 62.3
	scale_MCF = 60.6
	scale_MaxLysis = .12
	scale_AUC = 250.0
	hours = [168,36]

	sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]
	#weight all metrics equally for now
	weights = ones(size(sel_scales))
	#heavily weight CT and CFT
	weights[1]=100.0*10
	weights[2]=100.0*10
	weights[4]=1000.0*100
	count = 1
	threshold = 90.0
	for (idx,id) in enumerate(ids)
		for hour in hours
			#@show hour, id, count
			temp_params = parameter_array
			sel_target  = constructTarget(hour, id, all_data)
			dict = buildCompleteDictFromOneVector(temp_params)
			#update our time delay
			est_delay =sampleTimeDelay_UW_R()
			#@show est_delay
			dict["TIME_DELAY"][1]=est_delay
			initial_parameter_estimate = dict["INITIAL_CONDITION_VECTOR"]
			initial_condition_vector,dict = updateICsForCompletePatient(id, hour,more_data, initial_parameter_estimate,dict)
			#now use our actually measured ICs to update initial_parameter_estimate
			fbalances(y,p,t)= Balances(t,y,dict) 
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = .1,maxiters = 1e6)
			t =sol.t
			X = sol
			#A = convertToROTEM(t,X,tPA)
			tPA = initial_condition_vector[16]
			A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
			#plot(t,A)
			metrics = calculateCommonMetrics(A,t)
			CT = metrics[1]
			CFT = metrics[2]
			alpha = metrics[3]
			MCF = metrics[4]
			LI30 = calculateLysisAtTime(A,t,30.0)
			sel_metrics = [CT,CFT,alpha,MCF,LI30]
			#@show broadcast(abs, (sel_metrics-sel_target))
			#if any of the metrics are negative, something went wrong, penalize our parameters
			if(true in (sel_metrics .<0)) #if any of our metrics are negative
				calc_obj = 1E8
			elseif(MCF<30.0) #if we're not getting a large enough MCF'

				calc_obj = 1E8

			elseif(abs(CT-sel_metrics[1])>threshold || abs(CFT-sel_metrics[2])>threshold)
				#if our CT and CFT are too far from where they should be
				calc_obj=1E10
			else
				sum = 0.0
				#use MSE to determine how far away we are from hitting our target
				for j in 1:maximum(size(sel_metrics))
					sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2*weights[j]
				end
				#@show params
				calc_obj = sqrt(sum)
			end
			obj_array[count]=calc_obj
			count = count+1
		end
	end
	@show obj_array
	return obj_array

end

function objectiveUWHour0(parameter_array)
	obj_array =ones(6,1)
	ids = [855, 856,880]
	TSTART = 0.0
	TSTOP =60.0
	#since we don't actually know platelets, give a normal platelet count'
	curr_platelets =300
	#use the median as characteritic scaling
	scale_CT = 6.7*60
	scale_CFT = 2*60
	scale_alpha = 62.3
	scale_MCF = 60.6
	scale_MaxLysis = .12
	scale_AUC = 250.0
	hours = [0,36]

	sel_scales = [scale_CT, scale_MCF, scale_alpha, scale_MCF, scale_MaxLysis]
	#weight all metrics equally for now
	weights = ones(size(sel_scales))
	#heavily weight CT and CFT
	weights[1]=100.0*10
	weights[2]=100.0*10
	weights[4]=1000.0*100
	count = 1
	threshold = 90.0
	for (idx,id) in enumerate(ids)
		for hour in hours
			#@show hour, id, count
			temp_params = parameter_array
			sel_target  = constructTarget(hour, id, all_data)
			dict = buildCompleteDictFromOneVector(temp_params)
			#update our time delay
			est_delay =sampleTimeDelay_UW_R()
			@show est_delay
			dict["TIME_DELAY"][1]=est_delay
			initial_parameter_estimate = dict["INITIAL_CONDITION_VECTOR"]
			initial_condition_vector = updateICsForPatient(id, hour,all_data, initial_parameter_estimate)
			#now use our actually measured ICs to update initial_parameter_estimate
			fbalances(y,p,t)= Balances(t,y,dict) 
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = .1,maxiters = 1e6)
			t =sol.t
			X = sol
			#A = convertToROTEM(t,X,tPA)
			tPA = initial_condition_vector[16]
			A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
			#plot(t,A)
			metrics = calculateCommonMetrics(A,t)
			CT = metrics[1]
			CFT = metrics[2]
			alpha = metrics[3]
			MCF = metrics[4]
			LI30 = calculateLysisAtTime(A,t,30.0)
			sel_metrics = [CT,CFT,alpha,MCF,LI30]
			@show broadcast(abs, (sel_metrics-sel_target))
			#if any of the metrics are negative, something went wrong, penalize our parameters
			if(true in (sel_metrics .<0)) #if any of our metrics are negative
				calc_obj = 1E8
			elseif(MCF<30.0) #if we're not getting a large enough MCF'

				calc_obj = 1E8

			elseif(abs(CT-sel_metrics[1])>threshold || abs(CFT-sel_metrics[2])>threshold)
				#if our CT and CFT are too far from where they should be
				calc_obj=1E10
			else
				sum = 0.0
				#use MSE to determine how far away we are from hitting our target
				for j in 1:maximum(size(sel_metrics))
					sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2*weights[j]
				end
				#@show params
				calc_obj = sqrt(sum)
			end
			obj_array[count]=calc_obj
			count = count+1
		end
	end
	@show obj_array
	return obj_array

end

function attemptOptimization_UW()
	number_of_subdivisions = 40
	number_of_parameters = 77
	number_of_objectives = 6
	outputfile = "../parameterEstimation/POETS_info_07_10_19_PlateletContributionToROTEM_UW_Scaled.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	#allp = readdlm("../parameterEstimation/Best12_UW_Rescaled_19_09_19.txt")
	#allp = readdlm("../parameterEstimation/Best12_UW_Rescaled_20_09_19.txt")
	#allp = readdlm("../parameterEstimation/Best12_UW_Rescaled_20_09_19Round4.txt")
	#allp= readdlm("../parameterEstimation/Best12_UW_Rescaled_09_23_19.txt")
	#allp= readdlm("../parameterEstimation/Best12_UW_Rescaled_02_10_19.txt")
	#allp = readdlm("../parameterEstimation/Best12_UW_Rescaled_03_10_19.txt")
	initial_parameter_estimate=vec(mean(allp, dims=1))

	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*4, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/4, initial_parameter_estimate[47:end]/1000)
	EC = []
	PC = []
	RA = []
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		if(index==1)
			initial_parameter_estimate =initial_parameter_estimate+initial_parameter_estimate*rand()*.1
		else
			#use results from previous round to start
			sum_errs = sum(EC, dims=1)
			min_idx =  findmin(sum_errs)[2][2] #weird indexing because it returns a cartesian idx
			bestParams = PC[:, min_idx]
			initial_parameter_estimate = bestParams

		end

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble(objectiveUWHour0,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)

		# Package -
		@show (EC, PC, RA)
		ec_array = [ec_array EC]
		pc_array = [pc_array PC]
		f = open(outputfile, "a")
		write(f, string(EC, ",", PC, ",", RA, "\n"))
		close(f)
	end

	return (ec_array,pc_array)
end

function attemptOptimization_UW_MoreComplete()
	number_of_subdivisions = 40
	number_of_parameters = 77
	number_of_objectives = 4
	outputfile = "../parameterEstimation/POETS_info_09_10_19_PlateletContributionToROTEM_UW_Scaled_MoreComplete.txt"
	ec_array = zeros(number_of_objectives)
	pc_array = zeros(number_of_parameters)
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	initial_parameter_estimate=vec(mean(allp, dims=1))

	#bound thrombin generation parameters more tightly than fibrinolysis ones
	global up_arr = vcat(initial_parameter_estimate[1:46]*4, initial_parameter_estimate[47:end]*1000)
	global lb_arr = vcat(initial_parameter_estimate[1:46]/4, initial_parameter_estimate[47:end]/1000)
	EC = []
	PC = []
	RA = []
	for index in collect(1:number_of_subdivisions)

		# Grab a starting point -
		if(index==1)
			initial_parameter_estimate =initial_parameter_estimate+initial_parameter_estimate*rand()*.1
		else
			#use results from previous round to start
			sum_errs = sum(EC, dims=1)
			min_idx =  findmin(sum_errs)[2][2] #weird indexing because it returns a cartesian idx
			bestParams = PC[:, min_idx]
			initial_parameter_estimate = bestParams

		end

		# Run JuPOETs -
		(EC,PC,RA) = estimate_ensemble( objectiveUWMoreComplete,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_estimate;rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)

		# Package -
		@show (EC, PC, RA)
		ec_array = [ec_array EC]
		pc_array = [pc_array PC]
		f = open(outputfile, "a")
		write(f, string(EC, ",", PC, ",", RA, "\n"))
		close(f)
	end

	return (ec_array,pc_array)
end
