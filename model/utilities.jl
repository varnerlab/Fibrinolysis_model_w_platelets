using PyPlot
using ExcelReaders

@everywhere function calculateMSE(t,predictedCurve, experimentalData)
	num_points = size(t,1)
	interpolatedExperimentalData = Float64[]
	for j in collect(1:num_points)
		currt = t[j]
		upperindex = searchsortedfirst(experimentalData[:,1], currt, by=abs)
		if(upperindex ==1)
			lowerindex = 1
		else
			lowerindex = upperindex -1
		end

		if(upperindex >= size(experimentalData,1)&& upperindex !=1)
			upperindex = size(experimentalData,1)
			lowerindex = upperindex-1
		end
		val = linearInterp(experimentalData[lowerindex,2], experimentalData[upperindex,2],experimentalData[lowerindex,1], experimentalData[upperindex,1], currt)
		#@show currt, val,experimentalData[lowerindex,2], experimentalData[upperindex,2]
		push!(interpolatedExperimentalData,val)
	end
	sum = 0.0
	for j in collect(1:size(predictedCurve,1))
		sum = sum +(predictedCurve[j]-interpolatedExperimentalData[j])^2
	end
#	figure()
#	plot(t, interpolatedExperimentalData)
	@show sum, size(predictedCurve,1)
	return sum/size(predictedCurve,1), interpolatedExperimentalData#MSE
end

@everywhere function linearInterp(lowerVal, upperVal, tstart, tend,tdesired)
	
	val = lowerVal + (upperVal-lowerVal)/(tend-tstart)*(tdesired-tstart)
	if(isnan(val))
		val = 0.0
	end
	return val
end

@everywhere function calculateAUC(t,y)
	   local n = length(t)
    if (length(y) != n)
        error("Vectors 't', 'y' must be of same length")
    end

	sum = 0.0
	for j in collect(2:n)
		curr_rect = (t[j]-t[j-1])*(y[j]+y[j-1])/2
		sum = sum+curr_rect
	end
	return sum
end

@everywhere function buildDictFromOneVector(vector)
	kinetic_parameter_vector = vector[1:18]
	control_parameter_vector=vector[19:38]
	platelet_parameter_vector=vector[39:44]
	timing = vector[45:46]
	platelet_count=vector[47]
	dict = buildCoagulationModelDictionary(kinetic_parameter_vector, control_parameter_vector, platelet_parameter_vector, timing, platelet_count)
	return dict
end

@everywhere function buildCompleteDictFromOneVector(vector)
	kinetic_parameter_vector = vector[1:18]
	control_parameter_vector=vector[19:38]
	platelet_parameter_vector=vector[39:44]
	timing = vector[45:46]
	platelet_count=vector[47]
	fibrin_kinetic_parameters = vector[48:69]
	fibrin_control_parameters = vector[70:77]
	dict = buildCoagulationModelDictionary(kinetic_parameter_vector, control_parameter_vector, platelet_parameter_vector, timing, platelet_count, fibrin_kinetic_parameters, fibrin_control_parameters)
	return dict
end




# Generates new parameter array, given current array -
@everywhere function neighbor_function(parameter_array)
	outputfile="parameterEstimation/POETS_28_03_2017.txt"
	#@show size(parameter_array)
#	f = open(outputfile, "a")
#	write(f,string(parameter_array, "\n"))
#	close(f)

  SIGMA = 0.05
  number_of_parameters = length(parameter_array)

  # calculate new parameters -
  new_parameter_array = parameter_array.*(1+SIGMA*randn(number_of_parameters))

  # Check the bound constraints -
  LOWER_BOUND = 0
  UPPER_BOUND = 1E9
#  lb_arr= LOWER_BOUND*ones(number_of_parameters)
#  up_arr =UPPER_BOUND*ones(number_of_parameters)
	#lb_arr[9] = 10.0 #lower bound on k_inhibition_ATIII
	#lb_arr[45]= 3.0 #lower bound on time delay, 3 minutes
	#up_arr[46]= .01 #upper bound on scaling for tau
	#up_arr[3] = 10.0 #upper bound on the k_cat for self activation of thrombin

  # return the corrected parameter arrays -
  return parameter_bounds_function(new_parameter_array,lb_arr, up_arr)
end

@everywhere function acceptance_probability_function(rank_array,temperature)
  return (exp(-rank_array[end]/temperature))
end

@everywhere function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
	@show temperature
  return alpha*temperature
end


# Helper functions -
@everywhere function parameter_bounds_function(parameter_array,lower_bound_array,upper_bound_array)

  # reflection_factor -
  epsilon = 0.01

  # iterate through and fix the parameters -
  new_parameter_array = copy(parameter_array)
  for (index,value) in enumerate(parameter_array)

    lower_bound = lower_bound_array[index]
    upper_bound = upper_bound_array[index]

    if (value<lower_bound)
      new_parameter_array[index] = lower_bound
    elseif (value>upper_bound)
      new_parameter_array[index] = upper_bound
    end
  end

  return new_parameter_array
end

@everywhere function createCorrectDict(basic_dict, exp_index)
	if(exp_index==1)
		
	elseif(exp_index==2)
		basic_dict["FACTOR_LEVEL_VECTOR"][3] =basic_dict["FACTOR_LEVEL_VECTOR"][3]*1.06 
	elseif(exp_index==3)
		basic_dict["FACTOR_LEVEL_VECTOR"][3] =basic_dict["FACTOR_LEVEL_VECTOR"][3]*.39
	elseif(exp_index==4)
		basic_dict["FACTOR_LEVEL_VECTOR"][3] =basic_dict["FACTOR_LEVEL_VECTOR"][3]*.07  
	elseif(exp_index==5)
		basic_dict["FACTOR_LEVEL_VECTOR"][3] =basic_dict["FACTOR_LEVEL_VECTOR"][3]*.01
	elseif(exp_index==6)
		basic_dict["FACTOR_LEVEL_VECTOR"][3] =basic_dict["FACTOR_LEVEL_VECTOR"][3]*.00
	end
	return basic_dict
end

@everywhere function generateBestNparameters(n, ec_array, pc_array)
	#calculate error
	best_params = Array[]
	total_error = sum(ec_array[:,2:end],1)
	total_error= vec(total_error)
	for k in collect(1:n)
		min_index = indmin(total_error)
		curr_best_params = pc_array[:,min_index]
		@show size(curr_best_params)
		push!(best_params, curr_best_params)
		@show min_index
		@show curr_best_params
		#delete the best ones we've found
		@show size(pc_array)
		@show size(total_error)
		pc_array[1:size(pc_array,1) .!= min_index,: ]
		total_error=deleteat!((total_error),min_index)
		@show size(pc_array)
		@show size(total_error)
	end
	writedlm(string("parameterEstimation/Best", n, "OverallParameters_19_04_2017.txt"), best_params)
	return best_params

end

@everywhere function generateNbestPerObjective(n,ec_array, pc_array)
	num_objectives =size(ec_array,1)
	best_params=Array{Array}(num_objectives*n)
	counter = 1
	for j in collect(1:num_objectives)
		curr_error=ec_array[j,:]
		allidx = collect(1:size(pc,2))
		removed = allidx
		for k in collect(1:n)
			min_index = indmin(curr_error)
			curr_best_params = pc_array[:,min_index]
			best_params[counter] = curr_best_params
			removed = deleteat!(removed, min_index)
			@show min_index
			#@show curr_best_params
			#delete the best ones we've found
			#@show size(total_error)
			#pc_array=pc_array[:,removed]
			curr_error=deleteat!(vec(curr_error),min_index)
			@show size(pc_array)
			counter=counter+1
		end
	end	
	writedlm(string("parameterEstimation/Best", n, "PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"), best_params)
	return best_params
end

@everywhere function generateNbestGivenObjective(n,ec_array, pc_array,objectivenum)
	num_objectives =size(ec_array,1)
	best_params=Array{Array}(n)
	counter = 1
	curr_error=ec_array[objectivenum,:]
	allidx = collect(1:size(pc_array,2))
	removed = allidx
	for k in collect(1:n)
		min_index = indmin(curr_error)
		curr_best_params = pc_array[:,min_index]
		best_params[counter] = curr_best_params
		removed = deleteat!(removed, min_index)
		@show min_index, curr_error[min_index]
		#@show curr_best_params
		#delete the best ones we've found
		#@show size(total_error)
		#pc_array=pc_array[:,removed]
		curr_error=deleteat!(vec(curr_error),min_index)
		@show size(pc_array)
		counter=counter+1
	end
	return best_params
end


function analyzeParams()
	allparams = zeros(1,46)
	for j in collect(1:6)
		currparams = readdlm(string("parameterEstimation/LOOCVSavingAllParams_2016_12_23_Take2/bestParamSetsFromLOOCV",j,"excluded.txt"), ',')
		meancurrparams = mean(currparams,1)
		allparams = vcat(allparams, meancurrparams)
	end
	return allparams[2:end,:]
end

function dict_to_vec(d)
    v = Array(Float64, 0)
	selectedkeys = ["FACTOR_LEVEL_VECTOR","CONTROL_PARAMETER_VECTOR","PLATELET_PARAMS", "TIME_DELAY","KINETIC_PARAMETER_VECTOR", "ALEPH"] 
    for k in selectedkeys
	for j in collect(1:length(d[k]))
        	push!(v, d[k][j])
	end
    end
    return v
end

function extractValueFromDual(input)
	# to get the value of a dual as Float64 that can be used by min/max
	if(contains(string(typeof(input)), "Dual"))
		return input.value
	else
		return input
	end
end

function build_param_dict(problem_vec)
	params = Expr[:(p_1=>$(problem_vec[1])); :(p_2=>$(problem_vec[2])); :(p_3=>$(problem_vec[3])); :(p_4=>$(problem_vec[4])); :(p_5=>$(problem_vec[5])); :(p_6=>$(problem_vec[6])); :(p_7=>$(problem_vec[7])); :(p_8=>$(problem_vec[8])); :(p_9=>$(problem_vec[9])); :(p_10=>$(problem_vec[10])); :(p_11=>$(problem_vec[11])); :(p_12=>$(problem_vec[12])); :(p_13=>$(problem_vec[13])); :(p_14=>$(problem_vec[14])); :(p_15=>$(problem_vec[15])); :(p_16=>$(problem_vec[16])); :(p_17=>$(problem_vec[17])); :(p_18=>$(problem_vec[18])); :(p_19=>$(problem_vec[19])); :(p_20=>$(problem_vec[20])); :(p_21=>$(problem_vec[21])); :(p_22=>$(problem_vec[22])); :(p_23=>$(problem_vec[23])); :(p_24=>$(problem_vec[24])); :(p_25=>$(problem_vec[25])); :(p_26=>$(problem_vec[26])); :(p_27=>$(problem_vec[27])); :(p_28=>$(problem_vec[28])); :(p_29=>$(problem_vec[29])); :(p_30=>$(problem_vec[30])); :(p_31=>$(problem_vec[31])); :(p_32=>$(problem_vec[32])); :(p_33=>$(problem_vec[33])); :(p_34=>$(problem_vec[34])); :(p_35=>$(problem_vec[35])); :(p_36=>$(problem_vec[36])); :(p_37=>$(problem_vec[37])); :(p_38=>$(problem_vec[38])); :(p_39=>$(problem_vec[39])); :(p_40=>$(problem_vec[40])); :(p_41=>$(problem_vec[41])); :(p_42=>$(problem_vec[42])); :(p_43=>$(problem_vec[43])); :(p_44=>$(problem_vec[44])); :(p_45=>$(problem_vec[45])); :(p_46=>$(problem_vec[46])); :(p_47=>$(problem_vec[47])); :(p_48=>$(problem_vec[48])); :(p_49=>$(problem_vec[49])); :(p_50=>$(problem_vec[50])); :(p_51=>$(problem_vec[51])); :(p_52=>$(problem_vec[52])); :(p_53=>$(problem_vec[53]))]
	return params
end

function parsePOETsoutput(filename)
	f = open(filename)
	alltext = readall(f)
	close(f)

	outputname = "textparsing.txt"
	number_of_parameters = 77
  	number_of_objectives = 8#4
	ec_array = zeros(number_of_objectives)
  	pc_array = zeros(number_of_parameters)
	rank_array = zeros(1)	
	counter =1
	for grouping in matchall(r"\[([^]]+)\]", alltext)
		cleanedgrouping = replace(grouping, "[", "")
		nocommas = replace(cleanedgrouping, ","," ")
		allcleaned = replace(nocommas, "]", "")
		allcleaned = replace(allcleaned, ";", "\n")
		outfile = open(outputname, "w")
		write(outfile, allcleaned)
		close(outfile)
		formatted = readdlm(outputname)
		@show formatted	
		@show size(formatted), counter
		if(counter == 1)
			ec_array = [ec_array formatted]
			counter = counter +1
		elseif(counter == 2)
			pc_array = [pc_array formatted]
			counter = counter +1
		elseif(counter == 3)
			rank_array = [rank_array formatted]
			#@show formatted
			counter =1
		end
		
	end
	return ec_array[:,2:end], pc_array[:,2:end], rank_array[:,2:end]
end

function plotTradeOffCurve(ec_array, rank_array)
	figure()
	PyCall.PyDict(matplotlib["rcParams"])["font.sans-serif"] = ["Helvetica"]
	hold("on")
	#@show size(ec_array,2)
	for j in collect(1:size(ec_array,2))
		if(rank_array[j]==0)
			plot(ec_array[1,j], ec_array[2,j], linewidth = .3,"ko", markersize = 2.5,markeredgewidth=0.0)
		else
			plot(ec_array[1,j], ec_array[2,j], linewidth = .3,"o", color = ".75", markersize = 2.5,markeredgewidth=0.0)
		end
	end
	xlabel("Objective 1", fontsize=18)
	ylabel("Objective 2", fontsize=18)
	axis([0,4000,0,4000])
	savefig("figures/tradeoffCurve_07_04_2017.pdf")
end

function plotTradeOffCurve(ec_array, rank_array, obj1, obj2)
	close("all")
	figure()
	PyCall.PyDict(matplotlib["rcParams"])["font.sans-serif"] = ["Helvetica"]
	hold("on")
	#@show size(ec_array,2)
	for j in collect(1:size(ec_array,2))
		if(rank_array[j]==0)
			plot(ec_array[obj1,j], ec_array[obj2,j], linewidth = .3,"ko", markersize = 2.5,markeredgewidth=0.0)
		else
			plot(ec_array[obj1,j], ec_array[obj2,j], linewidth = .3,"o", color = ".75", markersize = 2.5,markeredgewidth=0.0)
		end
	end
	xlabel(string("Objective ", obj1), fontsize=18)
	ylabel(string("Objective ", obj2), fontsize=18)
	axis([0,4000,0,4000])
	savefig(string("figures/TradeOffCurves/tradeoffCurve_19_05_2017_obj_",obj1,"and_obj_",obj2, ".pdf"))
end

function plotAllTradeOffCurves(ec_array, ra_array, num_objectives)
	idx1 = collect(1:1:num_objectives)
	idx2= collect(1:1:num_objectives)
	for(j in idx1)
		for(k in idx2)
			if(j!=k)
				plotTradeOffCurve(ec_array, ra_array, j, k)
			end
		end
	end
end

function peturbIC(ICvec,seed)
	#peturbIC by 10% of nominal value
	srand(seed)
	selectedindices = [1,3]
	genrand = randn(1, size(ICvec,1))
	for j in selectedindices
		ICvec[j] = ICvec[j]+ICvec[j]*.05*(.5-genrand[j])*2
	end
	return ICvec
end

function setIC(IC, exp_index)
	if(exp_index==2)
		IC[1]= IC[1]*1.55
		#IC[3] =IC[3]*.95
		#IC[5] = IC[5]*.95
		IC[7] = IC[7]*1.5
		IC[8]= IC[8]*1.5
		IC[11] =.01
	elseif(exp_index== 3)
		IC[1]= IC[1]*1.15
		IC[6] = IC[6]*.75
		IC[7] = IC[7]*1.2
		IC[11] = .00
	elseif(exp_index ==4)
		IC[1] = IC[1]*.75
		IC[7] = IC[7]*.95
	elseif(exp_index ==5)
		IC[1] = IC[1]*.65
		IC[7] = IC[7]*.6
	elseif(exp_index ==6)
		IC[1] = IC[1]*.65
		IC[7] = IC[7]*.55
	end
	return IC
end

function setCompleteModelIC(IC, patient_id)
	# 1 FII
   	 # 2 FIIa
   	 # 3 PC
	 # 4 APC
	 # 5 ATIII
    	# 6 TM
	# 7 TRIGGER
	#8 Fraction of platelets activated
 	#9  FV+FX
	#10  FVa+FXa
	#11  prothombinase complex
	#12 Fibrin
	#13 plasmin
  	#14 Fibrinogen
	#15 plasminogen
  	#16 tPA
	#17 uPA
	#18 fibrin monomer
	#19 protofibril
	#20 antiplasmin
    	#21 PAI_1
	#22 Fiber
	println(string("Adjusting IC for patient", patient_id))
	if(patient_id==3)
		IC[1]=IC[1]*1.35
		IC[5] = IC[5]*.9
		IC[7] = IC[7]*1.01
		IC[14] = IC[14]*1.1
		IC[15] = IC[15]*.9
		IC[20] = IC[20]*1.1
	elseif(patient_id==4)
		IC[1]=IC[1]*1.05
		IC[15] = IC[15]*.8
		IC[14]=IC[14]*.5
		IC[16] = IC[16]*.90
		#IC[20]=IC[20]*1.35
		#IC[21] = IC[21]*.75
	elseif(patient_id==5)
		IC[1]=IC[1]*1.4
		IC[5]=IC[5]*.95
		#IC[20]=IC[20]*1.35
#	elseif(patient_id==6)
#		IC[1]=IC[1]*.95
#	elseif(patient_id==7)
#		IC[1] = IC[1]*.95
	elseif(patient_id==8)
		IC[1]= IC[1]*.95
		#IC[20]=IC[20]*1.35
#	elseif(patient_id==9)
#		IC[1]=IC[1]*1.05
#		IC[20]=IC[20]*1.35
#	elseif(patient_id==10)
#		#IC[1]=IC[1]*1.05
	end
	return IC
end

function checkForDynamics(alldata)
	num_param_sets = size(alldata,1)
	threshold = 10 #it has dynamics if it creates 10 thrombin
	hasdynamics = zeros(1, num_param_sets) #if it has dynamics, 1, else zero
	for j in collect(1:num_param_sets)
		ts = alldata[j,:] #get the time series data generated by this parameter set
		mid = ts[Int(floor(end/2))] #get the approximate midpoint
		if(mid>threshold)
			hasdynamics[j] = 1
		end
	end
	return hasdynamics
end

@everywhere function checkForDynamics(thrombin, t)
	threshold = 10 #it has dynamics if it creates 10 thrombin
	hasdynamics = false
	mid = thrombin[Int(floor(end/2))] #get the approximate midpoint
	if(mid>threshold || maximum(thrombin)> threshold)
		hasdynamics = true
	end
	return hasdynamics
end

@everywhere function convertToROTEM(t,x, tPA)
	F = [a[12] for a in x]+ [a[18] for a in x]+ [a[19] for a in x]+ [a[22] for a in x] # fibrin related species 12,18,19,22
	A0 = .01 #baseline ROTEM signal
	K = 5.0-3.75*tPA
	#K = 1
	n = 2
	if(tPA ==2)
		S = 60
		#S = 3.5
	else
		S = 60
		#S = 1.5
	end
	A1 = S
	A = A0+A1.*F.^n./(K.^n+F.^n)
	return A
end

@everywhere function setROTEMIC(tPA, ID)
	@show ID
	pathToData="../data/Viscoelasticmeasurements.xlsx"
	all_platelets = Dict("3"=>189, "4"=>208, "5"=>210, "6"=>263, "7"=>194, "8"=>190, "9"=>149, "10"=>195)
	look_ups_0_tPA = Dict("3"=>"Timecourse!BA3:BC1387", "4"=>"Timecourse!AT3:AV1835", "5"=>"Timecourse!AM3:AO1901", "6"=>"Timecourse!AF3:AH988", "7"=>"Timecourse!Y3:AA1447", "8"=>"Timecourse!R3:T2069", "9"=>"Timecourse!K3:M1544", "10"=>"Timecourse!D3:F789")
	look_ups_2_tPA = Dict("3"=>"Timecourse!AX3:AZ1528", "4"=>"Timecourse!AQ3:AS1682", "5"=>"Timecourse!AJ3:AL1111", "6"=>"Timecourse!AC3:AE1097", "7"=>"Timecourse!V3:X998", "8"=>"Timecourse!O3:Q1154", "9"=>"Timecourse!H3:J1166", "10"=>"Timecourse!A3:C1301")
	currPlatelets = all_platelets[string(ID)]
	datastr = ""
	if(tPA ==0)
		datastr=look_ups_0_tPA[string(ID)]
	else
		datastr=look_ups_2_tPA[string(ID)]
	end
	@show datastr
	data=readxl(pathToData, datastr)
	data=Array{Float64}(data)
	time = data[:,1]
	avg_run = mean(data[:,2:3],2);
	exp_data = hcat(time/60, avg_run) #convert to minute from seconds
	return currPlatelets, exp_data
end

function plotAverageROTEMWData(t,meanROTEM,stdROTEM,expdata, savestr)
	fig = figure(figsize = (15,15))
	ylabel("ROTEM")
	xlabel("Time, in minutes")
	plot(t, transpose(meanROTEM), "k")
	axis([0, t[end], 0, 100])
	@show size(meanROTEM)
	@show size(stdROTEM)
	@show size(t)
	upper = transpose(meanROTEM+stdROTEM)
	lower = transpose(meanROTEM-stdROTEM)
	@show size(vec(upper))
	@show size(vec(lower))
	fill_between((t), vec(upper), vec(lower), color = ".5", alpha =.5)
	plot(expdata[:,1], expdata[:,2], ".k")
	savefig(savestr)
end

function plotAverageROTEMWDataSubplot(fig,t,meanROTEM,stdROTEM,expdata)
	plot(t, transpose(meanROTEM), "k")
	#axis([0, t[end], 0, 150])
	@show size(meanROTEM)
	@show size(stdROTEM)
	@show size(t)
	upper = transpose(meanROTEM+stdROTEM)
	lower = transpose(meanROTEM-stdROTEM)
	@show size(vec(upper))
	@show size(vec(lower))
	fill_between((t), vec(upper), vec(lower), color = ".5", alpha =.5)
	plot(expdata[:,1], expdata[:,2], ".k")
	return fig
end

function makeAllPredictions()
	pathToParams="parameterEstimation/Best11OverallParameters_19_04_2017.txt"
	ids = [3,4,9,10]
	tPAs = [0,2]
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("figures/UsingBest11PredictingPatient", ids[j], "_tPA=", tPAs[k], "_19_04_2017.pdf")
			testROTEMPredicition(pathToParams, ids[j], tPAs[k], savestr)
		end
	end
end

function makeAllEstimatedCurves()
	pathToParams="parameterEstimation/Best11OverallParameters_19_04_2017.txt"
	ids = [5,6,7,8]
	tPAs = [0,2]
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("figures/Patient", ids[j], "_tPA=", tPAs[k], "_18_04_2017.pdf")
			alldata, meanROTEM, stdROTEM=testROTEMPredicition(pathToParams, ids[j], tPAs[k], savestr)
		end
	end
end

function makeTrainingFigure()
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	#pathToParams="parameterEstimation/Best11OverallParameters_19_04_2017.txt"
	POETs_data = "parameterEstimation/POETS_info_24_05_2017maxstep1_OriginalShapeFunction.txt"
	ec,pc,ra=parsePOETsoutput(POETs_data)
	ids = [5,6,7,8]
	tPAs = [0,2]
	close("all")
	fig=figure(figsize = [15,15])
	counter = 1
	numParamSets = 15
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("figures/Patient", ids[j], "_tPA=", tPAs[k], "_18_04_2017.pdf")
			bestparams=generateNbestGivenObjective(numParamSets,ec, pc,counter)
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParams(bestparams, ids[j], tPAs[k], savestr)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			@show counter
			plt[:subplot](4,2,counter)
			fig=plotAverageROTEMWDataSubplot(fig,TSIM,meanROTEM,stdROTEM,expdata)
			if(mod(counter,2)==1)
				axis([0, TSIM[end], 0, 70])
				ylabel(string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 90])
			end

			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	annotate("tPA = 0 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	savefig(string("figures/trainingFigureUsing",numParamSets, "ParameterSets_22_05_17OriginalShapeFunction.pdf"))
end

function makeTrainingFigureBestOveralParams()
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
	ids = [5,6,7,8]
	tPAs = [0,2]
	close("all")
	fig=figure(figsize = [15,15])
	counter = 1
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "_18_04_2017.pdf")
			datasavestr = string("../generatedData/Patient", ids[j],"_tPA=", tPAs[k], "_31_05_2017.txt" )
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicition(pathToParams, ids[j], tPAs[k], savestr)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			writetodisk= hcat(TSIM, transpose(meanROTEM), transpose(stdROTEM))
			writedlm(datasavestr, writetodisk)
			@show counter
			plt[:subplot](size(ids,1),size(tPAs,1),counter)
			fig=plotAverageROTEMWDataSubplot(fig,TSIM,meanROTEM,stdROTEM,expdata)
			if(mod(counter,2)==1)
				axis([0, TSIM[end], 0, 80])
				ylabel(string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 80])
			end
			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	annotate("tPA = 0 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	savefig("../figures/TrainingFigureUsingBest2ParamSetPerObj_25_05_17OriginalShapeFunctionOnlyFittPA2.pdf")
end

function makePredictionsFigure()
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
	ids = [3,4,9,10]
	tPAs = [0,2]
	close("all")
	fig=figure(figsize = [15,15])
	counter = 1
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "_31_05_2017.pdf")
			datasavestr = string("../generatedData/Patient", ids[j],"_tPA=", tPAs[k], "_31_05_2017.txt" )
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicition(pathToParams, ids[j], tPAs[k], savestr)
			@show size(TSIM), meanROTEM, stdROTEM
			writetodisk= hcat(TSIM, transpose(meanROTEM), transpose(stdROTEM))
			writedlm(datasavestr, writetodisk)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			@show counter
			plt[:subplot](size(ids,1),size(tPAs,1),counter)
			fig=plotAverageROTEMWDataSubplot(fig,TSIM,meanROTEM,stdROTEM,expdata)
			if(mod(counter,2)==1)
				axis([0, TSIM[end], 0, 80])
				ylabel(string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 80])
			end
			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	annotate("tPA = 0 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	savefig("../figures/PredictionsFigureUsingBest2ParamSetPerObj_31_05_17OriginalShapeFunctionOnlyFittPA2.pdf")
end

function testROTEMPredicitionGivenParams(allparams,patient_id,tPA,savestr)
	numparams = 77
	pathToThrombinData="../data/fromOrfeo_Thrombin_HT_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 60.0
	end
	TSIM = collect(TSTART:Ts:TSTOP)
	platelets,usefuldata = setROTEMIC(tPA, patient_id)
	platelet_count =platelets
	alldata = zeros(1,size(TSIM,1))
	if(size(allparams,1)==numparams) #deal with parameters being stored either vertically or horizontally
		itridx = 2
	else
		itridx = 1
	end
	
	for j in collect(1:size(allparams,itridx))
		if(itridx ==2)
			currparams = vec(allparams[:,j])
		else
			currparams = vec(allparams[j,:])
		end
		#@show currparams
		if(typeof(currparams)==Array{Array,1}) #deal with params being inside an extra layer of array
			currparams= currparams[1]
		end
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		initial_condition_vector=setCompleteModelIC(initial_condition_vector,patient_id)
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(t,y)= BalanceEquations(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		#@show size([a[2] for a in X])
		A = convertToROTEM(t,X,tPA)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	#plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
	return alldata, meanROTEM, stdROTEM, TSIM
end



function testROTEMPredicition(pathToParams,patient_id,tPA,savestr)
	numparams = 77
	allparams = readdlm(pathToParams, '\t')
	pathToThrombinData="../data/fromOrfeo_Thrombin_HT_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 60.0
	end
	TSIM = collect(TSTART:Ts:TSTOP)
	platelets,usefuldata = setROTEMIC(tPA, patient_id)
	platelet_count =platelets
	alldata = zeros(1,size(TSIM,1))
#	@show size(alldata)
#	@show size(allparams)
#	@show allparams
	if(size(allparams,1)==numparams) #deal with parameters being stored either vertically or horizontally
		itridx = 2
	else
		itridx = 1
	end
	
	for j in collect(1:size(allparams,itridx))
		if(itridx ==2)
			currparams = vec(allparams[:,j])
		else
			currparams = vec(allparams[j,:])
		end
		@show currparams
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		initial_condition_vector=setCompleteModelIC(initial_condition_vector,patient_id)
		#@show dict
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(t,y)= BalanceEquations(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		#@show size([a[2] for a in X])
		A = convertToROTEM(t,X,tPA)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	#plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function testROTEMPredicitionAndPlot(pathToParams,patient_id,tPA,savestr)
	close("all")
	numparams = 77
	allparams = readdlm(pathToParams, '\t')
	pathToThrombinData="../data/fromOrfeo_Thrombin_HT_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 60.0
	end
	TSIM = collect(TSTART:Ts:TSTOP)
	platelets,usefuldata = setROTEMIC(tPA, patient_id)
	fig1 = figure(figsize = (15,15))
	fig2 = figure(figsize = (15,15))
	fig3 = figure(figsize = (15,15))
	platelet_count =platelets
	alldata = zeros(1,size(TSIM,1))
#	@show size(alldata)
#	@show size(allparams)
#	@show allparams
	if(size(allparams,1)==numparams) #deal with parameters being stored either vertically or horizontally
		itridx = 2
	else
		itridx = 1
	end
	
	for j in collect(1:size(allparams,itridx))
		if(itridx ==2)
			currparams = vec(allparams[:,j])
		else
			currparams = vec(allparams[j,:])
		end
		@show currparams
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		@show dict
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(t,y)= BalanceEquations(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		figure(1)
		plotThrombinWData(t,X,pathToThrombinData)
		figure(2)
		makeLoopPlots(t,X)
		#@show alldata
		#@show size([a[2] for a in X])
		A = convertToROTEM(t,X,tPA)
		figure(3)
		plot(t, A)
		@show size(A)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function generateSobolParams()
	#let's use the averaged top 8 parameter sets
	startingpt =  readdlm("parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	outputfn="sensitivity/sobolboundspm50percentOnlyParams_05_30_17.txt"
	best = mean(startingpt,1)
	data_dictionary=buildCompleteDictFromOneVector(best)
	names = data_dictionary["parameter_name_mapping_array"]
	str = ""
	j = 1
	for name in names
		currstr = string(name, " ", best[j]*.5, " ", best[j]*1.5, "\n")
		str = string(str,currstr)
		j=j+1
	end
	touch(outputfn)
	f = open(outputfn, "a")
	write(f,str)
	close(f)
	
end

function generateSobolParamsForPerturbIC()
	#let's use the averaged top 8 parameter sets
	startingpt =  readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	outputfn="../sensitivity/sobolboundspm50percent_05_31_17.txt"
	meanparams = mean(startingpt,1)
	data_dictionary=buildCompleteDictFromOneVector(meanparams)
	names = data_dictionary["parameter_name_mapping_array"]
	initial_conditions = data_dictionary["INITIAL_CONDITION_VECTOR"]
	IC_names = ["FII", "FIIa", "PC", "APC", "ATIII", "TM", "TRIGGER", "Fraction_Platelets_Activated", "FV+FX", "[FV+FX]a", "Prothrombinase_Complex", "Fibrin", "Plasmin", "Fibrinogen", "Plasminogen", "tPA", "uPA", "Fibrin_Monomer", "Protofibril", "Antiplasmin", "PAI1", "Fiber"]
	str = ""
	j = 1
	for name in names
		currstr = string(name, " ", meanparams[j]*.5, " ", meanparams[j]*1.5, "\n")
		str = string(str,currstr)
		j=j+1
	end
	j = 1
	for name in IC_names
		if(initial_conditions[j]==0)
			lb = 0.0
			up = 1.0
		else
			lb = initial_conditions[j]*.5
			up = initial_conditions[j]*1.5
		end
		currstr = string("initial_", name, " ", lb, " ", up, "\n")
		str = string(str, currstr)
		j = j+1
	end

	touch(outputfn)
	f = open(outputfn, "a")
	write(f,str)
	close(f)
end

function generateSobolParamsForOnlyPerturbIC()
	#let's use the averaged top 8 parameter sets
	startingpt =  readdlm("parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	outputfn="sensitivity/sobolboundspm50percentOnlyIC_05_30_17.txt"
	meanparams = mean(startingpt,1)
	data_dictionary=buildCompleteDictFromOneVector(meanparams)
	names = data_dictionary["parameter_name_mapping_array"]
	initial_conditions = data_dictionary["INITIAL_CONDITION_VECTOR"]
	IC_names = ["FII", "FIIa", "PC", "APC", "ATIII", "TM", "TRIGGER", "Fraction_Platelets_Activated", "FV+FX", "[FV+FX]a", "Prothrombinase_Complex", "Fibrin", "Plasmin", "Fibrinogen", "Plasminogen", "tPA", "uPA", "Fibrin_Monomer", "Protofibril", "Antiplasmin", "PAI1", "Fiber"]
	str = ""
	j = 1
	for name in IC_names
		if(initial_conditions[j]==0)
			lb = 0.0
			up = 1.0
		else
			lb = initial_conditions[j]*.5
			up = initial_conditions[j]*1.5
		end
		currstr = string("initial_", name, " ", lb, " ", up, "\n")
		str = string(str, currstr)
		j = j+1
	end

	touch(outputfn)
	f = open(outputfn, "a")
	write(f,str)
	close(f)
end

function concatSobolResults()
	filestr1 = "sensitivity/05_30_17_AUCForSobolPM50PercentOnlyICN2000_"
	filestr2="_of_8.txt"
	str = ""
	for j in collect(1:8)
		fn = string(filestr1, j, filestr2)
		currstr = replace(readstring(fn), ",", " ")
		str = string(str, currstr)
	end
	write("sensitivity/AllSobol_05_30_17_AUCForSobolOnlyICPM50PercentN2000.txt", str)
end

function checkThermoFeasability(params)
	data_dictionary = buildCompleteDictFromOneVector(params)
	kinetic_parameter_vector = data_dictionary["KINETIC_PARAMETER_VECTOR"]
	fibrin_kinetic_parameter_vector = data_dictionary["FIBRIN_KINETIC_PARAMETER_VECTOR"]
	upper_limit = 1E8 #kcat/Km 1/(M s)
	num_rates = 21
	eff_rate = zeros(num_rates,1)
	k_trigger = kinetic_parameter_vector[1]
	K_trigger = kinetic_parameter_vector[2]
	k_amplification = kinetic_parameter_vector[3]
	K_FII_amplification = kinetic_parameter_vector[4]
	k_APC_formation = kinetic_parameter_vector[5]
	K_PC_formation = kinetic_parameter_vector[6]
	k_inhibition = kinetic_parameter_vector[7]
	K_FIIa_inhibition = kinetic_parameter_vector[8]
	k_inhibition_ATIII = kinetic_parameter_vector[9]
	k_FV_X_activation = kinetic_parameter_vector[10]
	K_FV_X_actiation = kinetic_parameter_vector[11]
	#k_FX_activation = kinetic_parameter_vector[12]
	#K_FX_activation = kinetic_parameter_vector[13]
	k_complex = kinetic_parameter_vector[14]
	k_amp_prothombinase = kinetic_parameter_vector[15]
	K_FII_amp_prothombinase=kinetic_parameter_vector[16]
	k_amp_active_factors= kinetic_parameter_vector[17]
	K_amp_active_factors = kinetic_parameter_vector[18]

	k_cat_Fibrinogen = fibrin_kinetic_parameter_vector[1]
	Km_Fibrinogen=fibrin_kinetic_parameter_vector[2]
    	k_fibrin_monomer_association=fibrin_kinetic_parameter_vector[3]
	k_protofibril_association=fibrin_kinetic_parameter_vector[4]
	k_protofibril_monomer_association=fibrin_kinetic_parameter_vector[5]
	k_cat_plasminogen_Fibrin_tPA=fibrin_kinetic_parameter_vector[6]
	Km_plasminogen_Fibrin_tPA=fibrin_kinetic_parameter_vector[7]
	k_cat_plasminogen_Fibrin_uPA=fibrin_kinetic_parameter_vector[8]
	Km_plasminogen_Fibrin_uPA=fibrin_kinetic_parameter_vector[9]
	k_cat_Fibrin=fibrin_kinetic_parameter_vector[10]
	Km_Fibrin=fibrin_kinetic_parameter_vector[11]
	k_inhibit_plasmin=fibrin_kinetic_parameter_vector[12]
	k_inhibit_PAI_1_APC=fibrin_kinetic_parameter_vector[13]
	k_inhibit_tPA_PAI_1=fibrin_kinetic_parameter_vector[14]
	k_inhibit_uPA_PAI_1=fibrin_kinetic_parameter_vector[15]
	k_fibrin_formation=fibrin_kinetic_parameter_vector[16]
	k_cat_fiber=fibrin_kinetic_parameter_vector[17]
	Km_fiber=fibrin_kinetic_parameter_vector[18]
	Ka_plasminogen_Fibrin_tPA=fibrin_kinetic_parameter_vector[19]
	K_plasminogen_Fibrin_tPA = fibrin_kinetic_parameter_vector[20]
	k_cat_fibrinogen_deg=fibrin_kinetic_parameter_vector[21]
	Km_fibrinogen_deg=fibrin_kinetic_parameter_vector[22]

	eff_rate[1] = k_trigger/K_trigger
	eff_rate[2]=k_amplification/K_FII_amplification
	eff_rate[3] = k_APC_formation/K_PC_formation
	eff_rate[4]=#k_inhibition_ATIII
	eff_rate[5]=#k_complex
	eff_rate[6]= k_amp_prothombinase/(K_FII_amp_prothombinase)
	eff_rate[7] = k_amp_active_factors/K_amp_active_factors
	eff_rate[8]=k_cat_Fibrinogen/Km_Fibrinogen
	eff_rate[9]=#k_fibrin_monomer_association
	eff_rate[10]=#k_protofibril_association
	eff_rate[11]=#k_protofibril_monomer_association
	eff_rate[12]=#k_fibrin_formation
	eff_rate[13]=k_cat_plasminogen_Fibrin_tPA/Km_plasminogen_Fibrin_tPA
	eff_rate[14]=k_cat_plasminogen_Fibrin_uPA/Km_plasminogen_Fibrin_uPA
	eff_rate[15]=k_cat_Fibrin/Km_Fibrin
	eff_rate[16]=#k_inhibit_plasmin
	eff_rate[17]=#k_inhibit_PAI_1_APC
	eff_rate[18]=#k_inhibit_tPA_PAI_1
	eff_rate[19]=#k_inhibit_uPA_PAI_1
	eff_rate[20]=k_cat_fiber/Km_fiber
	eff_rate[21]=k_cat_fibrinogen_deg/Km_fibrinogen_deg

	return eff_rate
	
end

