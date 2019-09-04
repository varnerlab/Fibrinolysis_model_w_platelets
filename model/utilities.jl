#using PyPlot
using ExcelReaders
#for helvetical
#using PyCall
#PyCall.PyDict(matplotlib["rcParams"])["font.sans-serif"] = ["Helvetica"]
using Distances #for calculated eucliedian distance
using LinearAlgebra
using DelimitedFiles
using Statistics
using DifferentialEquations

 function calculateMSE(t,predictedCurve, experimentalData)
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
	for j in collect(1:maximum(size(predictedCurve)))
		sum = sum +(predictedCurve[j]-interpolatedExperimentalData[j])^2
	end
#	figure()
#	plot(t, interpolatedExperimentalData, "b")
#	plot(experimentalData[:,1], experimentalData[:,2], "--", linewidth = 2.0)
#	plot(t, predictedCurve', "k")
	#@show sum, maximum(size(predictedCurve))
	return sum/maximum(size(predictedCurve)), interpolatedExperimentalData#MSE
end

 function linearInterp(lowerVal, upperVal, tstart, tend,tdesired)
	
	val = lowerVal + (upperVal-lowerVal)/(tend-tstart)*(tdesired-tstart)
	if(isnan(val))
		val = 0.0
	end
	return val
end

 function calculateAUC(t,y)
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

 function buildDictFromOneVector(vector)
	kinetic_parameter_vector = vector[1:18]
	control_parameter_vector=vector[19:38]
	platelet_parameter_vector=vector[39:44]
	timing = vector[45:46]
	platelet_count=vector[47]
	dict = buildCoagulationModelDictionary(kinetic_parameter_vector, control_parameter_vector, platelet_parameter_vector, timing, platelet_count)
	return dict
end

 function buildCompleteDictFromOneVector(vector)
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


 function createCorrectDict(basic_dict, exp_index)
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

 function generateBestNparameters(n, ec_array, pc_array)
	#calculate error
	best_params = Array[]
	total_error = sum(ec_array[:,1:end],1)
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
	writedlm(string("../parameterEstimation/Best", n, "OverallParameters_29_04_2018.txt"), best_params)
	return best_params

end

 function generateNbestPerObjective(n,ec_array, pc_array)
	num_objectives =size(ec_array,1)
	#best_params=Array{Array}(num_objectives*n)
	best_params=zeros(num_objectives*n,size(pc_array,1))
	#best_params =Array{Array, num_objectives*n}
	counter = 1
	for j in collect(1:num_objectives)
		curr_error=ec_array[j,:]
		allidx = collect(1:size(pc_array,2))
		removed = allidx
		for k in collect(1:n)
			allidx = collect(1:size(pc_array,2))
			removed = allidx
			#min_index = indmin(curr_error)
			min_index = argmin(curr_error)
			curr_best_params = pc_array[:,min_index]
			best_params[counter,:] = curr_best_params
			removed = deleteat!(removed, min_index)
			@show min_index, ec_array[:,min_index]
			#@show curr_best_params
			#delete the best ones we've found
			#@show size(total_error)
			@show size(pc_array)
			pc_array=pc_array[:,removed]
			curr_error=deleteat!(vec(curr_error),min_index)
			@show size(pc_array)
			counter=counter+1
		end
	end	
	writedlm(string("../parameterEstimation/Best", n, "PerObjectiveParameters_01_02_19PlateletContributionToROTEM.txt"), best_params)
	#writedlm(string("../parameterEstimation/Best", n, "PerObjectiveParameters_05_12_18PlateletContributionToROTEM.txt"), best_params)
	return best_params
end

 function generateNbestPerObjective(n,ec_array, pc_array,savestring)
	num_objectives =size(ec_array,1)
	best_params=Array{Array}(num_objectives*n)
	counter = 1
	for j in collect(1:num_objectives)
		curr_error=ec_array[j,:]
		allidx = collect(1:size(pc_array,2))
		removed = allidx
		for k in collect(1:n)
			allidx = collect(1:size(pc_array,2))
			removed = allidx
			min_index = indmin(curr_error)
			curr_best_params = pc_array[:,min_index]
			best_params[counter] = curr_best_params
			removed = deleteat!(removed, min_index)
			@show min_index, ec_array[:,min_index]
			#@show curr_best_params
			#delete the best ones we've found
			#@show size(total_error)
			pc_array=pc_array[:,removed]
			curr_error=deleteat!(vec(curr_error),min_index)
			@show size(pc_array)
			counter=counter+1
		end
	end	
	writedlm(savestring, best_params)
	return best_params
end

 function generateNbestGivenObjective(n,ec_array, pc_array,objectivenum)
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
	#alltext = readstring(f)
	alltext = read(f, String)
	close(f)

	outputname = "textparsing.txt"
	number_of_parameters = 77
  	number_of_objectives = 4
	ec_array = zeros(number_of_objectives)
  	pc_array = zeros(number_of_parameters)
	rank_array = zeros(1)	
	counter =1
	#for grouping in matchall(r"\[([^]]+)\]", alltext)
	for grouping in collect(m.match for m in eachmatch(r"\[([^]]+)\]", alltext))
		cleanedgrouping = replace(grouping, "["=>"")
		nocommas = replace(cleanedgrouping, ","=>" ")
		allcleaned = replace(nocommas, "]"=>"")
		allcleaned = replace(allcleaned, ";"=>"\n")
		outfile = open(outputname, "w")
		write(outfile, allcleaned)
		close(outfile)
		formatted = readdlm(outputname)
		#@show formatted	
		#@show size(formatted), counter
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

function parsePOETsoutput(filename, number_of_objectives)
	f = open(filename)
	#alltext = readstring(f)
	alltext = read(f, String)
	close(f)

	outputname = "textparsing.txt"
	number_of_parameters = 77
	ec_array = zeros(number_of_objectives)
  	pc_array = zeros(number_of_parameters)
	rank_array = zeros(1)	
	counter =1
	for grouping in collect(m.match for m in eachmatch(r"\[([^]]+)\]", alltext))
		cleanedgrouping = replace(grouping, "["=>"")
		nocommas = replace(cleanedgrouping, ","=>" ")
		allcleaned = replace(nocommas, "]"=>"")
		allcleaned = replace(allcleaned, ";"=>"\n")
		outfile = open(outputname, "w")
		write(outfile, allcleaned)
		close(outfile)
		formatted = readdlm(outputname)
		#@show formatted	
		#@show size(formatted), counter
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
	for j in idx1
		for k in idx2
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
		IC[1]=IC[1]*1.45
		IC[3]=IC[3]*.7
		IC[5] = IC[5]*.6
		IC[6]=IC[6]*.6
		IC[7] = IC[7]*1.05
##		IC[9] = IC[9]*1.15
		IC[14] = IC[14]*1.3
#		IC[15] = IC[15]*.93
	#	IC[20] = IC[20]*1.2
		IC[21]=IC[21]*.5
	elseif(patient_id==4)
		IC[1]=IC[1]*.9
		IC[3]=IC[3]*1.2
		IC[5] = IC[5]*1.2
		IC[6]=IC[6]*1.2
		#IC[7] = IC[7]*1.05
		#IC[14] = IC[14]*1.2
#		IC[15] = IC[15]*1.15
	elseif(patient_id==5)
		IC[1]=IC[1]*.8
		IC[5]=IC[5]*1.2
		IC[6] = IC[6]*1.1
#		IC[7] = IC[7]*1.05
#		IC[14]=IC[14]*1.2
		IC[20]=IC[20]*.85
	elseif(patient_id==6)
		#IC[1]=IC[1]*.65				
		#IC[3]=IC[3]*1.2
	#	IC[5] = IC[5]*1.2
#		IC[6]=IC[6]*1.2
#		IC[9] = IC[9]*.8
	#	IC[14] = IC[14]*.95
#		IC[15] = IC[15]*1.25
	elseif(patient_id==7)
#		IC[1] = IC[1]*.8
#		IC[3] = IC[3]*1.2
#		IC[5] = IC[5]*1.2
#		IC[14] = IC[14]*1.1
#		IC[20]=IC[20]*.75
	elseif(patient_id==8)
#		IC[1]= IC[1]*.8
#		IC[3] = IC[3]*1.2
#		IC[5] = IC[5]*1.2
#		IC[15] = IC[15]*1.2
#		IC[20]=IC[20]*.75
#		IC[3] = IC[3]*1.2
#		IC[5] = IC[5]*1.2
		IC[14] = IC[14]*.95
#		IC[20]=IC[20]*.75
	elseif(patient_id==9)
#		IC[1]=IC[1]*1.05
#		IC[5] = IC[5]*.8
#		IC[7] = IC[7]*1.03
#		IC[20]=IC[20]*1.35
	elseif(patient_id==10)
		IC[15] = IC[15]*1.2
		IC[1]=IC[1]*1.05
		IC[13]=IC[13]*1.2
		IC[20]=IC[20]*.8
	end
	return IC
end

function setICsBasedOnNM(initial_condition_vector, patient_id)
	basestr = "../parameterEstimation/ICEstimation/MinimumICsForPatient"
	if(patient_id==6 || patient_id==7 ||patient_id==8)
		#set manually
		params =setCompleteModelIC(initial_condition_vector, patient_id)
		return params
	else
		filestr = string(basestr, patient_id, "usingParametersFrom_05_12_18UsingNM_W15Evals.txt")
	end
	@show filestr
	params = [] #appearently, need to define params before try catch statement
	try 
		params = readdlm(filestr, '\n')
	catch 
		println("No file found-will use unadjusted ICS")
		params = initial_condition_vector
	end
	#@show size(params)
	return params
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

 function checkForDynamics(thrombin, t)
	threshold = 10 #it has dynamics if it creates 10 thrombin
	hasdynamics = false
	mid = thrombin[Int(floor(end/2))] #get the approximate midpoint
	if(mid>threshold || maximum(thrombin)> threshold)
		hasdynamics = true
	end
	return hasdynamics
end

 function convertToROTEM(t,x, tPA)
	#to work with Differential Equations
	if(contains(string(typeof(x)), "DiffEqBase"))
		F = x[12,:]+x[18,:]+ x[19,:]+x[22,:]
	else
		F = [a[12] for a in x]+ [a[18] for a in x]+ [a[19] for a in x]+ [a[22] for a in x] # fibrin related species 12,18,19,22
	end
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

 function convertToROTEMPlateletContribution(t,x, tPA,platelet_count)
	#to work with Differential Equations
	#@show typeof(x)
	if(occursin("DiffEq",string(typeof(x))))
		F = x[12,:]+x[18,:]+ x[19,:]+x[22,:]
	else
		F = [a[12] for a in x]+ [a[18] for a in x]+ [a[19] for a in x]+ [a[22] for a in x] # fibrin related species 12,18,19,22
	end
	normal_platelet_count = 300 #*10^6 #/mL
	A0 = fill(.01, size(F)) #baseline ROTEM signal
	#K = 2000-375*tPA
	#K = 2000-200*tPA
	#K = 5000-375*tPA#-use me on Monday to check predictions
	#K = 5000-1000*tPA
	#K = 2000+1125*tPA
	K = fill(1000+100*tPA, size(F))
	#K = 1
	n = 2
	#for infomration about weights
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4568902/ Assessing the Methodology for Calculating Platelet Contribution to Clot Strength (Platelet Component) in Thromboelastometry and Thrombelastography
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4238905/ FIBRINOGEN AND PLATELET CONTRIBUTIONS TO CLOT FORMATION: IMPLICATIONS FOR TRAUMA RESUSCITATION AND THROMBOPROPHYLAXIS

	wp = .5 #platelet contribution
	wf = .5 #fibrin related species contribution
	if(tPA ==2)
		S=65
		#S = 60
		#S = 3.5
	else
		S = 65
		#S = 60
		#S = 1.5
	end
	A1 = fill(S, size(F))
	if(occursin("DiffEq", string(typeof(x))))
		#print("here!")
		P=platelet_count/normal_platelet_count.*vec(x[8,:]).*vec(x[12,:])
	else
		P=platelet_count/normal_platelet_count.*[a[8] for a in x].*[a[12] for a in x]
	end
#	figure()
#	plot(t,wf.*F, "r")
#	plot(t, wp.*P, "g")
#	legend(["Fibrin contribution", "Polymerized Platelet Contribution"])
	#x[8]-fraction activated platelets
	R = wf.*F+wp.*P
	A = A0+A1.*R.^n./(K.^n+R.^n)
	return A
end

 function convertToROTEMPlateletContributionScaledF(t,x, tPA,platelet_count,F0)
	#to work with Differential Equations
	if(contains(string(typeof(x)), "DiffEqBase"))
		F = x[12,:]+x[18,:]+ x[19,:]+x[22,:]
	else
		F = [a[12] for a in x]+ [a[18] for a in x]+ [a[19] for a in x]+ [a[22] for a in x] # fibrin related species 12,18,19,22
	end
	normal_platelet_count = 300 #*10^6 #/mL
	A0 = .01 #baseline ROTEM signal
	K = 5000.0-375*tPA
	#K = 1
	n = 2
	#for infomration about weights
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4568902/ Assessing the Methodology for Calculating Platelet Contribution to Clot Strength (Platelet Component) in Thromboelastometry and Thrombelastography
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4238905/ FIBRINOGEN AND PLATELET CONTRIBUTIONS TO CLOT FORMATION: IMPLICATIONS FOR TRAUMA RESUSCITATION AND THROMBOPROPHYLAXIS

	wp = 1/30 #platelet contribution
	wf =0.0 #fibrin related species contribution
	if(tPA ==2)
		S = 60
		#S = 3.5
	else
		S = 60
		#S = 1.5
	end
	A1 = S
	#x[8]-fraction activated platelets
	@show minimum(F)
	@show maximum(F)
	P = platelet_count/normal_platelet_count.*[a[8] for a in x].*F
	@show minimum(P)
	@show maximum(P)
	figure()
	plot(t,F*wf)
	plot(t,wp*P)

	#R = wf.*F+wp.*P
	#A = A0+A1.*R.^n./(K.^n+R.^n)
	A = P*wp
	return A
end

 function setROTEMIC(tPA, ID)
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
	avg_run = mean(data[:,2:3],dims=2);
	exp_data = hcat(time/60, avg_run) #convert to minute from seconds
	return currPlatelets, exp_data
end

function plotAverageROTEMWData(t,meanROTEM,stdROTEM,expdata, savestr)
	fig = figure(figsize = (15,15))
	ylabel("ROTEM")
	xlabel("Time, in minutes")
	plot(t, vec(meanROTEM), "k")
	axis([0, t[end], 0, 100])
	@show size(meanROTEM)
	@show size(stdROTEM)
	@show size(t)
	upper = vec(meanROTEM+stdROTEM)
	lower = vec(meanROTEM-stdROTEM)
	@show size(vec(upper))
	@show size(vec(lower))
	fill_between((t), vec(upper), vec(lower), color = ".5", alpha =.5)
	plot(expdata[:,1], expdata[:,2], ".k")
	savefig(savestr)
end

function plotAverageROTEMWDataSubplot(ax,t,meanROTEM,stdROTEM,expdata)
	ax[:plot](t, transpose(meanROTEM), "k")
	#axis([0, t[end], 0, 150])
	@show size(meanROTEM)
	@show size(stdROTEM)
	@show size(t)
	upper = transpose(meanROTEM+1.96*stdROTEM)
	lower = transpose(meanROTEM-1.96*stdROTEM)
	@show size(vec(upper))
	@show size(vec(lower))
	ax[:fill_between]((t), vec(upper), vec(lower), color = "lightskyblue", alpha =.5)
	ax[:plot](expdata[:,1], expdata[:,2], ".k")
	return ax
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

function makeTrainingFigurePlatletContributionToROTEM()
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	#POETs_data = "../parameterEstimation/POETS_info_28_03_18_PlateletContributionToROTEM.txt"
	#POETs_data="../parameterEstimation/POETS_info_05_04_18_PlateletContributionToROTEM.txt"
	#POETs_data ="../parameterEstimation/POETS_info_11_04_18_PlateletContributionToROTEM.txt"
	#POETs_data = "../parameterEstimation/POETS_info_25_04_18_PlateletContributionToROTEM.txt"
	#POETs_data = "../parameterEstimation/POETS_info_27_04_18_PlateletContributionToROTEM.txt"
	#POETs_data = "../parameterEstimation/POETS_info_29_04_18_PlateletContributionToROTEM.txt"
	#POETs_data = "../parameterEstimation/POETS_info_03_05_18_PlateletContributionToROTEM.txt"
	#POETs_data = "../parameterEstimation/POETS_info_02_05_18_PlateletContributionToROTEMFlatness1.txt" #decent set, reuse 
	#POETs_data ="../parameterEstimation/POETS_info_12_10_18_PlateletContributionToROTEMFlatness1SmallerConversion.txt"
	#POETs_data = "../parameterEstimation/POETS_info_15_10_18_PlateletContributionToROTEMFlatness1SmallerConversion.txt"
	#POETs_data = "../parameterEstimation/POETS_info_19_10_18_PlateletContributionToROTEMFlatness1SmallerConversion.txt"''
	#POETs_data = "../parameterEstimation/POETS_info_05_12_18_PlateletContributionToROTEMFlatness1SmallerConversion.txt"
	POETs_data ="../parameterEstimation/POETS_info_02_01_19_PlateletContributionToROTEMFlatness1SmallerConversion.txt"

	ec,pc,ra=parsePOETsoutput(POETs_data)
	ids = [5,6,7,8]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	numParamSets = 4
	adjustICs = false
	for j in collect(1:size(ids,1))
		@show ids[j]
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "_05_12_2018.pdf")
			bestparams=generateNbestPerObjective(numParamSets,ec,pc)
			#bestparams=generateBestNparameters(8,ec,pc)
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParamsPlatetContributionToROTEM(bestparams, ids[j], tPAs[k], savestr,adjustICs)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			@show norm(meanROTEM), norm(expdata[:,2])
			@show counter
			curraxis=axarr[j,k]
			#plotAverageROTEMWDataSubplot(curraxis,TSIM,meanROTEM,stdROTEM,expdata)
			curraxis[:plot](TSIM, vec(meanROTEM), "k")
			upper = vec(meanROTEM+stdROTEM)
			lower = vec(meanROTEM-stdROTEM)
			curraxis[:fill_between]((TSIM), vec(upper), vec(lower), color = ".5", alpha =.5)
			curraxis[:plot](expdata[:,1], expdata[:,2], ".k")
			if(mod(counter,2)==1)
				curraxis[:set_ylim](0, 70)
				curraxis[:set_xlim](0,TSIM[end])
				curraxis[:set_ylabel](string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 90])
			end

			if(counter==7 || counter ==8)
				#xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	figure(1)
	annotate("tPA = 0 nanomolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 nanomolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	fig[:text](0.5, 0.04, "Time (minutes)", ha="center", va="center", fontsize=40)
	fig[:text](0.06, 0.5, "Ampltiude (mm)", ha="center", va="center", rotation="vertical",fontsize=40)

	savefig(string("../figures/trainingFigureUsing",numParamSets, "ParameterSets_02_01_19PlatletContributionToROTEMICAdjustment=", string(adjustICs),"Best2PerObjNewConversionUsingDE.pdf"))
end

function makePredictionFigurePlatletContributionToROTEM()
	#use me
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	#POETs_data = "../parameterEstimation/POETS_info_28_03_18_PlateletContributionToROTEM.txt"
	#POETs_data="../parameterEstimation/POETS_info_05_04_18_PlateletContributionToROTEM.txt"
	#POETs_data ="../parameterEstimation/POETS_info_19_04_18_PlateletContributionToROTEM.txt"
	#POETs_data="../parameterEstimation/POETS_info_09_04_18_PlateletContributionToROTEMDiffROTEM.txt"
	#POETs_data= "../parameterEstimation/POETS_info_02_05_18_PlateletContributionToROTEMFlatness1.txt"
	#POETs_data = "../parameterEstimation/POETS_info_05_12_18_PlateletContributionToROTEMFlatness1SmallerConversion.txt"
	POETs_data ="../parameterEstimation/POETS_info_02_01_19_PlateletContributionToROTEMFlatness1SmallerConversion.txt"
	ec,pc,ra=parsePOETsoutput(POETs_data)
	ids = [3,4,9,10]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	numParamSets = 2
	adjustICs=false
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "WithoutICAdjustment_05_12_18.pdf")
			bestparams=generateNbestPerObjective(numParamSets,ec,pc)
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParamsPlatetContributionToROTEM(bestparams, ids[j], tPAs[k], savestr,adjustICs)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			@show counter
			@show norm(meanROTEM), norm(expdata[:,2])
			curraxis=axarr[j,k]
			#plotAverageROTEMWDataSubplot(curraxis,TSIM,meanROTEM,stdROTEM,expdata)
			curraxis[:plot](TSIM, vec(meanROTEM), "k")
			upper = vec(meanROTEM+stdROTEM)
			lower = vec(meanROTEM-stdROTEM)
			curraxis[:fill_between]((TSIM), vec(upper), vec(lower), color = ".5", alpha =.5)
			curraxis[:plot](expdata[:,1], expdata[:,2], ".k")
			if(mod(counter,2)==1)
				curraxis[:set_ylim](0, 70)
				curraxis[:set_xlim](0,TSIM[end])
				curraxis[:set_ylabel](string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 90])
			end

			if(counter==7 || counter ==8)
				figure(1)
				#xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	figure(1)
	annotate("tPA = 0 nanomolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 nanomolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	fig[:text](0.5, 0.04, "Time (minutes)", ha="center", va="center", fontsize=40)
	fig[:text](0.06, 0.5, "Ampltiude (mm)", ha="center", va="center", rotation="vertical",fontsize=40)

	savefig(string("../figures/PredictionsFigureUsing",numParamSets, "ParameterSets_01_02_19PlatletContributionToROTEMWithICAdjustment", adjustICs,"ChangedConversionUsingDE.pdf"))
end


function testROTEMPredicitionGivenParams(allparams,patient_id,tPA,savestr)
	numparams = 77
	pathToThrombinData="../data/fromOrfeo_Thrombin_BL_PRP.txt"
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
		fbalances(t,y)= Balances(t,y,dict)
		#tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		#toc()	
		#@show size([a[2] for a in X])
		#A = convertToROTEM(t,X,tPA)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function testROTEMPredicitionGivenParamsPlatetContributionToROTEM(allparams,patient_id,tPA,savestr, adjustICs::Bool)
	numparams = 77
	pathToThrombinData="../data/fromOrfeo_Thrombin_BL_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 100.0
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
		if(adjustICs==true)
			initial_condition_vector = setICsBasedOnNM(initial_condition_vector, patient_id)
		end
		initial_condition_vector[16]=tPA #set tPA level
		#initial_condition_vector=setCompleteModelIC(initial_condition_vector,patient_id)
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(y,p,t)= Balances(t,y,dict) 
		#fbalances(t,y)= Balances(t,y,dict) 
		#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
		prob = DifferentialEquations.ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
		@time sol = DifferentialEquations.solve(prob, saveat=.02)
		t = sol.t
		X=sol	
		#@show size([a[2] for a in X])
		A = convertToROTEMPlateletContribution(t,X,tPA,platelet_count)
		@show size(A), size(TSIM)
		while(size(A)!=size(TSIM))
			fbalances(y,p,t)= Balances(t,y,dict) 
			#fbalances(t,y)= Balances(t,y,dict) 
			#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP), alg_hints=[:stiff])
			@time sol = solve(prob,saveat=.02)
			t = sol.t
			X=sol	
			#@show size([a[2] for a in X])
			A = convertToROTEMPlateletContribution(t,X,tPA,platelet_count)
		end

		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,dims=1)
	stdROTEM = std(alldata,dims=1)
	plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function testROTEMPredicitionGivenParamsPolymerizedPlatelets(allparams,patient_id,tPA,savestr)
	numparams = 77
	pathToThrombinData="../data/fromOrfeo_Thrombin_BL_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 100.0
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
		fbalances(t,y)= Balances(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		#@show size([a[2] for a in X])
		A =  convertToROTEMPlateletContributionScaledF(t,X,tPA,platelet_count, initial_condition_vector[14])
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function testROTEMPredicition(pathToParams,patient_id,tPA,savestr)
	numparams = 77
	allparams = readdlm(pathToParams, '\t')
	pathToThrombinData="../data/fromOrfeo_Thrombin_BL_PRP.txt"
	TSTART = 0.0
	Ts = .02
	if(tPA==0)
		TSTOP =180.0
	else
		TSTOP = 100.0
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
		#@show currparams
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		initial_condition_vector=setCompleteModelIC(initial_condition_vector,patient_id)
		#@show dict
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		#fbalances(t,y)= BalanceEquations(t,y,dict)
		fbalances(t,y)= Balances(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		#@show size([a[2] for a in X])
		A = convertToROTEM(t,X,tPA)
		#plot(t,A)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
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
		#@show currparams
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		#@show dict
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
		#@show size(A)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	plotAverageROTEMWData(TSIM, meanROTEM, stdROTEM, usefuldata,savestr)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function generateAvgROTEMCurve(patient_id,tPA, IC_to_alter)
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
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
		currparams[47]=platelet_count
		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		initial_condition_vector =alterIC(initial_condition_vector,IC_to_alter)
		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(t,y)= BalanceEquations(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		A = convertToROTEM(t,X,tPA)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function generateAvgROTEMCurvePeturbPairwise(patient_id,tPA, IC_to_alter, params_to_alter)
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
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
		currparams[47]=platelet_count
		#we're going to perturb 5x params-params to alter should be 2 long
		currparams[params_to_alter[1]] = currparams[params_to_alter[1]]*5
		currparams[params_to_alter[2]] = currparams[params_to_alter[2]]*5	

		dict = buildCompleteDictFromOneVector(currparams)
		initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
		initial_condition_vector[16]=tPA #set tPA level
		initial_condition_vector =alterIC(initial_condition_vector,IC_to_alter)

		reshaped_IC = vec(reshape(initial_condition_vector,22,1))
		fbalances(t,y)= BalanceEquations(t,y,dict)
		tic() 
		t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0, points=:specified)
		toc()	
		A = convertToROTEM(t,X,tPA)
		alldata=vcat(alldata,transpose(A))
	end
	alldata = alldata[2:end, :] #remove row of zeros
	alldata = map(Float64,alldata)
	meanROTEM = mean(alldata,1)
	stdROTEM = std(alldata,1)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function characterizeCurve(patient_id,tPA)
	tic()
	IC_to_alter = Dict()
	alldata, meanROTEM, stdROTEM, TSIM=generateAvgROTEMCurve(patient_id,tPA,IC_to_alter)
	toc()
	calculateCommonMetrics(meanROTEM,TSIM)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function doPairwisePeturb()
	IC_to_alter = Dict()
	#generate nominal case
	startingpt =  readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	best = mean(startingpt,1)
	@show size(best)
	t,ROTEM_n = runModelWithParamsSetICReturnROTEM(best)
	numparams = 77
	alldistances = zeros(numparams, numparams)
	peturbAmount = 5.0
	for j in collect(1:1:numparams)
		currparams = best
		currparams[j] = best[j]*peturbAmount
		for k in collect(1:1:numparams)
			currparams[k] = best[k]*peturbAmount
			println("On pair ", j, " and ", k)
			#alldata, meanROTEM, stdROTEM, TSIM=generateAvgROTEMCurvePeturbPairwise(patient_id,tPA,IC_to_alter, [j,k])
			t,ROTEM = runModelWithParamsSetICReturnROTEM(currparams)
			@show size(ROTEM)
			@show size(ROTEM_n)
			if(size(ROTEM,1)==size(ROTEM_n,1))	
				currdist = evaluate(Euclidean(), ROTEM_n, ROTEM)
			else
				#if we've had an error in solving ODEs, probably'
				currdist = -1
			end
			@show currdist
			alldistances[j,k]=currdist
			#unpeturb this param
			currparams[k] = best[k]*1/peturbAmount
		end
	end
	writedlm("../generatedData/peturbParams_12_12_17_5x.txt", alldistances)
	return alldistances
end

function characterizeCurve(patient_id,tPA,IC_to_alter)
	tic()
	alldata, meanROTEM, stdROTEM, TSIM=generateAvgROTEMCurve(patient_id,tPA,IC_to_alter)
	toc()
	calculateCommonMetrics(meanROTEM,TSIM)
	return alldata, meanROTEM, stdROTEM, TSIM
end

function compareFibrinogenSup(patient_id,tPA)
	close("all")
#	depleated=Dict("Fibrinogen"=>.5)
#	normal = Dict("Fibrinogen"=>1.0)
#	supplemented=Dict("Fibrinogen"=>1.5)
	f1=figure()		
	allCTs = []

	#cases = [depleated, normal,supplemented]
	cases = collect(.1:.25:2.0)
	for case in cases
		currdict = Dict("Fibrinogen"=>case)
		alldata, meanROTEM, stdROTEM, TSIM=characterizeCurve(patient_id, tPA, currdict)
		CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(meanROTEM,TSIM)
		plot(TSIM, vec(meanROTEM))
		push!(allCTs,CT)
	end
	f2=figure()
	plot(cases,allCTs,"kx")
	xlabel("Percent of Normal Fibrinogen Concentration")
	ylabel("CT (in seconds)")
	savefig("../figures/CT_by_Fibrinogen.pdf")

end

function comparePCCSup(patient_id,tPA)
	close("all")
	f1=figure()		
	allCTs = []
	cases = collect(.1:.1:1.3)
	for case in cases
		currdict = Dict("FII"=>case, "FV_FX"=>case)
		alldata, meanROTEM, stdROTEM, TSIM=characterizeCurve(patient_id, tPA, currdict)
		CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(meanROTEM,TSIM)
		plot(TSIM, vec(meanROTEM))
		push!(allCTs,CT)
	end
	f2=figure()
	plot(cases,allCTs,"kx")
	xlabel("Percent of Normal FII and FV_FX concentration")
	ylabel("CT (in seconds)")
	savefig("../figures/CT_by_PCC.pdf")
end

function comparePCCandFCCSup(patient_id,tPA)
	close("all")
	f1=figure()		
	casesFS = collect(.1:.1:2.0)
	casesPC = collect(.1:.1:2.0)
	allCTs = zeros(size(casesFS,1), size(casesPC,1))
	j =1
	k =1
	for casePC in casesPC
		j = 1
		for caseFS in casesFS
			currdict = Dict("FII"=>casePC, "FV_FX"=>caseFS, "Fibrinogen"=>caseFS)
			alldata, meanROTEM, stdROTEM, TSIM=characterizeCurve(patient_id, tPA, currdict)
			CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(meanROTEM,TSIM)
			plot(TSIM, vec(meanROTEM))
			allCTs[j,k]=CT
			j=j+1
		end
		k = k+1
	end
	f2=figure(figsize=(10,10))
	xgrid = repmat(casesFS,1,maximum(size(casesFS)))
	ygrid = repmat(casesPC',maximum(size(casesPC)),1)
	plot_surface(xgrid,ygrid,allCTs, alpha=.8, cmap=ColorMap("coolwarm"))
	ylabel("Percent of Normal FII and FV_FX concentration")
	xlabel("Percent of Normal Fibrinogen Concentration")
	zlabel("CT in seconds")
	savefig("../figures/CT_by_PCC_and_FS.pdf")
	return casesFS, casesPC, allCTs
end

function createSpeciesNameIdxDict()
	names = ["FII", "FIIa", "PC", "APC", "ATIII", "TM", "TRIGGER", "Fraction Activated Platelets", "FV_FX", "FV_FXa", "Prothombinase-Platelets",
	"Fibrin", "Plasmin", "Fibrinogen", "Plasminogen", "tPA", "uPA", "fibrin monomer", "protofibril", "antiplasmin", "PAI_1", "Fiber"]
	nameIdxDict =Dict(names[i]=> i for i = 1:maximum(size(names)))
	return nameIdxDict
	
end

function alterIC(IC, dict_to_alter)
	#will multiply selected IC in dict_to_alter by the amount specified in the dict
	nameIdxDict = createSpeciesNameIdxDict()
	#@show nameIdxDict
	for key in keys(dict_to_alter)
		IC[nameIdxDict[key]] = IC[nameIdxDict[key]]*dict_to_alter[key] 
	end
	return IC
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

 function calculateCT(ROTEM_curve,TSIM)
	#CT-when we first see an increase in diameter to 2mm
	#for TEG data, R value
	j = 1
	cutoff = 2
	CT = -1
	while(j<maximum(size(ROTEM_curve)))
		if(ROTEM_curve[j]>cutoff)
			CT = TSIM[j]
			break
		end
		j = j+1
	end
	#need to convert to seconds
	return CT*60
end

 function calculateCFT(ROTEM_curve,TSIM)
	#difference between CT and how long it takes to get to 20mm
	#for TEG data, no equivalent
	CT = calculateCT(ROTEM_curve,TSIM)
	cutoff = 20.0
	j =1
	T20 = -1
	if(maximum(ROTEM_curve)<20) #if we never got to 20 mm
		CFT = -1
		return CFT
	end
	while(j<maximum(size(ROTEM_curve)))
		if(ROTEM_curve[j]>cutoff)
			T20 = TSIM[j] #this is in minutes
			break
		end
		j = j+1
	end
	CFT = T20*60-CT
	#need to convert to seconds
	return CFT
end

 function calculateAlpha(ROTEM_curve,TSIM)
	#for TEG, equivalent to calculating angle/alpha
	CFT = calculateCFT(ROTEM_curve,TSIM)
	if(CFT ==-1) #if we've not reached 20 mm
		alpha = -1 #give a nonsense result
	end

	#since we know CT is at 2mm and CFT is at 20mm, we can calculate the slope-this slope is in minutes
	slope = (20-2)/(CFT)*360/(2*pi) #need to convert to degrees
	alpha = atand(slope)#draw a picture-this works
	return alpha
end

 function calculateFirmnessAtTime(ROTEM_curve,TSIM,tdesired)
	j =1
	firmness = -1
	while(j<maximum(size(ROTEM_curve)))
		if(TSIM[j]>=tdesired)
			firmness = ROTEM_curve[j]
			break
		end
		j = j+1
	end
	return firmness
end

 function calculateMCF(ROTEM_curve,TSIM)
	#for TEG, MA
	j = 1
	MCF = -1
	#@show TSIM
	MCF = maximum(ROTEM_curve)
	#@show MCF
	return MCF
end

 function calculateLysisAtTime(ROTEM_curve,TSIM,tdesired)
	MCF = calculateMCF(ROTEM_curve,TSIM)
	later_firmess = calculateFirmnessAtTime(ROTEM_curve,TSIM,tdesired)
	#@show MCF, later_firmess
	LI = 1-(later_firmess)/MCF #since this is expressed as a percent
	if(LI <0)
		LI = 0.0
	end
	return LI
end

 function calculateCommonMetrics(ROTEM_curve,TSIM)
	CT = calculateCT(ROTEM_curve,TSIM)
	CFT = calculateCFT(ROTEM_curve,TSIM)
	alpha = calculateAlpha(ROTEM_curve, TSIM)
	MCF = calculateMCF(ROTEM_curve,TSIM)
	A10 = calculateFirmnessAtTime(ROTEM_curve,TSIM,10)
	A20 = calculateFirmnessAtTime(ROTEM_curve,TSIM,20)
	LI30 = calculateLysisAtTime(ROTEM_curve,TSIM,30)
	LI60 = calculateLysisAtTime(ROTEM_curve,TSIM,60)
#	println("CT: ", CT, " seconds")
#	println("CFT: ", CFT, " seconds")
#	println("alpha: ",alpha, " degrees")
#	println("MCF: ", MCF, " mm")
#	println("A10: ", A10," mm")
#	println("A20: ",A20," mm")
#	println("LI30: ",LI30, " percent of the MCF remains")
#	println("LI60: ", LI60," percent of the MCF remains")
	return CT,CFT,alpha,MCF,A10,A20,LI30,LI60
end

 function checkForFlatness(t,A)
	tstartcheck =5.0
	tendcheck = 40.0
	tstart=findfirst(x -> x == tstartcheck,t)
	tend =findfirst(x -> x == tendcheck,t)
	#go through in 5 minute intervals, and check for flatness
	times = tstartcheck:5:tendcheck
	threshold = .1
	flats = [] #if an interval is flat, gets true. Else, false
	if (maximum(size(t))<=1 || t[end]<tendcheck) #if t is of size 1, we didn't solve properly, so this param set is bad
		AnyFlats = true
		return AnyFlats
	else	
		for j =2:size(times,1)
			tstart = times[j-1]
			tend = times[j]

			tstartidx = findfirst(x -> x >= tstart,t)
			tendidx = findfirst(x -> x >= tend,t)
			#@show tstartidx,tendidx, tstart, tend
			#deal with failure
#			if(maximum(size(tstartidx))==0 || maximum(size(tendidx))==0)
#				AnyFlats = true
#				return AnyFlats
#			end
			tstartidx = tstartidx[1]
			tendidx=tendidx[1]
			#@show tstartidx,tendidx, tstart, tend
			if(tstartidx ==0 && tendidx==0)
				AnyFlats = true
				return AnyFlats
			end

			Aslice = A[tstartidx:tendidx]
			#@show size(Aslice)
			#check for change
			Amin = minimum(Aslice)
			Amax = maximum(Aslice)
			#@show Amin, Amax, tstart, tend

			if(abs(Amax-Amin)>threshold)
				push!(flats, false)

			else
				push!(flats, true)
			end
		end
	#@show flats
		AnyFlats=any(x->x==true, flats)
	end
	return AnyFlats
end

