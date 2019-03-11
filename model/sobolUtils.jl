
function generateSobolParams()
	#let's use the averaged top 8 parameter sets
	#startingpt =  readdlm("parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	#startingpt = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	outputfn="../sensitivity/sobolboundspm50percentOnlyParams_3_08_19.txt"
	allids = collect(11:14)
	outputstr = "../LOOCV/POETS_info_14_02_19_PlateletContributionToROTEMFlatness1ToBeTestedOn" #only tPA 2
	batch_number=10
	numParams = 77
	numObjs = 3 #we had 3 different objectives
	numPatients = 4
	numPerObj = 2 #let's pick the best 2 parameter sets per objective
	currparams = zeros(numPerObj*numObjs*size(allids,1), numParams)
	count = 1
	offset = numPerObj*numObjs
	for id in allids				
		currstr = string(outputstr, id,"Round_1Max_iters10.txtRound_", batch_number, "outOf10.txt")
		#get results from LOOCV
		ec,pc,ra = parsePOETsoutput(currstr, numObjs)
		bestp, besterrs = generateNbestPerObjectiveAndErrors(numPerObj, ec, pc, "temp.txt")
		currparams[(count-1)*offset+1:count*offset,:]=transpose(reshape(collect(Iterators.flatten(bestp)),numParams,offset))
		count = count+1
	end


	#best = startingpt #mean(startingpt,1)
	best =mean(currparams,dims=1)
	data_dictionary=buildCompleteDictFromOneVector(best)
	names = data_dictionary["parameter_name_mapping_array"]
	str = ""
	j = 1
	for name in names
		currstr = string(name, " ", best[j]*.75, " ", best[j]*1.5, "\n")
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
	#startingpt =  readdlm("parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	startingpt = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	outputfn="../sensitivity/sobolboundspm50percentOnlyIC_02_05_18.txt"
	meanparams = startingpt#mean(startingpt,1)
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

function generateSobolParamsForOnlyPeturbFibrin()
	#let's use the averaged top 8 parameter sets
	startingpt =  readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	outputfn="../sensitivity/sobolboundspm50percentOnlyFibrin_05_12_17.txt"
	best = mean(startingpt,1)
	data_dictionary=buildCompleteDictFromOneVector(best)
	names = data_dictionary["parameter_name_mapping_array"]
	initial_conditions = data_dictionary["INITIAL_CONDITION_VECTOR"]
	IC_names = ["Fibrin", "Plasmin", "Fibrinogen", "Plasminogen", "tPA", "uPA", "Fibrin_Monomer", "Protofibril", "Antiplasmin", "PAI1", "Fiber"]
	str = ""
	j = 1
	selectedIdxs = vcat(collect(39:44), collect(48:77)) #params related to fibrin generation and degredation
	for name in names
		if j in selectedIdxs
			currstr = string(name, " ", best[j]*.5, " ", best[j]*1.5, "\n")
			str = string(str,currstr)
		end
		j=j+1
	end

	j = 12
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
	#filestr1 = "../sensitivity/02_05_18_MetricsForSobolPM50PercentOnlyParamsN1000_" #change me
	#filestr1 = "../sensitivity/01_10_18_MetricsForSobolPM50PercentOnlyParamsN_"
	filestr1 = "../sensitivity/05_10_18_MetricsForSobolPM50PercentOnlyICN100_"
	filestr2="_of_4.txt"
	str = ""
	for j in collect(1:4)
		fn = string(filestr1, j, filestr2)
		#currstr = readstring(fn)
		currstr = replace(readstring(fn), ",", " ")
		@show count(c->c=='\n', currstr)
		str = string(str, currstr)
	end
	write("../sensitivity/AllSobol_05_10_18_MetricsForSobolPM50PercentOnlyICN100.txt", str) #and me
end

