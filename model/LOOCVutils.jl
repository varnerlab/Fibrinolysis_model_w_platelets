function checkForPreviousAndLoad(outputfilestr)
	outputdir = "../LOOCV/"
	numObjs = 7
	allfiles = readdir(outputdir)
	#find out if we already have a file by this output name
	endfilestr = outputfilestr[10:end] #chop off first part
	f=[occursin(i,endfilestr) for i in allfiles]
	usefulfiles = allfiles[f]
	#if we don't find anything with the name, we start with the original guesses
	if(maximum(size(usefulfiles))==0)
		#initial_parameter_array = vec(readdlm("../parameterEstimation/startingPoint_02_05_18.txt"))
		#pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
		#allparams = readdlm(pathToParams, '\t')
		#allparams = readdlm("../parameterEstimation/Best2PerObj_09_04_18.txt")
		#allparams = readdlm("../parameterEstimation/Best8Overall_03_05_18.txt")
		allparams = readdlm("../parameterEstimation/Best2PerObjectiveParameters_05_12_18PlateletContributionToROTEM.txt")
		#allparams =readdlm("../parameterEstimation/Best2PerObjectiveParameters_01_02_19PlateletContributionToROTEM.txt")
		initial_parameter_array = allparams[1,:]
		round = 1
		return initial_parameter_array,round
	end
	#otherwise, we need to find the most recent file and use that
	editTimes = zeros(Float64,1)
	for f in usefulfiles
		#@show f
		currTime = Base.stat(string(outputdir,f)).mtime #get the time this file was last edited
		@show currTime
		#@show typeof(currTime)
		push!(editTimes,currTime)
	end
	editTimes = editTimes[2:end] #remove first element
	#this should be the index of the most recently edited file
	@show editTimes
	res = findmax(editTimes)
	@show res
	newestidx= res[2]
	@show usefulfiles
	startpointfn = usefulfiles[newestidx]
	@show startpointfn
	index1temp = search(startpointfn, "_", 36)
	index2temp = search(startpointfn, ".")
	round = parse(Int64, startpointfn[index1temp[1]+1:index2temp[1]-1])+1
	#we should now have the file name of the most recently edited file, so now we need to get the parameters out of it
	currPOETSpath = string("../LOOCV/", startpointfn)
	#something going wrong here-will fix
	ec,pc,ra = parsePOETsoutput(currPOETSpath, numObjs)
	#pull out the best overall parameters-n = 1
	bestp = generateBestNparameters(1, ec, pc)
	bestp = bestp[1]
	@show typeof(bestp)
	bestp = vec(bestp)
	@show bestp
	#bestp = convert(Array{Float64,1},bestp)
	@show typeof(bestp)
	return bestp,round
end

# Generates new parameter array, given current array -
function neighbor_function(parameter_array)
	outputfile="parameterEstimation/POETS_28_03_2017.txt"
	#@show size(parameter_array)
#	f = open(outputfile, "a")
#	write(f,string(parameter_array, "\n"))
#	close(f)

  SIGMA = 0.05
  number_of_parameters = length(parameter_array)

  # calculate new parameters -
  new_parameter_array = parameter_array.*(fill(1, size(parameter_array))+SIGMA*randn(number_of_parameters))

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

function acceptance_probability_function(rank_array,temperature)
  return (exp(-rank_array[end]/temperature))
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
	@show temperature
  return alpha*temperature
end


# Helper functions -
function parameter_bounds_function(parameter_array,lower_bound_array,upper_bound_array)

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

function generateNbestPerObjectiveAndErrors(n,ec_array, pc_array,savestring)
	num_objectives =size(ec_array,1)
	best_params=fill(Float64[], num_objectives*n)
	#best_params=Array{Array, (num_objectives*n)}
	best_errors =zeros(num_objectives, n)
	counter = 1
	for j in collect(1:num_objectives)
		curr_error=ec_array[j,:]
		allidx = collect(1:size(pc_array,2))
		removed = allidx
		for k in collect(1:n)
			min_index = argmin(curr_error)
			curr_best_params = pc_array[:,min_index]
			best_params[counter] = curr_best_params
			best_errors[j,k] = curr_error[min_index]
			removed = deleteat!(removed, min_index)
			#@show min_index, ec_array[:,min_index]
			#@show curr_best_params
			#delete the best ones we've found
			#@show size(total_error)
			#pc_array=pc_array[:,removed]
			curr_error=deleteat!(vec(curr_error),min_index)
			#@show size(pc_array)
			counter=counter+1
		end
	end	
	writedlm(savestring, best_params)
	return best_params,best_errors
end

function testLeftOutCase()
	#After running LOOCV, use me to test the parameter found on the other 7 cases on the 8th left out case
	numObjs = 7 #we had 7 different objectives
	numPatients = 8 #four patients used for test/validation
	numPerObj = 2 #let's pick the best 2 parameter sets per objective
	numCases = 1
	#change me to point to the correct files
	#outputstr = "../LOOCV/POETS_info_04_10_18_PlateletContributionToROTEMFlatness1ToBeTestedOn" #both tPa = 2 and tPa =0
	outputstr = "../LOOCV/POETS_info_31_10_18_PlateletContributionToROTEMFlatness1ToBeTestedOn" #only tPA 2
	allids = collect(9:16)
	countids = collect(1:8)
	countids2 = collect(9:16)

	#create storage. Fill with NaN to drop for means
	allerrors_train = NaN*ones(numPatients*numCases, numPerObj*numPatients)
	allerrors_test = NaN*ones(numPatients*numCases, numPerObj*numObjs)
	count = 1
	for id in allids
		println("Testing id=", id)
		currstr = string(outputstr, id,"Round_1.txt")
		#get results from LOOCV
		ec,pc,ra = parsePOETsoutput(currstr, numObjs)
		bestp, besterrs = generateNbestPerObjectiveAndErrors(numPerObj, ec, pc, "temp.txt")
		#store the errors
		upridx = count*numPerObj
		lwridx = (count-1)*numPerObj+1
		#adjust for 1 indexing
		if(lwridx ==0)
			lwridx =1
		end
		@show lwridx, upridx
		@show size(besterrs)
		leaveoutidx = count
		#trainidxs = [countids[1:leaveoutidx-1] ;countids[ leaveoutidx+1:end];countids2[1:leaveoutidx-1] ;countids2[ leaveoutidx+1:end]]
		trainidxs = [countids[1:leaveoutidx-1] ;countids[ leaveoutidx+1:end]]
		@show trainidxs
		allerrors_train[trainidxs,lwridx:upridx]= besterrs 

		#run the model with the best parameters
		#let's also save the timecourse data to file so we can plot it later
		@show size(bestp)
		for j in 1:size(bestp,1)
			println("On param set ", j, " out of ", size(bestp,1))
			temp_params = bestp[j]
			temp_params[47] = all_platelets[id] #set platelets to experimental value
			dict = buildCompleteDictFromOneVector(temp_params)
			initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
			#to adjust per patient. Comment out if we don't want to do that
			#initial_condition_vector=setCompleteModelIC(initial_condition_vector,id)
			#let's only consider tPA = 2 case for now
			tPA = 2.0
			initial_condition_vector[16]=tPA
			TSTART = 0.0
			Ts = .02
			TSTOP = 90.0
			#solve
			TSIM = collect(TSTART:Ts:TSTOP)
			fbalances(t,y)= Balances(t,y,dict) 
			#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
			t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0,points=:specified)
			FIIa = [a[2] for a in X]
			fibrinogen = [a[14] for a in X]
			A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[id])
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[id])
			allerrors_test[count,j]=MSE
			#store time series to disk so we can average and plot it
			writedlm(string("../LOOCV/validation/TrainsimulatedROTEMPatient", id,"paramset", j, "tPA=2.txt"),hcat(t,A))


			#let's run the tPA = 0 case, too
			tPA = 0.0
			initial_condition_vector[16]=tPA
			TSTART = 0.0
			Ts = .02
			TSTOP = 120.0
			#solve
			TSIM = collect(TSTART:Ts:TSTOP)
			fbalances(t,y)= Balances(t,y,dict) 
			#t,X = ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep=1E-9)
			t,X=ODE.ode23s(fbalances,(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.0,points=:specified)
			FIIa = [a[2] for a in X]
			fibrinogen = [a[14] for a in X]
			A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[id])
			#store time series to disk so we can average and plot it
			writedlm(string("../LOOCV/validation/TrainsimulatedROTEMPatient", id,"paramset", j, "tPA=0.txt"),hcat(t,A))
		end
		count = count+1	
	end
	return allerrors_test, allerrors_train
	#note to self: to get means out, since you have NaN, use the images package (using Images)
	#then meanfinite(allerrors_test,2)
	#meanfinite(allerrors_train,2)

	
end

function plotLeftOutCase()
	#after calculating errors, use this to plot the left out case
	numObjs = 7 #we had 7 different objectives
	numPatients = 8
	numPerObj = 2 #let's pick the best 2 parameter sets per objective
	#change me to point to the correct files
	ts_str = "../LOOCV/validation/TrainsimulatedROTEMPatient"
	allids = collect(9:16)
	countids = collect(1:8)
	#so we know dimensions
	TSTART = 0.0
	Ts = .02
	TSTOP = 90.0
	TSIM2 = collect(TSTART:Ts:TSTOP)
	TSIM0 = collect(TSTART:Ts:120.0)
	for id in allids
		#need to create array to store all of the ROTEM curves to average
		all_A2 = zeros(numObjs*numPerObj, maximum(size(TSIM2)))
		all_A0 = zeros(numObjs*numPerObj, maximum(size(TSIM0)))
		for j in 1:numObjs*numPerObj
			#read in data
			currstr = string(ts_str, id, "paramset", j, "tPA=0.txt")
			currdat = readdlm(currstr)
			t = currdat[:,1]
			A = currdat[:,2]
			all_A0[j, :]=A
			currstr = string(ts_str, id, "paramset", j, "tPA=2.txt")
			currdat = readdlm(currstr)
			t = currdat[:,1]
			A = currdat[:,2]
			all_A2[j, :]=A
		end
		#now, plot
		close("all")
		figure(figsize = (10,10))
		#plot tPA =0 case
		plot(allexperimentaldata[id-8][:,1], allexperimentaldata[id-8][:,2], "k", linewidth = 3)
		@show size(mean(all_A0,1))
		mean_arr = mean(all_A0,1)'
		plot(TSIM0, mean_arr, color="k", linewidth =.8)
		std_arr = std(all_A0, 1)'
		SF = 1.96/sqrt(numObjs*numPerObj)
		UB = mean_arr + SF*std_arr
		LB = mean_arr - SF*std_arr
		#to plot errors
		fill_between(TSIM0, vec(LB), vec(UB), color ="silver", alpha = .5)
		#plot tPA =2 case
		plot(allexperimentaldata[id][:,1], allexperimentaldata[id][:,2], "lightskyblue", linewidth = 3)
		mean_arr = mean(all_A2,1)'
		plot(TSIM2, mean_arr, color="lightskyblue", linewidth =.8)
		std_arr = std(all_A2, 1)'
		SF = 1.96/sqrt(numObjs*numPerObj)
		UB = mean_arr + SF*std_arr
		LB = mean_arr - SF*std_arr
		#to plot errors
		fill_between(TSIM2, vec(LB), vec(UB), color ="lightsteelblue", alpha = .5)
		xlabel("Time (minutes)", fontsize = 36)
		ylabel("Amplitude (mm)", fontsize=36)

		savefig(string("../LOOCV/validation/TestingOnID", id, "10_31_18.pdf"))
		
	end
end

function testLeftOutCaseValidation()
	#After running LOOCV, use me to test the parameter found on the other 7 cases on the 8th left out case
	numObjs = 3 #we had 3 different objectives
	numPatients = 4
	numPerObj = 2 #let's pick the best 2 parameter sets per objective
	numCases = 1
	#change me to point to the correct files
	#outputstr = "../LOOCV/POETS_info_04_10_18_PlateletContributionToROTEMFlatness1ToBeTestedOn" #both tPa = 2 and tPa =0
	outputstr = "../LOOCV/POETS_info_12_02_19_PlateletContributionToROTEMFlatness1ToBeTestedOn" #only tPA 2
	allids =collect(11:14)
	countids = collect(1:4)
	countids2 = collect(9:16)

	#create storage. Fill with NaN to drop for means
	allerrors_train = NaN*ones(numPatients*numCases, numPerObj*numPatients)
	allerrors_validate = NaN*ones(numPatients*numCases, numPerObj*numObjs)
	count = 1
	for id in allids
		println("Testing id=", id)
		currstr = string(outputstr, id,"Round_1Max_iters2.txt")
		#get results from LOOCV
		ec,pc,ra = parsePOETsoutput(currstr, numObjs)
		bestp, besterrs = generateNbestPerObjectiveAndErrors(numPerObj, ec, pc, "temp.txt")
		#store the errors
		upridx = count*numPerObj
		lwridx = (count-1)*numPerObj+1
		#adjust for 1 indexing
		if(lwridx ==0)
			lwridx =1
		end
		@show lwridx, upridx
		@show size(besterrs)
		leaveoutidx = count
		#trainidxs = [countids[1:leaveoutidx-1] ;countids[ leaveoutidx+1:end];countids2[1:leaveoutidx-1] ;countids2[ leaveoutidx+1:end]]
		trainidxs = [countids[1:leaveoutidx-1] ;countids[ leaveoutidx+1:end]]
		@show trainidxs
		@show size(allerrors_train[trainidxs,lwridx:upridx])
		allerrors_train[trainidxs,lwridx:upridx]= besterrs 

		#run the model with the best parameters
		#let's also save the timecourse data to file so we can plot it later
		@show size(bestp)
		for j in 1:size(bestp,1)
			println("On param set ", j, " out of ", size(bestp,1))
			temp_params = bestp[j]
			temp_params[47] = all_platelets[id] #set platelets to experimental value
			dict = buildCompleteDictFromOneVector(temp_params)
			initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
			#to adjust per patient. Comment out if we don't want to do that
			#initial_condition_vector=setCompleteModelIC(initial_condition_vector,id)
			#let's only consider tPA = 2 case for now
			tPA = 2.0
			initial_condition_vector[16]=tPA
			TSTART = 0.0
			Ts = .02
			TSTOP = 90.0
			#solve
			TSIM = collect(TSTART:Ts:TSTOP)
			fbalances(y,p,t)= Balances(t,y,dict) 
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			@time sol = solve(prob)
			t =sol.t
			X = sol
			FIIa = X[2,:]
			fibrinogen = X[14,:]
			A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[id])
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[id])
			allerrors_validate[count,j]=MSE
			#store time series to disk so we can average and plot it
			writedlm(string("../LOOCV/validation/ValidatesimulatedROTEMPatient", id,"paramset", j, "tPA=2.txt"),hcat(t,A))


			#let's run the tPA = 0 case, too
			tPA = 0.0
			initial_condition_vector[16]=tPA
			TSTART = 0.0
			Ts = .02
			TSTOP = 120.0
			#solve
			TSIM = collect(TSTART:Ts:TSTOP)
			fbalances(y,p,t)= Balances(t,y,dict) 
			#fbalances(t,y)= Balances(t,y,dict) 
			#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			@time sol = solve(prob)
			t =sol.t
			X = sol
			A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[id])
			#store time series to disk so we can average and plot it
			writedlm(string("../LOOCV/validation/ValidatesimulatedROTEMPatient", id,"paramset", j, "tPA=0.txt"),hcat(t,A))
		end
		count = count+1	
	end
	return allerrors_validate, allerrors_train
	#note to self: to get means out, since you have NaN, use the images package (using Images)
	#then meanfinite(allerrors_test,2)
	#meanfinite(allerrors_train,2)
	
end

function testLeftOutCaseValidation(batch_number::Int64)
	#After running LOOCV, use me to test the parameter found on the other 7 cases on the 8th left out case
	numObjs = 3 #we had 3 different objectives
	numPatients = 4
	numPerObj = 2 #let's pick the best 2 parameter sets per objective
	numCases = 1
	#change me to point to the correct files
	#outputstr = "../LOOCV/POETS_info_04_10_18_PlateletContributionToROTEMFlatness1ToBeTestedOn" #both tPa = 2 and tPa =0
	outputstr = "../LOOCV/POETS_info_14_02_19_PlateletContributionToROTEMFlatness1ToBeTestedOn" #only tPA 2
	allids =collect(11:14)
	countids = collect(1:4)
	countids2 = collect(9:16)

	#create storage. Fill with NaN to drop for means
	allerrors_train = NaN*ones(numPatients*numCases, numPerObj*numPatients)
	allerrors_validate =  NaN*ones(numPatients*numCases, numPerObj*numObjs)
	count = 1
	for id in allids
		println("Testing id=", id)
		currstr = string(outputstr, id,"Round_1Max_iters10.txtRound_", batch_number, "outOf10.txt")
		#get results from LOOCV
		ec,pc,ra = parsePOETsoutput(currstr, numObjs)
		bestp, besterrs = generateNbestPerObjectiveAndErrors(numPerObj, ec, pc, "temp.txt")
		#store the errors
		upridx = count*numPerObj
		lwridx = (count-1)*numPerObj+1
		#adjust for 1 indexing
		if(lwridx ==0)
			lwridx =1
		end
		#@show lwridx, upridx
		#@show size(besterrs)
		leaveoutidx = count
		#trainidxs = [countids[1:leaveoutidx-1] ;countids[ leaveoutidx+1:end];countids2[1:leaveoutidx-1] ;countids2[ leaveoutidx+1:end]]
		trainidxs = [countids[1:leaveoutidx-1] ;countids[ leaveoutidx+1:end]]
		#@show trainidxs
		#@show size(allerrors_train[trainidxs,lwridx:upridx])
		allerrors_train[trainidxs,lwridx:upridx]= besterrs 

		#run the model with the best parameters
		#let's also save the timecourse data to file so we can plot it later
		#@show size(bestp)
		for j in 1:size(bestp,1)
			println("On param set ", j, " out of ", size(bestp,1))
			temp_params = bestp[j]
			temp_params[47] = all_platelets[id] #set platelets to experimental value
			dict = buildCompleteDictFromOneVector(temp_params)
			initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
			#to adjust per patient. Comment out if we don't want to do that
			#initial_condition_vector=setCompleteModelIC(initial_condition_vector,id)
			#let's only consider tPA = 2 case for now
			tPA = 2.0
			initial_condition_vector[16]=tPA
			TSTART = 0.0
			Ts = .02
			TSTOP = 90.0
			#solve
			TSIM = collect(TSTART:Ts:TSTOP)
			fbalances(y,p,t)= Balances(t,y,dict) 
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			sol = solve(prob)
			t =sol.t
			X = sol
			FIIa = X[2,:]
			fibrinogen = X[14,:]
			A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[id])
			MSE, interpData = calculateMSE(t,A, allexperimentaldata[id])
			allerrors_validate[count,j]=MSE
			#store time series to disk so we can average and plot it
			writedlm(string("../LOOCV/validation/ValidatesimulatedROTEMPatient", id,"paramset", j, "tPA=2.txt"),hcat(t,A))


			#let's run the tPA = 0 case, too
			tPA = 0.0
			initial_condition_vector[16]=tPA
			TSTART = 0.0
			Ts = .02
			TSTOP = 120.0
			#solve
			TSIM = collect(TSTART:Ts:TSTOP)
			fbalances(y,p,t)= Balances(t,y,dict) 
			#fbalances(t,y)= Balances(t,y,dict) 
			#t,X=ODE.ode23s(fbalances,vec(initial_condition_vector),TSIM, abstol = 1E-6, reltol = 1E-6, minstep = 1E-8,maxstep = 1.00)
			prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
			sol = solve(prob)
			t =sol.t
			X = sol
			A = convertToROTEMPlateletContribution(t,X,tPA,all_platelets[id])
			#store time series to disk so we can average and plot it
			writedlm(string("../LOOCV/validation/ValidatesimulatedROTEMPatient", id,"paramset", j, "tPA=0.txt"),hcat(t,A))
		end
		count = count+1	
	end
	return allerrors_validate, allerrors_train
	#note to self: to get means out, since you have NaN, use the images package (using Images)
	#then meanfinite(allerrors_test,2)
	#meanfinite(allerrors_train,2)
	
end

function testAllBatches()
	close("all")
	numBatches =10
	mean_val_err = []
	mean_train_err = []
	std_val_err =[]
	std_train_err = []
	for j =1:numBatches
		println(string("On batch", j))
		allerrors_validate, allerrors_train= testLeftOutCaseValidation(j)
		@show allerrors_validate
		@show allerrors_train
		push!(mean_val_err, mean(allerrors_validate))
		push!(mean_train_err, mean(meanfinite(allerrors_train,2)))
		push!(std_val_err, std(allerrors_validate))
		push!(std_train_err, std(meanfinite(allerrors_train,2)))
		@show mean_train_err
		@show mean_val_err
	end
	writedlm("../LOOCV/compareTest_Validation/Train_mean_14_02_18.txt", mean_train_err)
	writedlm("../LOOCV/compareTest_Validation/Val_mean_14_02_18.txt", mean_val_err)
	writedlm("../LOOCV/compareTest_Validation/Train_std_14_02_18.txt", std_train_err)
	writedlm("../LOOCV/compareTest_Validation/Val_std_14_02_18.txt", std_val_err)
	figure(figsize = [10,10])
	batches = 1:numBatches
	@show std_val_err
	@show std_train_err
	errorbar(batches, mean_val_err,yerr=std_val_err)
	errorbar(batches, mean_train_err, yerr=std_train_err)
	ax = gca()
	ax[:tick_params]("both",labelsize=24) 
	legend(["Validation Error", "Training Error"])
	ylabel("Average Error", fontsize=36)
	xlabel("Number of POETs Iterations",fontsize=36)
	savefig("../figures/ConvergenceCurves_Feb14.pdf")
end

function testAllBatchesAndMakeFigures()
	close("all")
	numBatches =10
	numObjs = 3 #we had 3 different objectives
	numPatients = 4
	numPerObj = 2 #let's pick the best 2 parameter sets per objective
	numCases = 1
	allids =collect(11:14)
	outputstr = "../LOOCV/POETS_info_14_02_19_PlateletContributionToROTEMFlatness1ToBeTestedOn" #only tPA 2
	numParams = 77
	#for j =1:numBatches
	for j =10:numBatches
		println(string("On batch ", j))
		currparams = zeros(numPerObj*numObjs*size(allids,1), numParams)
		count = 1
		offset = numPerObj*numObjs
		batch_number=j
		for id in allids				
			currstr = string(outputstr, id,"Round_1Max_iters10.txtRound_", batch_number, "outOf10.txt")
			#get results from LOOCV
			ec,pc,ra = parsePOETsoutput(currstr, numObjs)
			bestp, besterrs = generateNbestPerObjectiveAndErrors(numPerObj, ec, pc, "temp.txt")
			currparams[(count-1)*offset+1:count*offset,:]=transpose(reshape(collect(Iterators.flatten(bestp)),numParams,offset))
			count = count+1
		end
		writedlm(string("../LOOCV/bestparamsForBatch_", batch_number, "_14_02_19.txt"),currparams)
		makePredictionFigurePlatletContributionToROTEM(currparams, string("../figures/LOOCV_14_02_19_Predicitions_after_batch", j, "final.pdf"),string("../figures/LOOCV_14_02_19_Predicitions_after_batch", j, "subfig.pdf"))
		makeTrainingFigurePlatletContributionToROTEM(currparams, string("../figures/LOOCV_14_02_19_Training_after_batch", j, "final.pdf"),string("../figures/LOOCV_14_02_19_Training_after_batch", j, "subfig.pdf"))
	end
end

function makeTrainingFigurePlatletContributionToROTEM(best_params, savestr_final,savestr)
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")

	ids = [5,6,7,8]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	numParamSets = size(best_params,1)
	adjustICs = false
	for j in collect(1:size(ids,1))
		@show ids[j]
		for k in collect(1:size(tPAs,1))
			#bestparams=generateBestNparameters(8,ec,pc)
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParamsPlatetContributionToROTEM(best_params, ids[j], tPAs[k], savestr,adjustICs)
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

	savefig(savestr_final)
end

function makePredictionFigurePlatletContributionToROTEM(best_params, savestr_final,savestr)
	#use me
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	ids = [3,4,9,10]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	numParamSets = size(best_params,1)
	adjustICs=false
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParamsPlatetContributionToROTEM(best_params, ids[j], tPAs[k], savestr,adjustICs)
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

	savefig(savestr_final)
end

function plotTrainValErr()
	close("all")
	mean_train_err=readdlm("../LOOCV/compareTest_Validation/Train_14_02_18.txt")
	mean_val_err=readdlm("../LOOCV/compareTest_Validation/Val_14_02_18.txt")
	numBatches=size(mean_train_err,1)
	figure(figsize = [10,10])
	batches = 1:numBatches
	plot(batches, mean_val_err, "kx", markersize=13)
	plot(batches, mean_train_err, "ko", markersize=13)
	ax = gca()
	ax[:tick_params]("both",labelsize=24) 
	legend(["Validation Error", "Training Error"], fontsize=30)
	ylabel("Average Error", fontsize=36)
	xlabel("Number of POETs Iterations",fontsize=36)
	savefig("../figures/ConvergenceCurves_Feb14_nobars.pdf")
end
