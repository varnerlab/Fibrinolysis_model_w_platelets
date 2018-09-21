@everywhere function checkForPreviousAndLoad(outputfilestr)
	outputdir = "../LOOCV/"
	numObjs = 7
	allfiles = readdir(outputdir)
	#find out if we already have a file by this output name
	endfilestr = outputfilestr[10:end] #chop off first part
	f=[contains(i,endfilestr) for i in allfiles]
	usefulfiles = allfiles[f]
	#if we don't find anything with the name, we start with the original guesses
	if(maximum(size(usefulfiles))==0)
		initial_parameter_array = vec(readdlm("../parameterEstimation/startingPoint_02_05_18.txt"))
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
