using PyPlot
using POETs

function parsePOETsoutputmakegraphs(filename)
	close("all")
	f = open(filename)
	alltext = readall(f)
	close(f)

	outputname = "parameterEstimation/textparsing.txt"
	number_of_parameters = 46
  	number_of_objectives = 5
	ec_array = zeros(number_of_objectives)
  	pc_array = zeros(number_of_parameters)
	rank_array = zeros(1)	
	counter =1
	for grouping in matchall(r"\[([^]]+)\]", alltext)
		cleanedgrouping = replace(grouping, "[", "")
		nocommas = replace(cleanedgrouping, ","," ")
		allcleaned = replace(nocommas, "]", "")
		outfile = open(outputname, "w")
		write(outfile, allcleaned)
		close(outfile)
		formatted = readdlm(outputname)
#		@show formatted	
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
	#@show pc_array

#	@show size(ec_array)
#	@show size(pc_array)
#	@show size(rank_array)

	#@show ec_array[:,2:end]
	#rerank
	rerankedarray=rank_function(ec_array[:,2:end])
	#@show size(rerankedarray)
	#make trade off surface
	#plt.rc("font", family="")
#	figure()
#	PyCall.PyDict(matplotlib["rcParams"])["font.sans-serif"] = ["Helvetica"]
#	hold("on")
#	#@show size(ec_array,2)
#	for j in collect(2:size(ec_array,2))
#		if(rank_array[j]==0)
#			plot(ec_array[1,j], ec_array[2,j], linewidth = .3,"ko", markersize = 2.5,markeredgewidth=0.0)
#		else
#			plot(ec_array[1,j], ec_array[2,j], linewidth = .3,"o", color = ".75", markersize = 2.5,markeredgewidth=0.0)
#		end
#	end
#	xlabel("Objective 2", fontsize=18)
#	ylabel("Objective 1", fontsize=18)
#	savefig("parameterEstimation/TradeOffCurveNoPeturb.png")
#	generateBestNparameters(5, ec_array, pc_array)
	return ec_array[:,2:end], rerankedarray, pc_array[:,2:end]
end

function generateBestNparameters(n, ec_array, pc_array)

	#calculate error
	best_params = Array[]
	total_error = sum(ec_array[:,2:end],1)
	for k in collect(1:n)
		min_index = indmin(total_error)
		curr_best_params = pc_array[:,min_index]
		push!(best_params, curr_best_params)
		@show min_index
		@show curr_best_params
		#delete the best ones we've found
		pc_array[1:size(pc_array,1) .!= min_index,: ]
		deleteat!(vec(total_error),min_index)
	end
	writedlm(string(n, "ParameterSetsNoPeturb.txt"), best_params)

end
