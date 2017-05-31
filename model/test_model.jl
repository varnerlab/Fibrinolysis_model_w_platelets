#estimate parameters
using POETs
include("objective_function.jl")

function test_model(initial_parameter_array, selected_idxs)
	outputfile = string("parameterEstimation/LOOCVSavingAllParams_2017_02_06_low_acc/POETS_selectedindices",string(selected_idxs) ,".txt")
	number_of_subdivisions = 10
	number_of_parameters = 46
	number_of_objectives = 5
	ec_array = zeros(number_of_objectives)
  	pc_array = zeros(number_of_parameters)
	for index in collect(1:number_of_subdivisions)

	    # Run JuPOETs -
		objective_function_fed(initial_parameter_array) = objective_function(initial_parameter_array, selected_idxs)
	    (EC,PC,RA) = estimate_ensemble(objective_function_fed,neighbor_function,acceptance_probability_function,cooling_function,(vec(initial_parameter_array));rank_cutoff=4,maximum_number_of_iterations=10,show_trace=true)
		@show (EC, PC, RA)
		# Package -
		ec_array = [ec_array EC]
		pc_array = [pc_array PC]
		#ra_array = [ra_array RA]
	
		#write outputs to file
		@show (EC, PC, RA)
		f = open(outputfile, "a")
		write(f, string(EC, ",", PC, ",", RA, "\n"))
		close(f)
	  end

  return (ec_array,pc_array)

end

function doLOOCV()
	initial_parameter_array = readdlm("parameterEstimation/handfit.txt", ',')
	experimental_indices = collect(1:6)
	for j in collect(1: size(experimental_indices,1))
		experimental_indices = collect(1:6)
		selected_idxs = deleteat!(experimental_indices,j)
		ec_array, pc_array=test_model(initial_parameter_array, selected_idxs)
		#create new initial parameter array best on best for this objective-take the set that gives lowest error on all 5 objectives
		initial_parameter_array_top10 = generateBestNparameters(10,ec_array, pc_array)
		bestparams_output = string("parameterEstimation/LOOCVSavingAllParams_2017_02_06_low_acc/bestParamSetsFromLOOCV",j, "excluded.txt")
		f = open(bestparams_output, "a+")
		writedlm(f,initial_parameter_array_top10, ',')
		close(f)
		initial_parameter_array = initial_parameter_array_top10[1]
	end
end

function regenerateBestParams()
	experimental_indices = collect(1:6)
	for j in collect(1:6)
		bestparams_output = string("parameterEstimation/LOOCVSavingAllParams_2017_1_31_No_specific_bounds/bestParamSetsFromLOOCV",j, "excluded.txt")
		experimental_indices = collect(1:6)
		selected_idxs = deleteat!(experimental_indices,j)
		outputfile = string("parameterEstimation/LOOCVSavingAllParams_2017_1_31_No_specific_bounds/POETS_selectedindices",string(selected_idxs) ,".txt")
		ec, pc, ra = parsePOETsoutput(outputfile)
		initial_parameter_array_top10 = generateBestNparameters(10,ec, pc)
		bestparams_output = string("parameterEstimation/LOOCVSavingAllParams_2017_1_31_No_specific_bounds/bestParamSetsFromLOOCV",j, "excluded.txt")
		f = open(bestparams_output, "a+")
		writedlm(f,initial_parameter_array_top10, ',')
		close(f)
	end
end
