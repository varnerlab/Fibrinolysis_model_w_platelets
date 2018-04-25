include("Kinetics.jl")
include("Control.jl")

@everywhere function Balances(t,x,data_dictionary)
	idx_small = find(x.<0)
	x[idx_small] = 0.0
	#call the kinetics function
	rate_array = Kinetics(t,x,data_dictionary)
	#call the control function
	control_array = Control(t,x,data_dictionary)

	#compute the modified rate_array
	rate_array = rate_array.*control_array
	# calculate dxdt_reaction -
	dxdt_total = similar(x)
	dxdt_total[1] = -1*rate_array[2] -rate_array[7]-rate_array[6]			# 1 FII
	dxdt_total[2] = rate_array[2] - rate_array[4]+rate_array[7]+rate_array[6]	 # 2 FIIa
	dxdt_total[3] = -1*rate_array[3] 						# 3 PC
	dxdt_total[4] = 1*rate_array[3]  						# 4 APC
	dxdt_total[5] = -rate_array[4]							# 5 ATIII
	dxdt_total[6] = 0.0 							# 6 TM (acts as enzyme, so no change)
	dxdt_total[7] = -0.0							# 7 TRIGGER
	dxdt_total[8] = rate_array[22]-rate_array[23]#-rate_array[24]		 #8 frac active platelets
	dxdt_total[9] = -1*rate_array[1] 					#9 FV_FX
	dxdt_total[10] = rate_array[1] -rate_array[5]				#10 FV_FXa
	dxdt_total[11] = rate_array[5] 					#11 Prothromibase_platelets
	dxdt_total[12] = rate_array[12]-rate_array[15]                     #12 Fibrin
  	dxdt_total[13] = rate_array[13]+rate_array[14]-rate_array[16]     #13 Plasmin
	dxdt_total[14] = -rate_array[8]                 # 14 Fibrinogen
	dxdt_total[15] = -rate_array[13] - rate_array[14]                # 15 Plasminogen
	dxdt_total[16] = -rate_array[18]                                 # 16 tPA
	dxdt_total[17] = -rate_array[19]                                 # 17 uPA
	dxdt_total[18] = 1.0(rate_array[8]-rate_array[9]- rate_array[21])                 # 18 Fibrin monomer
	dxdt_total[19] = rate_array[9]+rate_array[11]-rate_array[10]       # 19 Protofibril
	dxdt_total[20] = -rate_array[16]                                 # 20 Antiplasmin
	dxdt_total[21] = -rate_array[17]-rate_array[18]-rate_array[19]   # 21 PAI_1
	dxdt_total[22] = 1.0*(rate_array[10]-rate_array[12]-rate_array[20])# 22 Fibers

	idx = find(x->(x<0),x);
	x[idx] = 0.0;
	dxdt_total[idx]= 0.0


	#adjust for mixing delay
	FIIa = x[2]
	#timing
	timing = data_dictionary["TIME_DELAY"]
	time_delay = timing[1]
	time_coeff = timing[2]
	aleph = data_dictionary["ALEPH"]

	tau = time_coeff*(1-FIIa/aleph)
	time_scale =1-1*exp(-tau*(t-time_delay))

	if(t<time_delay)
		time_scale = 0.0
		#@show t, time_scale, time_delaye
	end

	dxdt_total = dxdt_total.*time_scale

	return dxdt_total
end
