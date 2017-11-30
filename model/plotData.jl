using PyPlot
using ExcelReaders

function plotAverageThrombin(pathToData, savestr)
	close("all")
	data = readdlm(pathToData)
	time = data[:,1]
	avg_run = mean(data[:,2:3],2);
	plot(time, avg_run, "k.");
	xlabel("Time in seconds");
	ylabel("Thrombin Concentration (micormolar)")
	savefig(savestr)
end

function plotAverageROTEM()
	close("all")
	pathToData="../data/Viscoelasticmeasurements.xlsx"
	p9caseWTPA="Timecourse!H3:J1166"
	p9caseNoTPA="Timecourse!K3:M1544"
	p10caseWTPA="Timecourse!A3:C1301"
	p10caseNoTPA="Timecourse!D3:F789"
	p6caseWTPA = "Timecourse!AC3:AE1097"
	p6caseNoTPA="Timecourse!AF3:AH988"
	cases = [p6caseWTPA,p6caseNoTPA,p9caseWTPA, p9caseNoTPA, p10caseWTPA, p10caseNoTPA]
	j = 1
	widths = [.5,.5, 1,1,2, 2]
	figure()
	for case in cases
		if(j%2 ==0)
			colorstr = "0.1"
		else
			colorstr=".5"
		end
		data=readxl(pathToData, case)
		data=Array{Float64}(data)
		time = data[:,1]
		avg_run = mean(data[:,2:3],2);
		calculateCommonMetrics(avg_run,time./60)
		plot(time, avg_run, "-", color= colorstr, linewidth = widths[j]);
		j = j+1
	end
	legend(["With tPA=2nM p6", "Without tPA p6","With tPA=2nM p9", "Without tPA p9","With tPA=2nM p10", "Without tPA p10"], loc="best")
	xlabel("Time in seconds");
	ylabel("ROTEM Curve")
	savefig("figures/ROTEMcurvesForP_9_10_6.pdf")
end

function plotAllROTEM()
	close("all")
	figure(figsize=[20,20])
	leg_arr=AbstractString[]
	poss_tPA = [0,2]
	colors = [".3",".7"]
	markers = ["x", ".", "o", "h", "+", "p", "3", "8"]
	ids = ["3", "4", "5", "6", "7", "8", "9", "10"]
	for j in collect(1:size(poss_tPA,1))
		for k in collect(1:size(ids,1))
			platelets,currdata = setROTEMIC(poss_tPA[j], ids[k])
			plot(currdata[collect(1:10:end),1], currdata[collect(1:10:end),2], color=colors[j], marker = markers[k])
			leg_str = string("Id = ", ids[k], " tPA = ", poss_tPA[j])
			push!(leg_arr, leg_str)
		end
	end
	xlabel("Time, in minutes")
	ylabel("Diameter")
	legend(leg_arr)
	savefig("figures/allROTEMcurves.pdf")
end

function plotROTEM_given_tPA(tPA)
	#close("all")
	ids = ["3", "4", "5", "6", "7", "8", "9", "10"]
	for k in collect(1:size(ids,1))
		platelets,currdata = setROTEMIC(tPA, ids[k])
		plot(currdata[:,1], currdata[:,2], "k")
	end
	xlabel("Time, in minutes")
	ylabel("Diameter")
	savefig(string("figures/allROTEMcurves_tPA=",tPA, ".pdf"))
end
