#!/bin/sh

cd /fs/home/ril34/Documents/FibrinolysisWPlatelets/model
 #if we're not in the right place, move there
#for i in `seq 1 25`;
for i in `seq 1 50`
do 
	echo "#!/bin/sh" >> currjob.nbs
	echo '#PBS -M ril34@cornell.edu' >> currjob.nbs
	echo '##NBS-email: ril34@cornell.edu' >> currjob.nbs #for getting update emails
	echo '##NBS-unique: yes'>> currjob.nbs
	echo '##NBS-queue: medium'>> currjob.nbs
	echo '#PBS -p 1024' >> currjob.nbs
	echo 'cd /fs/home/ril34/Documents/FibrinolysisWPlatelets/model/' >>currjob.nbs

	echo $'~/./Documents/julia/bin/julia -p 1 solveInverseProbGenerateProbSA.jl' `$(i)` ''>> currjob.nbs
	mv currjob.nbs submitFindICNominaliter$i.nbs
	jsub submitFindICNominaliter$i.nbs -nproc `
	rm submitFindICNominaliter$i.nbs
done

