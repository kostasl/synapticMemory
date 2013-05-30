#!/bin/sh
echo "Filter Size $1 (1) To $2 (2) Repetitions : $3 (3)"
echo "cAMP Decay $4 (4)"
echo "PKA AllocThreshold $5 (5)"
echo "TRIALS:  ${6:-1000}(6)"
echo "ThresholdCycles: ${7:-1}"
logc="Release/simSATOut$1-${2:?NoSize}R${3:?NoReps}F${4:?NoFc}P${5:?NoPKAThres}c${7:-1}.log"

echo "std output in $logc"
cd ../
nohup ./Release/simAllocation-Sat-PKA --simulation=AllocSignalVsRepetitionTime --model=synapseSingleFilterUnifiedWithDecay --trials=${6:-1000} --startSize=$1 --endSize=$2 --synapsesSize=10001 --initPeriod=0 --cSimTimeSecs=1000 --Timestep=1.00 --repPatIndex=0 --repPatCount=$3 --AllocRefrac=${7:-1} --metaSampleTime=0 --cAMPDecay=$4 --PKAAllocThres=$5 > $logc 2>&1 &

echo "execute tail -f $log"

