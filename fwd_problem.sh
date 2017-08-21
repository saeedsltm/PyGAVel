#!/bin/bash

# RUN HYPOELLIPSE
hypoell-loc.sh ${1}

# GET MEAN RMS
#gawk 'BEGIN { rmssum=0; neqs=0 }{ rmssum=rmssum+substr($0,48,4)/100.; neqs=neqs+1 } END { printf("%10.5f\n", rmssum/neqs ) }' ${1}.sum  > ${1}misfit.val
awk '/average rms of all events =/{print substr($0,29,11)}' ${1}.out > ${1}misfit.val
