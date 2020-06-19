function maxDeviation = ppmDeviation(mz_ref,tolerance)

mz_max = mz_ref+(mz_ref*tolerance)/1000000;
maxDeviation = abs(mz_max-mz_ref);