from flashover import firecalc
import pdb
configlocation = 'modelconfig.yaml'
A = firecalc(configlocation, 'dalmarnock','dalmarnock',8)
print(A[0],A[1])
