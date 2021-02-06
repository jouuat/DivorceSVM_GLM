# MODEL SECTION
model    

# PARAMETERS
param n     # num of X
param col     # num of columns for each X
param v     # nu
param x{i in 1..n, j in 1..col}     #X
param y{i in 1..n}     #Y


# VARIABLE
var lambda {i in 1..n}    
var w{i in 1..col}    
#var gamma{i in 1..n}    
var gamma    
#var bin_s{i in 1..n} binary     # slack variables (tantes com files)
var y_pred{i in 1..n}    
var errors    


# OBJECTIVE FUNCTION (we will maximaze the Lagrange fucntion)
maximize obj_func :
  sum{i in 1..n} lambda[i] -
  0.5*sum{i in 1..n, j in 1..n}lambda[i]*lambda[j]*y[i]*y[j]*
  (sum{k in 1..col}x[i,k]*x[j,k])    

# CONSTRAINT
subject to constraint:
  sum{i in 1..n}lambda[i]*y[i]=0    
subject to lambda_constraint_1{i in 1..n}:
  lambda[i]>=0    
subject to lambda_constraint_2{i in 1..n}:
  lambda[i]<=v    
subject to retrieve_w{i in 1..col}:
  w[i]=sum{j in 1..n}lambda[j]*y[j]*x[j,i]    


# DATA
data    
data instancesOTDM.dat    
#data divorce_sample.dat    

# SOLVE
model    
option solver './cplex'    
solve    

for{i in 1..n} {
	if ((0.0000001) < lambda[i] < (v-0.0000001)) then let gamma := 1/y[i] - sum{j in 1..col}w[j]*x[i,j]    
}    


# ANALIZE THE PREDICTION
for {i in 1..n} {
	if (sum{j in 1..col}(x[i,j]*w[j]) + gamma >= 0) then let y_pred[i] := 1    
	if (sum{j in 1..col}(x[i,j]*w[j]) + gamma < 0) then let y_pred[i] := -1    
	}    
for {i in 1..n} {
  if (y[i]>=0) then let errors := errors + (y[i] - y_pred[i])/2    
	if (y[i]<0) then let errors := errors + (y_pred[i] - y[i])/2    
	}    

#DISPLAY

display w, gamma, y, y_pred, errors    
