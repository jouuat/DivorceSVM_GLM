# PARAMETERS
param n     # num of X
param col     # num of columns for each X
param v     # nu
param x{i in 1..n, j in 1..col}     #X
param y{i in 1..n}     #Y


# VARIABLES
var w{i in 1..col}    
var gamma     # intercept of hyperplan
var s{i in 1..n}     # slack variables (tantes com files)
var y_pred{i in 1..n}    
var errors    


# OBJECTIVE FUNCTION
minimize obj_func: 0.5*sum{i in 1..col} (w[i]^2) + v*sum{i in 1..n}s[i]    

# CONSTRAINTS
subject to set_stack{i in 1..n}: # for the points misclassified set s = 1
  y[i] * ( sum{j in 1..col} ( w[j] * x[i,j] ) + gamma) >= 1-s[i]    
subject to positiveness{i in 1..n}: #impose that the slack must be positive
  s[i]>=0    

  # DATA
  data    
  data instancesOTDM.dat    
  #data divorce_sample.dat    

  # SOLVE
  model    
  option solver './cplex'    
  solve    


  #AMNALIZE THE PREDICTION
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
