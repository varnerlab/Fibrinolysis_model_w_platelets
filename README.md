# Fibrinolysis_model_w_platelets
## The Model
This model describes both coagulation and fibrinolysis, modeled in a logical rules/ODE frame work. The resulting ODEs are solved with the Julia DifferentialEquations package. 

## Dependencies
You will need DifferentialEquations, ExcelReaders, and PyPlot to run the model. JuPOETs is neccessary if you wish to estimate parameters in the multi-objective framework, and Optim is used to solve the inverse problem.

## Execution
The model can be run by calling runModelWithParams(params), which takes a 1 by 77 parameter vector. Many parameter estimates can be found in the parameterEstimation directory.
