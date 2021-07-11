using MAT

file = matopen("mesh.mat")
varnames = names(file)
close(file)

print(varnames)

vars = matread("mesh.mat")

# print(vars)