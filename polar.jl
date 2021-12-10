# Thomas Simcox
# MATH 466 - Project 2
#
# polar decomposition algorithm

using LinearAlgebra
using Plots

A = [-0.49 -0.21 -0.40 -0.21;
     0.36 0.10 0.29 -0.04;
     0.12 -0.01 0.48 -0.47;
     0.09 0.09 -0.41 0.22]

# function that calculates next iteration
pol! = (A::Matrix{Float64}) -> (A + (inv(A))') / 2;

# driver function that finds the n'th iteration
function pol_driver(A::Matrix{Float64}, n::Int64)::Matrix{Float64}
    B = copy(A)
    for i = 0:n
        B = pol!(B)
    end
    return B;
end

# test the algorithm
X_0 = copy(A);
X_1 = pol_driver(X_0, 0);
X_test = norm(X_0 - X_1);
println("||X_0 - X_1|| = $X_test");

# compute delta_k 
delta_k = (X, k) -> pol_driver(X, k) - X;

# store all the delta_k's in a list for k = 1:9
delta_ks = [delta_k(copy(A), k) for k = 0:9];

# display the norm of each delta_k
map(x -> "||Delta_k|| = $(norm(x))", delta_ks);


#######################################
# 1b.

log_delta_k = [log(norm(delta_ks[i])) for i = 1:9]
#######################################


#######################################
# 1c.

# plot log of delta_ks 
plot(log_delta_k, markershape = :circle, markercolors = :white, label="||Δ_k||");

# find slope of last two points
alpha = log_delta_k[9] - log_delta_k[8] / log_delta_k[8] - log_delta_k[7]
println("Slope between last two points: α = $alpha")
#######################################



#######################################
# 1d.

# calculate X_k for k = 8, 9, 10
X_8_9_10 = [pol_driver(copy(A), k) for k = 8:10]

# helper method to calculate X^T * X
check_orth = (X) -> X' * X;

# display the resulting matrices
map(x -> display(check_orth(x)), X_8_9_10)



#######################################
# 1e.

W = X_8_9_10[3]
P = inv(W) * copy(A)
println("Eigenvalues of P = $(eigvals(P))");
println("Eigenvalues of A = $(eigvals(copy(A)))");
#######################################
