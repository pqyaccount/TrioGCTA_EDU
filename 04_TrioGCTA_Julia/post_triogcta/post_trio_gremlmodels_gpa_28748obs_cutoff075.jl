using GREMLModels
using JLD
using LinearAlgebra
using StatsModels, StatsBase, Statistics
using DataFrames, CSV

# Define original size of trio families, focused phenotype variable and threshold of GRM to get the filename
sz = "86244"
println(sz)

vari = "gpa_normalized"
println(vari)

threshold = "075"
println(threshold)

# Derive the JLD dataset to compare different StatsModels

# Full model for TrioGCTA
if Sys.iswindows()
    m_full = JLD.load("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_full_" * sz * ".jld")["model"]
else
	m_full = JLD.load("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_full_" * sz * ".jld")["model"]
end
println("Full model")
println(m_full)

# MO
if Sys.iswindows()
    m_mo = JLD.load("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_mo2_" * sz * ".jld")["model"]
else
    m_mo = JLD.load("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_mo2_" * sz * ".jld")["model"]
end
println("MO model")
println(m_mo)

# FO
if Sys.iswindows()
    m_fo = JLD.load("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_fo3_" * sz * ".jld")["model"]
else
    m_fo = JLD.load("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_fo3_" * sz * ".jld")["model"]
end
println("FO model")
println(m_fo)

# MF
if Sys.iswindows()
    m_mf = JLD.load("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_mf4_" * sz * ".jld")["model"]
else
    m_mf = JLD.load("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_mf4_" * sz * ".jld")["model"]
end
println("MF model")
println(m_mf)

# Indep
if Sys.iswindows()
    m_indep = JLD.load("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_indep5_" * sz * ".jld")["model"]
else
    m_indep = JLD.load("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_indep5_" * sz * ".jld")["model"]
end
println("Indep model")
println(m_indep)

# Direct
if Sys.iswindows()
    m_direct = JLD.load("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_direct6_" * sz * ".jld")["model"]
else
    m_direct = JLD.load("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_direct6_" * sz * ".jld")["model"]
end
println("Direct model")
println(m_direct)

# RDR
if Sys.iswindows()
    m_rdr = JLD.load("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_rdr7_" * sz * ".jld")["model"]
else
    m_rdr = JLD.load("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_rdr7_" * sz * ".jld")["model"]
end
println("RDR model")
println(m_rdr)

# Compare fits
# -------------------------------------------
ms = [m_full, m_mo, m_fo, m_mf, m_indep, m_direct, m_rdr]


round.(aic.(ms), digits = 2)
println("AIC:")
println(round.(aic.(ms), digits = 2))

round.(bic.(ms), digits = 2)
println("BIC:")
println(round.(bic.(ms), digits = 2))

round.(deviance.(ms), digits = 2)
println("Deviance:")
println(round.(deviance.(ms), digits = 2))

for i in 2:length(ms)
    println("===== Model ", i, " vs Full Model (ms[1]) =====")
    println(lrtest(ms[i], ms[1]))  # reduced vs full
end

# Full model
# --------------------------------------------
#[Symmetric(A_m), Symmetric(A_f), Symmetric(A_c), Symmetric(D_fm), Symmetric(D_cm), Symmetric(D_cf), R]

# Variance Components
δ_full = m_full.δ
Δ_full = [δ_full[1] δ_full[4] δ_full[5] 0;
          δ_full[4] δ_full[2] δ_full[6] 0;
          δ_full[5] δ_full[6] δ_full[3] 0;
          0 0 0 δ_full[7]]

# Rescale to % variance
v_pos = [1, 2, 3, 5, 6, 7]
var_tot_full = sum(δ_full[v_pos])
sc_full = diagm(fill(inv(sqrt(var_tot_full)), 4))
Δ_std_full = sc_full * Δ_full * sc_full

# M, F, O, FM, OM, OF, R
δ_std_full = [Δ_std_full[1, 1], Δ_std_full[2, 2], Δ_std_full[3, 3], Δ_std_full[2, 1], Δ_std_full[3, 1], Δ_std_full[3, 2], Δ_std_full[4, 4]]
sum(δ_std_full)
print(δ_std_full)

# Standard Error
function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[4]^2
    δ[3] = θ[3]^2 + θ[5]^2 + θ[6]^2
    δ[4] = θ[2] * θ[1]
    δ[5] = θ[3] * θ[1]
    δ[6] = θ[5] * θ[4] + θ[3] * θ[2]
    δ[7] = θ[7]
    δ
end

hessian!(m_full.opt.H, m_full)

function GREMLModels.transform(θ::Vector)
     δ = similar(θ)
     δ[1] = θ[1]^2
     δ[2] = θ[2]^2 + θ[4]^2
     δ[3] = θ[3]^2 + θ[5]^2 + θ[6]^2
     δ[4] = θ[2] * θ[1]
     δ[5] = θ[3] * θ[1]
     δ[6] = θ[5] * θ[4] + θ[3] * θ[2]
     δ[7] = θ[7]
     totvar = sum(δ[[1, 2, 3, 5, 6, 7]])
     δ ./ totvar
end

J_full = jacobian(m_full)
C_full = vcovvc(m_full)
Ct_full = J_full * C_full * J_full'
set_full = sqrt.(diag(Ct_full))
GREMLModels.transform(m_full.θ)

println("Variance components (Full TrioGCTA):")
println(round.(δ_std_full, digits = 3))

println("Standard errors:")
println(round.(set_full, digits = 3))


# Correlations
# s = sqrt.(diag(Δ))
# R = diagm(1 ./ s) * Δ * diagm(1 ./ s)


# MO
# -------------------------------------------
# [rs[1], rs[3], rs[5], rs[7]]

δ_mo = m_mo.δ

# M, O, OM, R
δ_std_mo = δ_mo ./ sum(δ_mo)
sum(δ_std_mo)
print(δ_std_mo)

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ
end

hessian!(m_mo.opt.H, m_mo)

function GREMLModels.transform(θ::Vector)
    δ = similar(θ)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ ./ sum(δ)
end

J_mo = jacobian(m_mo)
C_mo = vcovvc(m_mo)
Ct_mo = J_mo * C_mo * J_mo'
set_mo = sqrt.(diag(Ct_mo))
GREMLModels.transform(m_mo.θ)

println("Variance components (MO Model):")
println(round.(δ_std_mo, digits = 3))

println("Standard errors:")
println(round.(set_mo, digits = 3))

# FO
# -------------------------------------------
# [rs[2], rs[3], rs[6], rs[7]]

δ_fo = m_fo.δ

# F, O, OF, R
δ_std_fo = δ_fo ./ sum(δ_fo)
sum(δ_std_fo)
print(δ_std_fo)

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ
end

hessian!(m_fo.opt.H, m_fo)

function GREMLModels.transform(θ::Vector)
    δ = similar(θ)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ ./ sum(δ)
end

J_fo = jacobian(m_fo)
C_fo = vcovvc(m_fo)
Ct_fo = J_fo * C_fo * J_fo'
set_fo = sqrt.(diag(Ct_fo))
GREMLModels.transform(m_fo.θ)
# works: m_fo.δ/(sum(m_fo.δ))

println("Variance components (FO Model):")
println(round.(δ_std_fo, digits = 3))

println("Standard errors:")
println(round.(set_fo, digits = 3))

# MF
# -------------------------------------------
# [rs[1], rs[2], rs[4], rs[7]]

δ_mf = m_mf.δ
Δ_mf = [δ_mf[1] δ_mf[3] 0;
     δ_mf[3] δ_mf[2] 0;
     0 0 δ_mf[4]]

# Rescale to % variance
v_mf_pos = [1,2,4]
var_tot_mf = sum(δ_mf[v_mf_pos])
sc_mf = diagm(fill(inv(sqrt(var_tot_mf)), 3))
Δ_std_mf = sc_mf * Δ_mf * sc_mf

# M, F, FM, R
δ_std_mf = [Δ_std_mf[1, 1], Δ_std_mf[2, 2], Δ_std_mf[2, 1], Δ_std_mf[3, 3]]
sum(δ_std_mf[[1,2,4]])
print(δ_std_mf[[1,2,4]])

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ
end

hessian!(m_mf.opt.H, m_mf)

function GREMLModels.transform(θ::Vector)
    δ = similar(θ)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ ./ sum(δ[[1, 2, 4]])
end

J_mf = jacobian(m_mf)
C_mf = vcovvc(m_mf)
Ct_mf = J_mf * C_mf * J_mf'
set_mf = sqrt.(diag(Ct_mf))
GREMLModels.transform(m_mf.θ)

println("Variance components (MF Model):")
println(round.(δ_std_mf, digits = 3))

println("Standard errors:")
println(round.(set_mf, digits = 3))


# Indep model
# --------------------------------------------

δ_indep = m_indep.δ

# M, F, O, R
δ_std_indep = δ_indep ./ sum(δ_indep)
sum(δ_std_indep)

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2
    δ[3] = θ[3]^2
    δ[4] = θ[4]
    δ
end

hessian!(m_indep.opt.H, m_indep)

function GREMLModels.transform(θ::Vector)
    δ = similar(θ)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2
    δ[3] = θ[3]^2
    δ[4] = θ[4]
    δ ./ sum(δ)
end

J_indep = jacobian(m_indep)
C_indep = vcovvc(m_indep)
Ct_indep = J_indep * C_indep * J_indep'
set_indep = sqrt.(diag(Ct_indep))
GREMLModels.transform(m_indep.θ)

println("Variance components (Indep Model):")
println(round.(δ_std_indep, digits = 3))

println("Standard errors:")
println(round.(set_indep, digits = 3))

# Direct model
# --------------------------------------------

δ_direct = m_direct.δ

# O, R
δ_std_direct = δ_direct ./ sum(δ_direct)

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]
    δ
end

hessian!(m_direct.opt.H, m_direct)

function GREMLModels.transform(θ::Vector)
    δ = similar(θ)
    δ[1] = θ[1]^2
    δ[2] = θ[2]
    δ ./ sum(δ)
end

J_direct = jacobian(m_direct)
C_direct = vcovvc(m_direct)
Ct_direct = J_direct * C_direct * J_direct'
set_direct = sqrt.(diag(Ct_direct))
GREMLModels.transform(m_direct.θ)

println("Variance components (Direct Model):")
println(round.(δ_std_direct, digits = 3))

println("Standard errors:")
println(round.(set_direct, digits = 3))

# RDR
# -------------------------------------------
# [rs[8], rs[3], rs[9], rs[7]]

δ_rdr = m_rdr.δ

# PP, O, OP, R
δ_std_rdr = δ_rdr ./ sum(δ_rdr)
sum(δ_std_rdr)

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ
end

hessian!(m_rdr.opt.H, m_rdr)

function GREMLModels.transform(θ::Vector)
    δ = similar(θ)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ ./ sum(δ)
end

J_rdr = jacobian(m_rdr)
C_rdr = vcovvc(m_rdr)
Ct_rdr = J_rdr * C_rdr * J_rdr'
set_rdr = sqrt.(diag(Ct_rdr))
GREMLModels.transform(m_rdr.θ)

println("Variance components (RDR Model):")
println(round.(δ_std_rdr, digits = 3))

println("Standard errors:")
println(round.(set_rdr, digits = 3))

# Write results
# ----------------------------------------------

max_len = 7  
fill_nan(vec, len) = vcat(vec, fill(NaN, max_len - length(vec))) 

#  --- Variance and Covariance ---
δ_std_full   = fill_nan(δ_std_full, length(δ_std_full))
δ_std_mo     = fill_nan(δ_std_mo, length(δ_std_mo))
δ_std_fo     = fill_nan(δ_std_fo, length(δ_std_fo))
δ_std_mf     = fill_nan(δ_std_mf, length(δ_std_mf))
δ_std_indep  = fill_nan(δ_std_indep, length(δ_std_indep))
δ_std_direct = fill_nan(δ_std_direct, length(δ_std_direct))
δ_std_rdr    = fill_nan(δ_std_rdr, length(δ_std_rdr))

# --- Standard Errors ---
set_full     = fill_nan(set_full, length(set_full))
set_mo       = fill_nan(set_mo, length(set_mo))
set_fo       = fill_nan(set_fo, length(set_fo))
set_mf       = fill_nan(set_mf, length(set_mf))
set_indep    = fill_nan(set_indep, length(set_indep))
set_direct   = fill_nan(set_direct, length(set_direct))
set_rdr      = fill_nan(set_rdr, length(set_rdr))

par = hcat(δ_std_full, set_full,
           δ_std_mo,   set_mo,
           δ_std_fo,   set_fo,
           δ_std_mf,   set_mf,
           δ_std_indep,set_indep,
           δ_std_direct,set_direct,
           δ_std_rdr,  set_rdr)

datres = DataFrame(par, ["full","full_se","mo","mo_se","fo","fo_se","mf","mf_se","indep","indep_se","direct","direct_se","rdr","rdr_se"])
datres.model = ["m", "f", "o", "mf", "mo", "fo", "e"]

if Sys.iswindows()
    CSV.write("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * vari * "_" * threshold * "_" * sz * "_results.csv", round.(datres[:,1:14], digits = 3))
else
    CSV.write("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * vari * "_" * threshold * "_" * sz * "_results.csv", round.(datres[:,1:14], digits = 3))
end
