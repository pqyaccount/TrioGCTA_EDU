
using GREMLModels
using SnpArrays, DataFrames, CSV, LinearAlgebra
using StatsBase, StatsModels
using JLD, HDF5

function vec2ugrm(v::AbstractVector{T}) where {T <: AbstractFloat}
    k = length(v) # unique elements / triangle + diagonal
    p = Int((sqrt(8k + 1) - 1) / 2) # rows/cols in matrix
    G = zeros(T, p, p)
    l = 1
    @inbounds for i in 1:p # cols
        @inbounds for j in 1:i # rows
            G[j, i] = v[l]
            l += 1
        end
    end
    
    G
end

function read_grm(filename::AbstractString) 
    nme = "grm_u"
    g = h5open(filename * ".h5", "r") do file
        read(file, nme)
    end
    vec2ugrm(g)
end

function findGRMs(A, m_inds_use, f_inds_use, c_inds_use) 
    A_m = A[m_inds_use, m_inds_use]
    A_f = A[f_inds_use, f_inds_use]
    A_c = A[c_inds_use, c_inds_use]
    D_fm = A[f_inds_use, m_inds_use] + A[m_inds_use, f_inds_use]
    D_cm = A[c_inds_use, m_inds_use] + A[m_inds_use, c_inds_use]
    D_cf = A[c_inds_use, f_inds_use] + A[f_inds_use, c_inds_use]
    R = Diagonal(ones(size(A_c, 1)))
    [Symmetric(A_m), Symmetric(A_f), Symmetric(A_c), Symmetric(D_fm), Symmetric(D_cm), Symmetric(D_cf), R]
end

function sel_trio(A, m_inds, f_inds, c_inds)
     # Check for relatives, threshold for GRM = 0.05
    mf_inds = vcat(m_inds, f_inds)
    A_mf = A[mf_inds, mf_inds]
    keep_A_mf = kinship_pruning(A_mf, method=:bottom_up, cutoff=0.05)
    m_inds_keep = keep_A_mf[m_inds]
    f_inds_keep = keep_A_mf[f_inds]
    mf_inds_keep = findall((m_inds_keep + f_inds_keep) .== 2)
    m_inds_use = m_inds[mf_inds_keep]
    f_inds_use = f_inds[mf_inds_keep]
    c_inds_use = c_inds[mf_inds_keep]
    (m_inds_use, f_inds_use, c_inds_use)
end

# Choose working folder based on OS
if Sys.iswindows()
    cd("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia")
else
    cd("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia")
end
println(pwd())
println("Direct the folder to Step05_TrioGCTA_Julia")

# Load moba data file with phenotypes, covariates and PCA
if Sys.iswindows()
	mobadat = DataFrame(CSV.File("M:\\p805-qiyuanp\\TrioGCTA\\data\\MoBaPsychGen_v1_1m_gpa_28748obs.moba", missingstring = "NA"))
else
	mobadat = DataFrame(CSV.File("/tsd/p805/home/p805-qiyuanp/TrioGCTA/data/MoBaPsychGen_v1_1m_gpa_28748obs.moba", missingstring = "NA"))
end
println("Load mobadat")

if Sys.iswindows()
	pcdat = DataFrame(CSV.File("N:\\durable\\data\\genetics\\MoBaPsychGen_v1\\MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt", missingstring="NA"))
else
	pcdat = DataFrame(CSV.File("/tsd/p805/data/durable/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt", missingstring="NA"))
end
println("Load PCs")

# Define original size of trio families, focused phenotype variable and threshold of GRM
sz = size(mobadat, 1)
println(sz)

vari = "gpa_normalized"
println(vari)

threshold = "05"
println(threshold)

A = read_grm("grm$(string(sz))")
println("Read grm.h5")

mobadat.order = 1:size(mobadat, 1)
println("Get mobadat order")

mobadat = leftjoin(mobadat, pcdat, on = [:iid => :IID])
println("Leftjoin PCs")

sort!(mobadat, :order)
println("Sort mobadat by order")

# Subset grm blocks, after possibly filtering high relatedness coefficients
A = Symmetric(A, :U)

# Get positions of trios (mother, father, child)
m_inds = findall(mobadat[!,:ROLE] .== "Mother")
f_inds = findall(mobadat[!,:ROLE] .== "Father")
c_inds = findall(mobadat[!,:ROLE] .== "Child")
println("Get positions for trios")

# Check positions for a specific trio (ii=100) and print the GRM matrix for this trio before filtering
ii = 100
trio_inds_i = vcat(m_inds[ii], f_inds[ii], c_inds[ii])
println(A[trio_inds_i, trio_inds_i])

# Remove trios with high relatedness using kinship pruning
m_inds_sel, f_inds_sel, c_inds_sel = sel_trio(A, m_inds, f_inds, c_inds)
println("Remove trios with high relatedness using kinship pruning")

# Extract the new indices for the selected trio after filtering
trio_inds_sel_i = vcat(m_inds_sel[ii], f_inds_sel[ii], c_inds_sel[ii])

# Print the updated GRM matrix for the same trio after filtering
println(A[trio_inds_sel_i, trio_inds_sel_i])

# Fit model
# --------------------------------------

	# Set covariates in right format
	mobadat[!, :sex] .= string.(mobadat[!, :sex])
	mobadat[!, :BATCH] .= string.(mobadat[!, :BATCH])
	
	# Remove families with missing phenotype data
    pos_use = findall(!ismissing, mobadat[c_inds_sel, vari])
	println(length(pos_use))
	println("Trios number after removing relatives and families with missing phenotype data")
	
	# Define children and their phenotype as the focused outcome
    mobadat[c_inds_sel[pos_use],vari]

	# Create mother, father and child dataframe without missing data
    m_inds_sel_nomiss = m_inds_sel[pos_use]
    f_inds_sel_nomiss = f_inds_sel[pos_use]
    c_inds_sel_nomiss = c_inds_sel[pos_use]

	# Get first 20 PCs for mother and father
    m_pc = Matrix(pcdat[m_inds_sel_nomiss, Symbol.(:PC, 1:20)])
    f_pc = Matrix(pcdat[f_inds_sel_nomiss, Symbol.(:PC, 1:20)])

	# Define TrioGCTA Model
    mf = ModelFrame(term(Symbol(vari)) ~ ConstantTerm(1) + term(:BATCH) + term(:sex), mobadat[c_inds_sel_nomiss, :])
    X = ModelMatrix(mf).m
    X = hcat(X, m_pc, f_pc)
    y = response(mf)
  	y = float(y)
	
    rs = findGRMs(A, m_inds_sel_nomiss, f_inds_sel_nomiss, c_inds_sel_nomiss)

# Full model - Mother, Father and Child - TrioGCTA
# --------------------------------------
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
    
	dat = GREMLData(y, X, rs)
    lbs = [0.0, -Inf, -Inf, 0.0, -Inf, 0.0, 0.0]
    vd = var(y)
    ini = vd .* [0.2, 0.0, 0.0, 0.1, 0.0, 0.1, 0.7]
	
    @time m_full = GREMLModel(dat, ini, lbs, false)
	# The optimization must converge completely:
	m_full.opt.ftol_abs = 0
	m_full.opt.ftol_rel = 0
    @time GREMLModels.fit!(m_full)
	println("Fitting Full")
	m_full.δ
	
	filename = vari * "_" * threshold * "_full_$(size(mobadat, 1)).jld"
	println(filename)

	if Sys.iswindows()
	JLD.save("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * filename, "model", m_full)
	else
    JLD.save("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * filename, "model", m_full)
	end

# MO - Mother and Child
# Indirect genetic effects from mother and direct genetic effects from offspring
# -------------------------------------------
	function GREMLModels.transform!(δ::Vector, θ::Vector)
		δ[1] = θ[1]^2
		δ[2] = θ[2]^2 + θ[3]^2
		δ[3] = θ[1] * θ[2]
		δ[4] = θ[4]
		δ
	end
	
	dat = GREMLData(y, X, [rs[1], rs[3], rs[5], rs[7]])
	lbs = [0.0, -Inf, 0.0, 0.0]
	vd = var(y)
	ini = vd .* [0.3, 0.0, 0.0, 0.7]
	
	@time m_mo = GREMLModel(dat, ini, lbs, false)
	# The optimization must converge completely:
	m_mo.opt.ftol_abs = 0
	m_mo.opt.ftol_rel = 0
	@time fit!(m_mo)
	println("Fitting MO")
	m_mo.δ
	
	filename = vari * "_" * threshold * "_mo2_$(size(mobadat, 1)).jld"
	println(filename)

	if Sys.iswindows()
	JLD.save("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * filename, "model", m_mo)
	else
    JLD.save("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * filename, "model", m_mo)
	end

# FO - Father and Child
# Indirect genetic effects from father and direct genetic effects from offspring
# -------------------------------------------
	function GREMLModels.transform!(δ::Vector, θ::Vector)
		δ[1] = θ[1]^2
		δ[2] = θ[2]^2 + θ[3]^2
		δ[3] = θ[1] * θ[2]
		δ[4] = θ[4]
		δ
	end
	
	dat = GREMLData(y, X, [rs[2], rs[3], rs[6], rs[7]])
	lbs = [0.0, -Inf, 0.0, 0.0]
	vd = var(y)
	ini = vd .* [0.3, 0.0, 0.0, 0.7]
	
	@time m_fo = GREMLModel(dat, ini, lbs, false)
	@time fit!(m_fo)
	println("Fitting FO")
	m_fo.δ
	
	filename = vari * "_" * threshold * "_fo3_$(size(mobadat, 1)).jld"
	println(filename)

	if Sys.iswindows()
	JLD.save("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * filename, "model", m_fo)
	else
    JLD.save("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * filename, "model", m_fo)
	end

# MF - Mother and Father
# Indirect genetic effects from mother and father
# -------------------------------------------
	function GREMLModels.transform!(δ::Vector, θ::Vector)
		δ[1] = θ[1]^2
		δ[2] = θ[2]^2 + θ[3]^2
		δ[3] = θ[1] * θ[2]
		δ[4] = θ[4]
		δ
	end
	
	dat = GREMLData(y, X, [rs[1], rs[2], rs[4], rs[7]])
	lbs = [0.0, -Inf, 0.0, 0.0]
	vd = var(y)
	ini = vd .* [0.3, 0.0, 0.0, 0.7]
	
	@time m_mf = GREMLModel(dat, ini, lbs, false)
	@time fit!(m_mf)
	println("Fitting MF")
	m_mf.δ
	
	filename = vari * "_" * threshold * "_mf4_$(size(mobadat, 1)).jld"
	println(filename)

	if Sys.iswindows()
	JLD.save("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * filename, "model", m_mf)
	else
    JLD.save("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * filename, "model", m_mf)
	end

# Indep -  No covariances
# No covariances between direct and indirect effects
# -------------------------------------------
    function GREMLModels.transform!(δ::Vector, θ::Vector)
        δ[1] = θ[1]^2
        δ[2] = θ[2]^2
        δ[3] = θ[3]^2
        δ[4] = θ[4]
        δ
    end
	
    dat = GREMLData(y, X, [rs[1], rs[2], rs[3], rs[7]])
    lbs = [0.0, 0.0, 0.0, 0.0]
	vd = var(y)
    ini = vd .* [0.3, 0.0, 0.0, 0.7]
	
    @time m_indep = GREMLModel(dat, ini, lbs, false)
    @time fit!(m_indep)
	println("Fitting Indep")
	m_indep.δ
	
	filename = vari * "_" * threshold * "_indep5_$(size(mobadat, 1)).jld"
	println(filename)

	if Sys.iswindows()
	JLD.save("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * filename, "model", m_indep)
	else
    JLD.save("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * filename, "model", m_indep)
	end

# Direct - Child
# Direct genetic effects only (offspring)
# -------------------------------------------
    function GREMLModels.transform!(δ::Vector, θ::Vector)
        δ[1] = θ[1]^2
        δ[2] = θ[2]
        δ
    end
	
    dat = GREMLData(y, X, [rs[3], rs[7]])
    lbs = [0.0, 0.0]
	vd = var(y)
    ini = sqrt.(vd .* [0.10, 0.9^2])
	
    @time m_direct = GREMLModel(dat, ini, lbs, false)
    @time fit!(m_direct)
	m_direct.δ
	
    filename = vari * "_" * threshold * "_direct6_$(size(mobadat, 1)).jld"
	println(filename)

	if Sys.iswindows()
	JLD.save("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * filename, "model", m_direct)
	else
    JLD.save("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * filename, "model", m_direct)
	end

	
# RDR
# Average value of parents
# -------------------------------------------

function findGRMs(A, m_inds_use, f_inds_use, c_inds_use) 
    A_m = A[m_inds_use, m_inds_use]
    A_f = A[f_inds_use, f_inds_use]
    A_c = A[c_inds_use, c_inds_use]
    D_fm = A[f_inds_use, m_inds_use] + A[m_inds_use, f_inds_use]
    D_cm = A[c_inds_use, m_inds_use] + A[m_inds_use, c_inds_use]
    D_cf = A[c_inds_use, f_inds_use] + A[f_inds_use, c_inds_use]
    R = Diagonal(ones(size(A_c, 1)))
    [Symmetric(A_m), Symmetric(A_f), Symmetric(A_c), Symmetric(D_fm), Symmetric(D_cm), Symmetric(D_cf), R]
end
	
	rs = findGRMs(A, m_inds_sel_nomiss, f_inds_sel_nomiss, c_inds_sel_nomiss)
	
	A_m = rs[1]
	A_f = rs[2]
	A_c = rs[3]
	D_fm = rs[4]
	D_cm = rs[5]
	D_cf = rs[6]
	R = rs[7]
	
	A_pp = (A_m + A_f + D_fm) ./ 2
	D_op = (D_cm + D_cf) ./ 2
	
	function GREMLModels.transform!(δ::Vector, θ::Vector)
		δ[1] = θ[1]^2
		δ[2] = θ[2]^2 + θ[3]^2
		δ[3] = θ[1] * θ[2]
		δ[4] = θ[4]
		δ
	end

	dat = GREMLData(y, X, [A_pp, rs[3], D_op, rs[7]])
	lbs = [-Inf, -Inf, -Inf, 0.0]
	ini = [0.26, 0.28, 0.16, 0.77]

	@time m_rdr = GREMLModel(dat, ini, lbs, false)
    @time fit!(m_rdr)
	m_rdr.δ

	filename = vari * "_" * threshold * "_rdr7_$(size(mobadat, 1)).jld"
	println(filename)

	if Sys.iswindows()
	JLD.save("M:\\p805-qiyuanp\\TrioGCTA\\Step05_TrioGCTA_Julia\\results\\" * filename, "model", m_rdr)
	else
    JLD.save("/tsd/p805/home/p805-qiyuanp/TrioGCTA/Step05_TrioGCTA_Julia/results/" * filename, "model", m_rdr)
	end

