####################################
# Initialise simulation parameters #
####################################

# Import libraries
using DifferentialEquations, LinearAlgebra, Symbolics, StaticArrays, Random, JLD2, LsqFit

# Check whether input arguments have been passed to code
inputs = ( if isempty(ARGS) ; ["50.0", # LBuff
                               "0.2",  # beta
                               "1.0",  # delM
                               "1.0",  # delE
                               "0.5",  # thresh
                               "0",    # DFeedback
                               "1",    # numVary
                               "1",    # wVary
                               "0",    # pointSource
                               "0",    # dynamicScalingSim
                               "0",    # dynamicRobustnessSim
                               "0",    # chiSqRatCalc
                               "0",    # originalERModel
                               "0",    # readInData
                               "0",    # saveDists
                               "1"] ;  # ID
                            else ARGS ; end )

# Read in parameters from input arguments
const LBuff = parse(Float64, inputs[1])
const beta = parse(Float64, inputs[2])
const delM = parse(Float64, inputs[3])
const delE = parse(Float64, inputs[4])
const thresh = parse(Float64, inputs[5])
const DFeedback = parse(Int, inputs[6]) # 1 == num will be varied, else other parameters will be varied instead
const numVary = parse(Int, inputs[7]) # 1 == num will be varied, else other parameters will be varied instead
const wVary = parse(Int, inputs[8]) # 1 == source wwidth will increase proportionally with L, else will be kept constant
const pointSource = parse(Int, inputs[9]) # 1 == source width will be equal to a single grid point, else will be defined using beta (supercedes wVaryb definition)
const dynamicScalingSim = parse(Int, inputs[10]) # 1 == system will be simulated for multiple L values only (supercedes all other types of simulation) and dynamic scaling will be calculated
const dynamicRobustnessSim = parse(Int, inputs[11]) # 1 == system will be simulated for multiple num values only (supercedes all other types of simulation)
const chiSqRatCalc = parse(Int, inputs[12]) # 1 == fit morphogen profiles with exponential and power-law and compare chi-squared values
const originalERModel = parse(Int, inputs[13]) # 1 == solve the original form of the expansion-repression model without expander self-expansion
const readInData = parse(Int, inputs[14]) # 1 == seed data will be read in from file if numVary == 1
const saveDists = parse(UInt, inputs[15]) # 1 == steady-state morphogen and expander profiles will be saved for each L and num
const ID = parse(UInt, inputs[16])

# Since in the case of a point source/input flux the source width is undefined, only utilise the solver for wVary == 1 and not wVary == 0
if((pointSource == 1) && (wVary == 0))
    error("ERROR: Cannot use pointSource == 1 and wVary == 0")
end

# Due to changes in source indexing, can only reliably automate dynamic scaling calculations for scaling source width simulations
if((dynamicScalingSim == 1) && (wVary == 0))
    error("ERROR: Cannot use dynamicScalingSim == 1 and wVary == 0")
end

# Only interested in simulating the original ER model using a point source as in the original paper, and want to minimise number of solves so use numVary == 1
if(((originalERModel == 1) && (pointSource == 0)) || ((originalERModel == 1) && (numVary == 0)))
    error("ERROR: Cannot use originalERModel == 1 and either pointSource == 0 or numVary == 0")
end

# Define how many simulations (independent sets of input parameters) will be run
const nSims::Int = 10

# Initialise grid size and morphogen source region width
const N::Int = 201
const N1::Int = N - 1

if(pointSource == 1)
    const sIndex::Int = 1 # Set different values of sWidth and sIndex1Arr to calculate suitable values of bIndexes and posL1 (calculate scaling over the whole domain)
    const sIndex1::Int = sIndex + 1 # This definition ensures that only position 1 in the solver receives the input flux
elseif(wVary == 1)
    const sIndex::Int = Int(round(beta * N1))
    const sIndex1::Int = sIndex + 1
else
    const sIndex::Int = Int(round(((beta * 50.) / LBuff) * N1)) # Assume source width at L = 50 is to be kept constant
    const sIndex1::Int = sIndex + 1

    const sWConstIndex::Int = Int(round(sIndex / 2)) # Only valid for doubling of tissue size
    const sWConstIndex1::Int = sWConstIndex + 1
    const NWConst::Int = N - sWConstIndex
    const NWConst1::Int = NWConst - 1
end

# Initialise tissue sizes for simulations
const LBuff_2::Float64 = LBuff / 2.
if(dynamicScalingSim == 1)

    # Quantify dynamic scaling across a small range of tissue lengths (see Fig.(S5,S14a))
    LVec = SVector{6, Float64}(LBuff, 2. * LBuff, 4. * LBuff, 6. * LBuff, 8. * LBuff, 10. * LBuff)
    const nLVec::Int = length(LVec)
    dxVec::SVector{6, Float64} = LVec / N1
    sWidth::SVector{6, Float64} = (beta .* LVec[1], beta .* LVec[2], beta .* LVec[3], beta .* LVec[4], beta .* LVec[5], beta .* LVec[6])
    NArr::SVector{6, Int} = [N, N, N, N, N, N]
    sIndex1Arr::SVector{6, Int} = [sIndex1, sIndex1, sIndex1, sIndex1, sIndex1, sIndex1]

    # Quantify dynamic scaling across a large range of tissue lengths (see Fig.(S14b-d))
    #= const nLVec::Int = 101
    LVec::Vector{Float64} = [LBuff ^ (((1.765776 - 1.) * ((i - 1.) / (nLVec - 1.))) + 1.) for i in 1:nLVec]
    println(LVec)
    dxVec::Vector{Float64} = LVec / N1
    sWidth::Vector{Float64} = pointSource == 1 ? zeros(nLVec) : beta .* LVec
    NArr::Vector{Int} = [N for i in 1:nLVec]
    sIndex1Arr::Vector{Int} = pointSource == 1 ? ones(nLVec) : sIndex1 * ones(nLVec) =#

elseif(dynamicRobustnessSim == 1)

    # Quantify robustness across a small range of parameter values (see Fig.(S15a))
    LVec = SVector{15, Float64}(LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff)
    const nLVec::Int = length(LVec)
    dxVec::SVector{15, Float64} = LVec / N1
    sWidth::SVector{15, Float64} = pointSource == 1 ? (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.) : (beta .* LVec[1], beta .* LVec[2], beta .* LVec[3], beta .* LVec[4], beta .* LVec[5], beta .* LVec[6], beta .* LVec[7], beta .* LVec[8], beta .* LVec[9], beta .* LVec[10], beta .* LVec[11], beta .* LVec[12], beta .* LVec[13], beta .* LVec[14], beta .* LVec[15])
    NArr::SVector{15, Int} = [N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
    sIndex1Arr::SVector{15, Int} = pointSource == 1 ? ones(15) : sIndex1 * ones(15)

    # Quantify robustness across a large range of parameter values (see Fig.(S15b-d))
    #= const nLVec::Int = 81
    LVec::Vector{Float64} = LBuff .* ones(Float64, nLVec)
    println(LVec)
    dxVec::Vector{Float64} = LVec / N1
    sWidth::Vector{Float64} = pointSource == 1 ? zeros(nLVec) : beta .* LVec
    NArr::Vector{Int} = [N for i in 1:nLVec]
    sIndex1Arr::Vector{Int} = pointSource == 1 ? ones(nLVec) : sIndex1 * ones(nLVec) =#

elseif(numVary == 1)
    if(wVary == 1)
		LVec = SVector{4, Float64}(LBuff, LBuff, LBuff, 2. * LBuff)
		const nLVec::Int = length(LVec)
		dxVec::SVector{4, Float64} = LVec / N1
        sWidth::SVector{4, Float64} = pointSource == 1 ? (0., 0., 0., 0.) : (beta .* LVec[1], beta .* LVec[2], beta .* LVec[3], beta .* LVec[4])
        NArr::SVector{4, Int} = [N, N, N, N]
        sIndex1Arr::SVector{4, Int} = pointSource == 1 ? [1., 1., 1., 1.] : [sIndex1, sIndex1, sIndex1, sIndex1]
	else
		LVec = SVector{4, Float64}(LBuff, LBuff, LBuff, (2 * LBuff) - (beta * 50.))
		const nLVec::Int = length(LVec)
		dxVec::SVector{4, Float64} = [LVec[i] == LBuff ? LVec[i] / N1 : LVec[i] / NWConst1 for i in 1:nLVec]
        sWidth::SVector{4, Float64} = (beta * 50., beta * 50., beta * 50., beta * 50.)
        NArr::SVector{4, Int} = [N, N, N, NWConst]
        sIndex1Arr::SVector{4, Int} = [sIndex1, sIndex1, sIndex1, sWConstIndex1]
	end
else
	LVec = SVector{11, Float64}(LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff, LBuff)
	const nLVec::Int = length(LVec)
	dxVec::SVector{nLVec, Float64} = LVec / N1
    sWidth::SVector{nLVec, Float64} = (beta .* LVec[1], beta .* LVec[2], beta .* LVec[3], beta .* LVec[4], beta .* LVec[5], beta .* LVec[6], beta .* LVec[7], beta .* LVec[8], beta .* LVec[9], beta .* LVec[10], beta .* LVec[11])
    NArr::SVector{nLVec, Int} = [N, N, N, N, N, N, N, N, N, N, N]
    sIndex1Arr::SVector{nLVec, Int} = [sIndex1, sIndex1, sIndex1, sIndex1, sIndex1, sIndex1, sIndex1, sIndex1, sIndex1, sIndex1, sIndex1]
end

xL1::Vector{Float64} = [(i - 1) * dxVec[1] for i in 1:N]
xL1_sIndex::Vector{Float64} = xL1[sIndex1Arr[1]:N]
const xL1End::Float64 = xL1[end]
const xL1SWidth::Float64 = xL1[sIndex1Arr[1]]
const xL1End_SWidth::Float64 = xL1End - xL1SWidth
const dxL1::Float64 = dxVec[1]

if(wVary == 1)
    xL2::Vector{Float64} = [(i - 1) * dxVec[end] for i in 1:N]
    const xL2End::Float64 = xL2[end]
    const xL2SWidth::Float64 = xL2[sIndex1Arr[end]]
else
    xL2::Vector{Float64} = [(i - 1) * dxVec[end] for i in 1:NWConst]
    const xL2End::Float64 = xL2[end]
    const xL2SWidth::Float64 = xL2[sWConstIndex1]
end
const xL2End_SWidth::Float64 = xL2End - xL2SWidth
const dxL2::Float64 = dxVec[4]

const dr::Float64 = 1. / (N1 - sIndex)
const dr2::Float64 = 2. * dr
const LCell::Float64 = 1.

# Initialise indexes and positions corresponding to gene expression boundary positions
const nBIndexes::Int = 19
bIndexes = SVector{nBIndexes, Int}(sIndex1Arr[1] .+ [i * Int((N - sIndex1Arr[1]) / (nBIndexes + 1)) for i in 1:nBIndexes])
posL1 = SVector{nBIndexes, Float64}((xL1[bIndexes] .- xL1SWidth) ./ xL1End_SWidth)
# println(bIndexes)
# println(xL1[bIndexes] ./ LVec[1])
# println(posL1)

# Initialise constants for expander production
const h::Int = 4
const threshPow::Float64 = thresh ^ h

# Initialise factors that define changes in morphogen production rate
# Store information about whether the system has a point source or not in delNumVec (need production term / dx compared to case without input flux)
if(dynamicScalingSim == 1)

    delNumVec::Vector{Float64} = pointSource == 1 ? 1 ./ dxVec : ones(nLVec)
	delThreshVec::Vector{Float64} = ones(nLVec)

elseif(dynamicRobustnessSim == 1)

    # Quantify robustness across a small range of parameter values (see Fig.(S15a))
    delNumVec = pointSource == 1 ? SVector{nLVec, Float64}(1. / dxVec[1], 1. / dxVec[2], 1. / dxVec[3], 1. / dxVec[4], 1. / dxVec[5], 1. / dxVec[6], 1. / dxVec[7], 1. / dxVec[8], 1. / dxVec[9], 1. / dxVec[10], 1. / dxVec[11], 1. / dxVec[12], 1. / dxVec[13], 1. / dxVec[14], 1. / dxVec[15]) : SVector{nLVec, Float64}(0.2, 0.225, 0.25, 0.4, 0.425, 0.405, 1., 1.025, 1.05, 2.5, 2.525, 2.55, 5., 5.025, 5.05)
	delThreshVec = SVector{nLVec, Float64}(1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)

    # Quantify robustness across a large range of parameter values (see Fig.(S15b-d))
    #= dNumVec::Float64 = 2. / (nLVec - 1.)
    delNumVec::Vector{Float64} = pointSource == 1 ? 1 ./ dxVec : [10 ^ (((dNumVec * (i - 1))) - 1) for i in 1:nLVec]
    println(delNumVec)
	delThreshVec::Vector{Float64} = ones(nLVec) =#

elseif(numVary == 1)

	delNumVec = pointSource == 1 ? SVector{nLVec, Float64}(1. / dxVec[1], 0.666 / dxVec[2], 1.5 / dxVec[3], 1. / dxVec[4]) : SVector{nLVec, Float64}(1., 0.666, 1.5, 1.)
	delThreshVec = SVector{nLVec, Float64}(1., 1., 1., 1.)

else

	delDmVec = SVector{nLVec, Float64}(1., 0.666, 1.5, 1., 1., 1., 1., 1., 1., 1., 1.)
	delDeVec = SVector{nLVec, Float64}(1., 1., 1., 0.666, 1.5, 1., 1., 1., 1., 1., 1.)
    delNumVec = pointSource == 1 ? SVector{nLVec, Float64}(1. / dxVec[1], 1. / dxVec[2], 1. / dxVec[3], 1. / dxVec[4], 1. / dxVec[5], 1. / dxVec[6], 1. / dxVec[7], 1. / dxVec[8], 1. / dxVec[9], 1. / dxVec[10], 1. / dxVec[11]) : SVector{nLVec, Float64}(1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)
	delKVec = SVector{nLVec, Float64}(1., 1., 1., 1., 1., 0.666, 1.5, 1., 1., 1., 1.)
	delMuVec = SVector{nLVec, Float64}(1., 1., 1., 1., 1., 1., 1., 0.666, 1.5, 1., 1.)
	delThreshVec = SVector{nLVec, Float64}(1., 1., 1., 1., 1., 1., 1., 1., 1., 0.666, 1.5)
    currThresh::Float64 = thresh
    currThreshPow::Float64 = thresh ^ h

end

# Inititalise tolerances for steady-state detection
const absTol::Float64 = 1e-8
const relTol::Float64 = 1e-6

# Initialise ranges and minimum values for diffusivity, degradation and production parameters
range = SVector{6, Float64}(4., 4., 4., 4. , 4., 4.) # Use same minimum and maximum values as for evolution simulations
minVal = SVector{6, Float64}(-3., -3., -3., -3., -3., -3.)

# Initialise flag and time variables required to decide when system has reached steady state
flagSS::Int = 0
tAnySS::Float64 = 0. # Stores the time when derivatives fall below threshold closest to final simulation time
tPrevSS::Float64 = 0. # Stores the time when derivatives fall below threshold or system fails check which must be repeated closest to final simulation time
mTotPrevSS::Float64 = 0.
mTotDevMax::Float64 = 0.
const nSSCheck::Int = 1e2
const ssCheckTime = 5e2

# Initialise flag storing whether system reached a NESS and variables used to calculate when this occurs
nessFlag::Int = 0
tPeakMPrev::Float64 = 0.
uPeakMPrev::Float64 = 0.
typePeakMPrev::Int = 0
uPeakM1 = Vector{Float64}()
uPeakM2 = Vector{Float64}()

# Import functions defined in other files
include("ER_Sweep_Functions.jl")

# Initialise the folder name for saving data
folder = @__DIR__



####################################
# Initialise solvers               #
####################################

# Define initial conditions and time-scales for solver
u0::Matrix{Float64} = ones(N, 2) * 1e-10
tSpan = (0., Inf)
const dtMin::Float64 = 1e-6

# Initialise parameters for calculating the Jacobian
du0 = copy(u0)
p = ones(Float64, 7)

# Initialise Jacobian and problem for simulations where the production rate is varied to probe system robustness
if(numVary == 1)

    # Check which type of model is being used
    if(originalERModel == 1)

        if(DFeedback == 0)
            jacK = Symbolics.jacobian_sparsity((du, u) -> fKFeedbackOriginal!(du, u, p, 0.), du0, u0)
            fKSparse = ODEFunction(fKFeedbackOriginal! ; jac_prototype = float.(jacK))
        else
            jacD = Symbolics.jacobian_sparsity((du, u) -> fDFeedbackOriginal!(du, u, p, 0.), du0, u0)
            fDSparse = ODEFunction(fDFeedbackOriginal! ; jac_prototype = float.(jacD))
        end

    else

        if(DFeedback == 0)
            jacK = Symbolics.jacobian_sparsity((du, u) -> fKFeedback!(du, u, p, 0.), du0, u0)
            fKSparse = ODEFunction(fKFeedback! ; jac_prototype = float.(jacK))
        else
            jacD = Symbolics.jacobian_sparsity((du, u) -> fDFeedback!(du, u, p, 0.), du0, u0)
            fDSparse = ODEFunction(fDFeedback! ; jac_prototype = float.(jacD))
        end

    end

else # Initialise Jacobian and problem for system with feedback in other system parameters

    if(DFeedback == 0)
        jacK = Symbolics.jacobian_sparsity((du, u) -> fKFeedback_threshVary!(du, u, p, 0.), du0, u0)
        fKSparse = ODEFunction(fKFeedback_threshVary! ; jac_prototype = float.(jacK))
    else
        jacD = Symbolics.jacobian_sparsity((du, u) -> fDFeedback_threshVary!(du, u, p, 0.), du0, u0)
        fDSparse = ODEFunction(fDFeedback_threshVary! ; jac_prototype = float.(jacD))
    end

end

# Convert Jacobian into problem for solver
if(DFeedback == 0)
    probK = ODEProblem(fKSparse, u0, tSpan, p)
else
    probD = ODEProblem(fDSparse, u0, tSpan, p)
end

# Initialise Jacobian and problem for different types of sources, need both this and normal definition for simulations
if(wVary == 0)

    # Initialise Jacobian and problem for system with same feedback but a constant source width
    u0WConst::Matrix{Float64} = ones(NWConst, 2) * 1e-10
    du0WConst = copy(u0WConst)

    if(DFeedback == 0)
        jacKWConst = Symbolics.jacobian_sparsity((du, u) -> fKFeedbackWConst!(du, u, p, 0.), du0WConst, u0WConst)
        fKWConstSparse = ODEFunction(fKFeedbackWConst! ; jac_prototype = float.(jacKWConst))
        probKWConst = ODEProblem(fKWConstSparse, u0WConst, tSpan, p)
    else
        jacDWConst = Symbolics.jacobian_sparsity((du, u) -> fDFeedbackWConst!(du, u, p, 0.), du0WConst, u0WConst)
        fDWConstSparse = ODEFunction(fDFeedbackWConst! ; jac_prototype = float.(jacDWConst))
        probDWConst = ODEProblem(fDWConstSparse, u0WConst, tSpan, p)
    end

end

# Initialise combined callback functions for detecting steady state, maintaining positive concentrations, oscillations and linear expander growth
cbTerm = DiscreteCallback(newTerminateSteadyStateCondition, newTerminateSteadyStateAffect!, save_positions = (false, false))
cbPosDomain = PositiveDomain()
cbOscCheckM = ContinuousCallback(oscCheckMCondition, oscCheckMAffect!, save_positions = (false, false))
cbSet = CallbackSet(cbTerm, cbPosDomain, cbOscCheckM)



####################################
# Read in scaling seed data        #
####################################

# Initialise arrays to store simulation ID labels and seeds
currSimInds::Vector{Int} = []
seeds::Vector{UInt} = []

# Check if data should be read in or not
if(readInData == 1)

    # Set ID for reading in scaling seed data
    scalingFolder = @__DIR__
    scalingFiles = filter(x -> startswith(x, "GoodShape_Scaling_Seeds"), readdir(scalingFolder))
    # println(unfinishedFiles)

    # Read in scaling data in a new environment so that input data is deleted in the clean-up step when exiting
    # Data is stored as: (1) seeds, (2) morphogen scaling profiles, (3) expander scaling profiles, (4) morphogen robustness profiles, (5) expander robustness profiles
    let scalingData::Vector{Any} = reduce(vcat, [load_object(scalingFolder * "/" * scalingFiles[i]) for i in eachindex(scalingFiles)])

        # Count the number of input data elements
        nSeeds::Int = length(scalingData[1])
        maxSeedID::Int = Int(ceil(nSeeds / nSims))

        # Check whether number of simulations needs to be modified
        if(ID == maxSeedID)
            global currSimInds = ((ID - 1) * nSims) + 1:nSeeds
            global nSims = length(currSimInds)
        elseif(ID > maxSeedID)
            global currSimInds = []
            global nSims = length(currSimInds)
        else
            global currSimInds = ((ID - 1) * nSims) + 1:ID * nSims
            # Do not need to change number of simulations
        end

        # Extract data for any simulations to be run
        if(nSims > 0)
            println("$(nSims) simulations: $(currSimInds[1]):$(currSimInds[end]) / $(nSeeds)")

            # Extract seeds for simulations
            global seeds = scalingData[1][currSimInds]

            # Extract original scaling and robustness profiles to store with outputted data
            global SmOriginal::Vector{Vector{Float64}} = scalingData[2][currSimInds]
            global SmNormOriginal::Vector{Vector{Float64}} = scalingData[3][currSimInds]
            global SeOriginal::Vector{Vector{Float64}} = scalingData[4][currSimInds]
            global RmOriginal::Vector{Vector{Float64}} = [scalingData[5][i][1] for i in currSimInds]
            global ReOriginal::Vector{Vector{Float64}} = [scalingData[6][i][1] for i in currSimInds]

        end

    end

else

    # Set seeds to be simulated in current simulation
    global currSimInds = ((ID - 1) * nSims) + 1:ID * nSims
    # println("$(currSimInds[1]):$(currSimInds[end])")
    global seeds = [UInt(round(time())) + ID + (1e6i) for i in 1:nSims]
    # println(seeds)

    # Initialise empty original scaling and robustness profiles
    global SmOriginal::Vector{Vector{Float64}} = [[] for i in 1:nSims]
    global SmNormOriginal::Vector{Vector{Float64}} = [[] for i in 1:nSims]
    global SeOriginal::Vector{Vector{Float64}} = [[] for i in 1:nSims]
    global RmOriginal::Vector{Vector{Float64}} = [[] for i in 1:nSims]
    global ReOriginal::Vector{Vector{Float64}} = [[] for i in 1:nSims]

end



####################################
# Loop over simulations            #
####################################

# Loop over simulations
@inbounds for i in 1:nSims
    # if(mod(i, 50) == 0)
        println("\nParameter set $(i) / $(nSims):")
    # end

    # Seed random numbers for remaining input parameters
    nSeed::UInt = seeds[i]

    # nSeed::UInt = [0x000000006a5e07ae, 0x000000007138d4f0, 0x00000000735e29ad, 0x00000000695aa551][i] # Best scaling systems in Fig.(3c)

    # Use random number seed to generate random input parameters
    Random.seed!(nSeed)
    r::Vector{Float64} = rand(6)
    @. r = 10 ^ ((range * r) + minVal)

    # Store current set of input parameters that can be used to calculate scaling/robustness
    # In order: dm, de, num, nue, nueThreshPow, k1, mu
    pNew::Vector{Float64} = [r[1], r[2], r[3], r[4], r[4] * threshPow, r[5], r[6]]
    # println(pNew)

    # Save empty data set
    save_object(folder * "/$(currSimInds[i])_SweepData.jld2", SimData(nSeed, 
                                                                      pNew,
                                                                      Vector{Int}(), 
                                                                      Vector{Float64}(), 
                                                                      # Vector{Int}(), # Can optionally store the number of peaks following perturbation
                                                                      # Vector{Vector{Float64}}(), # Can optionally store final steady-state morphogen distributions
                                                                      # Vector{Vector{Float64}}(), # Can optionally store final steady-state expander distributions
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Vector{Float64}}(), 
                                                                      Vector{Vector{Float64}}(),
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Vector{Float64}}(),
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(), 
                                                                      Vector{Float64}(),
                                                                      Vector{Float64}(),
                                                                      Vector{Float64}(),
                                                                      -1e6))

    # Extract original scaling and robustness profiles to output with new simulation data
    global SmOriginalCurr = SmOriginal[i]
    global SmNormOriginalCurr = SmNormOriginal[i]
    global SeOriginalCurr = SeOriginal[i]
    global RmOriginalCurr = RmOriginal[i]
    global ReOriginalCurr = ReOriginal[i]

    # Solve system numerically and save results
    if(wVary == 1)

        if(DFeedback == 0)
            data = solveSystem(nSeed, [pNew[1:5]..., pNew[6] * delM, pNew[7] * delE], probK, probK, LVec, dxVec)
        else
            data = solveSystem(nSeed, [pNew[1:5]..., pNew[6] * delM, pNew[7] * delE], probD, probD, LVec, dxVec)
        end

    else

        if(DFeedback == 0)
            data = solveSystem(nSeed, [pNew[1] / delM, pNew[2] / delE, pNew[3:7]...], probK, probKWConst, LVec, dxVec)
        else
            data = solveSystem(nSeed, [pNew[1] / delM, pNew[2] / delE, pNew[3:7]...], probD, probDWConst, LVec, dxVec)
        end

    end
    # println(currSimInds[i])
    save_object(folder * "/$(currSimInds[i])_SweepData.jld2", data)

end
