#####################################
# List of functions for SweepXX.jl  #
#####################################

# Define struct for storing simulation data
struct SimData # Does not need to be mutable as only ever construct the entire struct at once

	seed::UInt # Store the seeds for each set of random number generation
	paras::Vector{Float64} # Store input parameters (p)
    ssReached::Vector{Int} # Store whether system reached steady-state for each system size and production rate
    tSS::Vector{Float64} # Store the time taken for each system to reach steady-state
    # nPeaks::Vector{Int} # Store the number of peaks following perturbation
	# ssMDist::Vector{Vector{Float64}} # Store final steady-state morphogen distributions for each system size and production rate
	# ssEDist::Vector{Vector{Float64}} # Repeat for final steady-state expander distributions
	Sm::Vector{Float64} # Store relative changes in gene expression boundary positions for the morphogen when tissue length is changed
	SmNorm::Vector{Float64} # Store relative changes in gene expression boundary positions for the normalised morphogen when tissue length is changed
	SmDyn::Vector{Vector{Float64}} # Store relative changes in gene expression boundary positions for the morphogen when tissue length is changed many times in a single simulation (dynamicScalingSim == 1)
	Rm::Vector{Vector{Float64}} # Store relative changes in gene expression boundary positions for the morphogen at each tissue length
	Pm::Vector{Float64} # Store the local precision of the morphogen at each gene expression boundary position for the base tissue length and production rate
	mRange::Vector{Float64} # Store fractional range for each morphogen distribution
    mHDL::Vector{Float64} # Store half-decay lengths for each morphogen distribution
    mAmp0::Vector{Float64} # Store amplitude for each morphogen distribution at x = 0
    mAmpSWidth::Vector{Float64} # Store amplitude for each morphogen distribution at x = w
    mAmpL::Vector{Float64} # Store amplitude for each morphogen distribution at x = L
	mMaxRatio::Vector{Float64} # Store ratio between the maximum morphogen concentration and the concentration at x = 0 to check for spikes
    Se::Vector{Float64} # Store relative changes in gene expression boundary positions for the expander when tissue length is changed
    Re::Vector{Vector{Float64}} # Store relative changes in gene expression boundary positions for the expander at each tissue length
    eRange::Vector{Float64} # Store fractional range for each expander distribution
    eHDL::Vector{Float64} # Store half-decay lengths for each expander distribution
    eAmp::Vector{Float64} # Store amplitude for each morphogen expander within the target tissue
	eWidth::Vector{Float64} # Store width of expander source region for each system size and production rate
	eHighFrac::Vector{Int} # Store number of grid points where expander concentration is higher than Hill constants
	muAv::Vector{Float64} # Store the average expander degradation rate for each system size and production rate
    chiSqRat::Float64 # Store the ratio of the chi-squared values for exponential and power-law fits to the morphogen gradient at base-line tissue length and production rate

end

# Define exponential and quadratic morphogen fit functions
@views mFitExp(x, p) = (p[1] .* exp.((x .- p[2]) ./ p[3])) .+ (p[4] .* exp.(-(x .- p[2]) ./ p[3]));
@views mFitPow(x, p) = (p[1] .* ((x .- p[2]) .^ p[3])) .+ (p[4] .* ((x .- p[2]) .^ (1 - p[3])));

# Define function that calculates the fitness of a system with input parameters defined in the vector pNew
function solveSystem(nSeed, pNew, probL1, probL2, LVec, dxVec)

	# Store solutions for current set of parameters to use to calculate scaling/robustness
    ssReached = zeros(Int, nLVec)
	tSS = -ones(Float64, nLVec)
    nPeaks = zeros(Int, nLVec)
	solMArr = Vector{Vector{Float64}}(undef, nLVec)
	solEArr = Vector{Vector{Float64}}(undef, nLVec)

	@inbounds for n in 1:nLVec
		println("\n   Loop $(n)")

		# Initialise grid for current loop
		dx2 = dxVec[n] ^ 2

		# Zero flags and times used to decide whether system has reached steady state
		global flagSS = 0
		global tAnySS = 0.
		global tPrevSS = 0.
		global mTotPrevSS = 0. # These two definitions are technically not required here
		global mTotDevMax = 0.

		# Zero flags and arrays storing whether system reached a NESS or a linear expander growth state
		global nessFlag = 0
		global tPeakMPrev = 0.1
		global uPeakMPrev = 0.
		global typePeakMPrev = 0
		empty!(uPeakM1)
		empty!(uPeakM2)



		####################################
		# Solve problem                    #
		####################################

		# Solve ODE to output results
		try

			# Define initial condition of solver using system solution before production rate perturbation
			if(((wVary == 1) || (numVary == 0)) && (n > 1) && (ssReached[1] == 1)) @views u0New = hcat(solMArr[1], solEArr[1])
			elseif((wVary == 0) && ((n == 2) || (n == 3)) && (ssReached[1] == 1)) @views u0New = hcat(solMArr[1], solEArr[1])
			elseif((wVary == 0) && (n == 4) && (ssReached[1] == 1)) @views u0New = hcat(solMArr[1][1:NWConst], solEArr[1][1:NWConst])
			else u0New = u0[1:NArr[n], :]
			end

			# Use different parameters depending on what is being varyied
			if(numVary == 1)

				# Can use different solvers depending on system size
				if(n <= 3)

					# Recall that information about whether the system has a point source or not is included in delNumVec
					@views sol = solve(remake(probL1, u0 = u0New, p = (pNew[1] / dx2, pNew[2] / dx2, pNew[3] * delNumVec[n], pNew[4:5]..., pNew[6], pNew[7])),
										Rodas4P(),
										maxiters = 1e6, 
										save_everystep = true, 
										verbose = false, 
										dtmin = dtMin, 
										force_dtmin = false,
										callback = cbSet, # isoutofdomain = (y, p, t) -> any(x -> x < 0, y)
										abstol = absTol, 
										reltol = relTol)

				else

					@views sol = solve(remake(probL2, u0 = u0New, p = (pNew[1] / dx2, pNew[2] / dx2, pNew[3] * delNumVec[n], pNew[4:5]..., pNew[6], pNew[7])),
										Rodas4P(),
										maxiters = 1e6, 
										save_everystep = true, 
										verbose = false, 
										dtmin = dtMin, 
										force_dtmin = false,
										callback = cbSet, # isoutofdomain = (y, p, t) -> any(x -> x < 0, y)
										abstol = absTol, 
										reltol = relTol)
				
				end

			else

				# Set current value of expander production threshold (thresh)
				global currThresh = thresh * delThreshVec[n]
				global currThreshPow = currThresh ^ h
				# println(currThreshPow)
				# println([(pNew[1] * delDmVec[n]) / dx2, (pNew[2] * delDeVec[n]) / dx2, pNew[3:4]..., pNew[4] * currThreshPow, pNew[6] * delKVec[n] * delM, pNew[7] * delMuVec[n] * delE])

				# Vary dynamical parameters other than the morphogen production rate (delNumVec only contains information about flux vs. source in this case)
				@views sol = solve(remake(probL1, u0 = u0New, p = ((pNew[1] * delDmVec[n]) / dx2, (pNew[2] * delDeVec[n]) / dx2, pNew[3] * delNumVec[n], pNew[4], pNew[4] * currThreshPow, pNew[6] * delKVec[n], pNew[7] * delMuVec[n])),
									Rodas4P(),
									maxiters = 1e6, 
									save_everystep = true, 
									verbose = false, 
									dtmin = dtMin, 
									force_dtmin = false,
									callback = cbSet, # isoutofdomain = (y, p, t) -> any(x -> x < 0, y)
									abstol = absTol, 
									reltol = relTol)

			end

			# Check if system was terminated before solver reached MaxIters
			if(sol.retcode != :MaxIters)

				if(nessFlag == 0)

                    # Indicate that steady-state was reached in finite time
                    ssReached[n] = 1

					# Extract the indexes corresponding to near-steady-state
					tSSInds::Vector{Int} = findall(x -> x >= tAnySS, sol.t)
					nTSSInds::Int = length(tSSInds)

					# Calculate the final total morphogen and expander concentrations
					@views mTots::Vector{Float64} = [trapezium(sol.u[i][:, 1]) for i in tSSInds]
					@views eTots::Vector{Float64} = [trapezium(sol.u[i][:, 2]) for i in tSSInds]

					# Interpolate to find an approximation for the time taken to reach steady-state
					finSSInd::Int = 1
					@inbounds for i in nTSSInds - 1:-1:1
						@views if((abs(mTots[i] - mTots[end]) / mTots[end] > 0.005) || (abs(eTots[i] - eTots[end]) / eTots[end] > 0.005))
							finSSInd = i + 1
							break
						end
					end
					# println("$(finSSInd) / $(nTSSInds)")

                    # Store solution
					# println("Reached steady-state at: $(sol.t[tSSInds[finSSInd]]) (not $(tAnySS), $(tPrevSS) or $(sol.t[end]))\n")
                    tSS[n] = sol.t[tSSInds[finSSInd]]
                    nPeaks[n] = Int(length(uPeakM1)) + Int(length(uPeakM2))
                    @views solMArr[n] = sol.u[end][:, 1]
                    @views solEArr[n] = sol.u[end][:, 2]
					# println()
					# println(sol.u[end])

					# DEBUG: Save entire track
					if(saveDists == 1)
						println(folder * "/$(nSeed)_$(n)_SweepData.jld2")
						save_object(folder * "/$(nSeed)_$(n)_SweepData.jld2", [tSS[n], sol.t, sol.u])
					end

				elseif(nessFlag == 1) # Indicate NESS was reached (NESS and linear growth states are mutually exclusive)
					println("System reached NESS...")
					ssReached[n] = -1
				end

			else
				# println("Could not complete simulation in allocated number of iterations")
				ssReached[n] = -3
				break
			end

		catch err1

			# Check type of error
			println("\nTypes of error:\n$(err1)")

			# Print current variables for debugging
			println("\nCurrent system parameters are:")
			println("n = $(n)")
			println(pNew)

            # Create outputting error to stop simulation
            error("ERROR")

		end

	end



    ####################################
    # Calculate order parameters       #
    ####################################

	# No systems that did not finish calculations can reach this point (no need to check for ssReached == -3)

    # Check whether many different tissue lengths were simulated
    if(dynamicScalingSim == 1)

        # Not interested in conventional definition of scaling (single vector)
        Sm = SmOriginalCurr
        SmNorm = SmNormOriginalCurr
        Se = SeOriginalCurr

        # Cannot calculate robustness
        Rm = [RmOriginalCurr]
        Re = [ReOriginalCurr]

        # Calculate dynamic morphogen scaling
        @views SmDyn = [isassigned(solMArr, i) && isassigned(solMArr, i + 1) ? 
                            calculateLocalScalingOriginal(solMArr[1], solMArr[i + 1], sIndex1, N1) : [] for i in 1:nLVec - 1]
        #= @inbounds for i in 1:nLVec - 1
            println(SmDyn[i])
        end =#

    else

        # Different scaling and robustness are calculated for different types of simulations
        if(numVary == 1)

            # Calculate scaling
            if(wVary == 1)

                # Calculate morphogen scaling
                @views Sm = isassigned(solMArr, 1) && isassigned(solMArr, 4) ? 
                                calculateLocalScalingOriginal(solMArr[1], solMArr[4], sIndex1, N1) : []
                @views SmNorm = isassigned(solMArr, 1) && isassigned(solMArr, 4) ? 
                                calculateLocalScalingOriginal(solMArr[1] ./ solMArr[1][sIndex1], solMArr[4] ./ solMArr[4][sIndex1], sIndex1, N1) : []

                # Calculate expander scaling
                @views Se = isassigned(solEArr, 1) && isassigned(solEArr, 4) ?
                                calculateLocalScalingOriginal(solEArr[1], solEArr[4], sIndex1, N1) : []

            else # Use different indexing if the source width does not change

                @views Sm = isassigned(solMArr, 1) && isassigned(solMArr, 4) ? 
                                calculateLocalScalingOriginal(solMArr[1], solMArr[4], sWConstIndex1, NWConst1) : []
                @views SmNorm = isassigned(solMArr, 1) && isassigned(solMArr, 4) ? 
                                calculateLocalScalingOriginal(solMArr[1] ./ solMArr[1][sIndex1], solMArr[4] ./ solMArr[4][sWConstIndex1], sWConstIndex1, NWConst1) : []

                @views Se = isassigned(solEArr, 1) && isassigned(solEArr, 4) ?
                                calculateLocalScalingOriginal(solEArr[1], solEArr[4], sWConstIndex1, NWConst1) : []
                                
            end

            # Calculate morphogen robustness
            @views Rm = [isassigned(solMArr, 1) && isassigned(solMArr, 2) && isassigned(solMArr, 3) ?
                            calculateLocalRobustnessOriginal(solMArr[1], solMArr[2], solMArr[3]) : []]

            # Calculate expander robustness
            @views Re = [isassigned(solEArr, 1) && isassigned(solEArr, 2) && isassigned(solEArr, 3) ?
                            calculateLocalRobustnessOriginal(solEArr[1], solEArr[2], solEArr[3]) : []]

        else

            # Cannot define scaling
            Sm = SmOriginalCurr
            SmNorm = SmNormOriginalCurr
            Se = SeOriginalCurr

            # Calculate morphogen robustness
            @views Rm = [RmOriginalCurr,
                            isassigned(solMArr, 1) && isassigned(solMArr, 2) && isassigned(solMArr, 3) ?
                            calculateLocalRobustnessOriginal(solMArr[1], solMArr[2], solMArr[3]) : [],
                            isassigned(solMArr, 1) && isassigned(solMArr, 4) && isassigned(solMArr, 5) ?
                            calculateLocalRobustnessOriginal(solMArr[1], solMArr[4], solMArr[5]) : [],
                            isassigned(solMArr, 1) && isassigned(solMArr, 6) && isassigned(solMArr, 7) ?
                            calculateLocalRobustnessOriginal(solMArr[1], solMArr[6], solMArr[7]) : [],
                            isassigned(solMArr, 1) && isassigned(solMArr, 8) && isassigned(solMArr, 9) ?
                            calculateLocalRobustnessOriginal(solMArr[1], solMArr[8], solMArr[9]) : [],
                            isassigned(solMArr, 1) && isassigned(solMArr, 10) && isassigned(solMArr, 11) ?
                            calculateLocalRobustnessOriginal(solMArr[1], solMArr[10], solMArr[11]) : []]

            # Calculate expander robustness
            @views Re = [ReOriginalCurr,
                            isassigned(solEArr, 1) && isassigned(solEArr, 2) && isassigned(solEArr, 3) ?
                            calculateLocalRobustnessOriginal(solEArr[1], solEArr[2], solEArr[3]) : [],
                            isassigned(solEArr, 1) && isassigned(solEArr, 4) && isassigned(solEArr, 5) ?
                            calculateLocalRobustnessOriginal(solEArr[1], solEArr[4], solEArr[5]) : [],
                            isassigned(solEArr, 1) && isassigned(solEArr, 6) && isassigned(solEArr, 7) ?
                            calculateLocalRobustnessOriginal(solEArr[1], solEArr[6], solEArr[7]) : [],
                            isassigned(solEArr, 1) && isassigned(solEArr, 8) && isassigned(solEArr, 9) ?
                            calculateLocalRobustnessOriginal(solEArr[1], solEArr[8], solEArr[9]) : [],
                            isassigned(solEArr, 1) && isassigned(solEArr, 10) && isassigned(solEArr, 11) ?
                            calculateLocalRobustnessOriginal(solEArr[1], solEArr[10], solEArr[11]) : []]

        end

		# Cannot define dynamic morphogen scaling
		SmDyn = []

    end

    # Calculate morphogen precision
    @views Pm = isassigned(solMArr, 1) ? calculateLocalPrecision(solMArr[1]) : []

	# Store fractional distribution ranges
	mRange::Vector{Float64} = [isassigned(solMArr, n) ? (solMArr[n][sIndex1Arr[n]] - solMArr[n][NArr[n]]) / solMArr[n][sIndex1Arr[n]] : -1. for n in 1:nLVec]
	eRange::Vector{Float64} = [isassigned(solEArr, n) ? (solEArr[n][NArr[n]] - solEArr[n][sIndex1Arr[n]]) / solEArr[n][NArr[n]] : -1. for n in 1:nLVec]

	# Calculate half-decay lengths
	@views mHDL::Vector{Float64} = [isassigned(solMArr, n) ? calculateMHDL(solMArr[n], dxVec[n], sIndex1Arr[n]) : -1. for n in 1:nLVec]
	@views eHDL::Vector{Float64} = [isassigned(solEArr, n) ? calculateEHDL(solEArr[n], dxVec[n], NArr[n]) : -1. for n in 1:nLVec]

	# Calculate amplitudes
	mAmp0::Vector{Float64} = [isassigned(solMArr, n) ? solMArr[n][1] : -1. for n in 1:nLVec]
	mAmpSWidth::Vector{Float64} = [isassigned(solMArr, n) ? solMArr[n][sIndex1Arr[n]] : -1. for n in 1:nLVec]
	mAmpL::Vector{Float64} = [isassigned(solMArr, n) ? solMArr[n][end] : -1. for n in 1:nLVec]
	mMaxRatio::Vector{Float64} = [isassigned(solMArr, n) ? maximum(solMArr[n]) / solMArr[n][1] : -1. for n in 1:nLVec]
	eAmp::Vector{Float64} = [isassigned(solEArr, n) ? solEArr[n][NArr[n]] : -1. for n in 1:nLVec]

	# Store expander source width size
	@views eWidth::Vector{Float64} = reduce(vcat, [isassigned(solEArr, n) ? calculateESourceWidth(Vector{Float64}(solMArr[n]), (thresh * delThreshVec[n]) ^ h, LVec[n], dxVec[n], NArr[n]) : [-1., -1.] for n in 1:nLVec])

	# Store number of grid points where expander concentration is higher than Hill constants
	@views eHighFrac::Vector{Int} = [isassigned(solEArr, n) ? Int(sum(solEArr[n][sIndex1Arr[n]:NArr[n]] .> max(delM, delE))) : -1. for n in 1:nLVec]

    # Calculate the average expander degradation rate
    @views muAv = [isassigned(solEArr, n) ? sum(pNew[7] ./ (delE .+ solEArr[n][sIndex1Arr[n]:NArr[n]])) / (NArr[n] - sIndex1Arr[n] + 1) : -1. for n in 1:nLVec];

    # Fit morphogen gradient in target tissue
    chiSqRat::Float64 = -1e6
    if((chiSqRatCalc == 1) && (isassigned(solMArr, 1)))

        # Extract normalised morphogen gradient in target tissue
        @views mNorm = solMArr[1][sIndex1:N] ./ solMArr[1][sIndex1]

        # Fit with exponential distribution and calculate chi-squared value
        expFit = mHDL[1] != Inf ? LsqFit.curve_fit(mFitExp, xL1_sIndex, mNorm, [0.1, sWidth[1], mHDL[1], 1.]) : LsqFit.curve_fit(mFitExp, xL1_sIndex, mNorm, [0.1, sWidth[1], LBuff_2, 1.])
		# println(expFit.param)
        chiSqExp::Float64 = sum(((mFitExp(xL1_sIndex, expFit.param) .- mNorm) .^ 2) ./ mNorm);

        # Fit with power-law distribution and calculate chi-squared value
        powFit = LsqFit.curve_fit(mFitPow, xL1_sIndex, mNorm, [8e-5, 5., 1.5, 1.]; lower = [0., -10., 0., 0.], upper = [10., 9., Inf, 10.])
		# println(powFit.param)
        chiSqPow::Float64 = sum(((mFitPow(xL1_sIndex, powFit.param) .- mNorm) .^ 2) ./ mNorm);

        # Calculate fit ratio
        if((!isinf(chiSqExp)) && (!isinf(chiSqPow)))
            chiSqRat = log10(chiSqExp / chiSqPow)
			# println(chiSqRat)
         end

	end

	return SimData(nSeed, 
					pNew, 
					ssReached, 
					tSS, 
					# nPeaks, 
					# solMArr, 
					# solEArr, 
					Sm,
					SmNorm,
					SmDyn, 
					Rm,
                    Pm,
					mRange, 
					mHDL,
					mAmp0,
					mAmpSWidth,
					mAmpL,
					mMaxRatio,
					Se, 
					Re,
					eRange, 
					eHDL,
					eAmp,
					eWidth,
					eHighFrac,
                    muAv,
					chiSqRat)

end

# Define the discretized PDE with feedback through the morphogen degradation rate, including within the morphogen source region, as an ODE function
function fKFeedback!(du, u, p, t)
	dm, de, num, nue, nueThreshPow, k1, mu = p

	du[1, 1] = (2 * dm * u[2, 1]) - (((2 * dm) + (k1 / (delM + u[1, 2]))) * u[1, 1]) + num
	du[1, 2] = (2 * de * u[2, 2]) - (((2 * de) + (mu / (delE + u[1, 2]))) * u[1, 2]) + (nueThreshPow / (threshPow + (u[1, 1] ^ h)))
	@inbounds for i in 2:sIndex
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1]) + num
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	@inbounds for i in sIndex1:N1
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1])
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	du[N, 1] = (2 * dm * u[N1, 1]) - (((2 * dm) + (k1 / (delM + u[N, 2]))) * u[N, 1])
	du[N, 2] = (2 * de * u[N1, 2]) - (((2 * de) + (mu / (delE + u[N, 2]))) * u[N, 2]) + (nueThreshPow / (threshPow + (u[N, 1] ^ h)))

end
function fKFeedback_threshVary!(du, u, p, t)
	dm, de, num, nue, nueThreshPow, k1, mu = p

	du[1, 1] = (2 * dm * u[2, 1]) - (((2 * dm) + (k1 / (delM + u[1, 2]))) * u[1, 1]) + num
	du[1, 2] = (2 * de * u[2, 2]) - (((2 * de) + (mu / (delE + u[1, 2]))) * u[1, 2]) + (nueThreshPow / (currThreshPow + (u[1, 1] ^ h)))
	@inbounds for i in 2:sIndex
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1]) + num
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (currThreshPow + (u[i, 1] ^ h)))
	end
	@inbounds for i in sIndex1:N1
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1])
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (currThreshPow + (u[i, 1] ^ h)))
	end
	du[N, 1] = (2 * dm * u[N1, 1]) - (((2 * dm) + (k1 / (delM + u[N, 2]))) * u[N, 1])
	du[N, 2] = (2 * de * u[N1, 2]) - (((2 * de) + (mu / (delE + u[N, 2]))) * u[N, 2]) + (nueThreshPow / (currThreshPow + (u[N, 1] ^ h)))

end
function fKFeedbackWConst!(du, u, p, t)
	dm, de, num, nue, nueThreshPow, k1, mu = p

	du[1, 1] = (2 * dm * u[2, 1]) - (((2 * dm) + (k1 / (delM + u[1, 2]))) * u[1, 1]) + num
	du[1, 2] = (2 * de * u[2, 2]) - (((2 * de) + (mu / (delE + u[1, 2]))) * u[1, 2]) + (nueThreshPow / (threshPow + (u[1, 1] ^ h)))
	@inbounds for i in 2:sWConstIndex
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1]) + num
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	@inbounds for i in sWConstIndex1:NWConst1
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1])
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	du[NWConst, 1] = (2 * dm * u[NWConst1, 1]) - (((2 * dm) + (k1 / (delM + u[NWConst, 2]))) * u[NWConst, 1])
	du[NWConst, 2] = (2 * de * u[NWConst1, 2]) - (((2 * de) + (mu / (delE + u[NWConst, 2]))) * u[NWConst, 2]) + (nueThreshPow / (threshPow + (u[NWConst, 1] ^ h)))

end
function fKFeedbackOriginal!(du, u, p, t)
	dm, de, num, nue, nueThreshPow, k1, mu = p

	du[1, 1] = (2 * dm * u[2, 1]) - (((2 * dm) + (k1 / (delM + u[1, 2]))) * u[1, 1]) + num
	du[1, 2] = (2 * de * u[2, 2]) - (((2 * de) + (mu / delE)) * u[1, 2]) + (nueThreshPow / (threshPow + (u[1, 1] ^ h)))
	@inbounds for i in 2:sIndex
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1]) + num
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / delE)) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	@inbounds for i in sIndex1:N1
		du[i, 1] = (dm * (u[i - 1, 1] + u[i + 1, 1])) - (((2 * dm) + (k1 / (delM + u[i, 2]))) * u[i, 1])
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / delE)) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	du[N, 1] = (2 * dm * u[N1, 1]) - (((2 * dm) + (k1 / (delM + u[N, 2]))) * u[N, 1])
	du[N, 2] = (2 * de * u[N1, 2]) - (((2 * de) + (mu / delE)) * u[N, 2]) + (nueThreshPow / (threshPow + (u[N, 1] ^ h)))

end

# Define the discretized PDE with feedback through the morphogen diffusivity, including within the morphogen source region, as an ODE function
function fDFeedback!(du, u, p, t)
	dm, de, num, nue, nueThreshPow, k1, mu = p

	diffX = [dm * (delM + u[i, 2]) for i in 1:N]

	# du[1, 1] = (2 * dm * (delM + u[1, 2]) * u[2, 1]) - (((2 * dm * (delM + u[1, 2])) + k1) * u[1, 1]) + num
	du[1, 1] = (2 * diffX[1] * u[2, 1]) - (((2 * diffX[1]) + k1) * u[1, 1]) + num
	du[1, 2] = (2 * de * u[2, 2]) - (((2 * de) + (mu / (delE + u[1, 2]))) * u[1, 2]) + (nueThreshPow / (threshPow + (u[1, 1] ^ h)))
	@inbounds for i in 2:sIndex
		# du[i, 1] = (((dm * (delM + u[i + 1, 2])) - (dm * (delM + u[i - 1, 2]))) * (u[i + 1, 1] - u[i - 1, 1])) + ((dm * (delM + u[i, 2])) * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1]) + num
		du[i, 1] = ((diffX[i + 1] - diffX[i - 1]) * (u[i + 1, 1] - u[i - 1, 1])) + (diffX[i] * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1]) + num
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	@inbounds for i in sIndex1:N1
		# du[i, 1] = (((dm * (delM + u[i + 1, 2])) - (dm * (delM + u[i - 1, 2]))) * (u[i + 1, 1] - u[i - 1, 1])) + ((dm * (delM + u[i, 2])) * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1])
		du[i, 1] = ((diffX[i + 1] - diffX[i - 1]) * (u[i + 1, 1] - u[i - 1, 1])) + (diffX[i] * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1])
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	# du[N, 1] = (2 * dm * (delM + u[N, 2]) * u[N1, 1]) - (((2 * dm * (delM + u[N, 2])) + k1) * u[N, 1])
	du[N, 1] = (2 * diffX[N] * u[N1, 1]) - (((2 * diffX[N]) + k1) * u[N, 1])
	du[N, 2] = (2 * de * u[N1, 2]) - (((2 * de) + (mu / (delE + u[N, 2]))) * u[N, 2]) + (nueThreshPow / (threshPow + (u[N, 1] ^ h)))

end
function fDFeedback_threshVary!(du, u, p, t)
	dm, de, num, nue, nueThreshPow, k1, mu = p

	diffX = [dm * (delM + u[i, 2]) for i in 1:N]

	# du[1, 1] = (2 * dm * (delM + u[1, 2]) * u[2, 1]) - (((2 * dm * (delM + u[1, 2])) + k1) * u[1, 1]) + num
	du[1, 1] = (2 * diffX[1] * u[2, 1]) - (((2 * diffX[1]) + k1) * u[1, 1]) + num
	du[1, 2] = (2 * de * u[2, 2]) - (((2 * de) + (mu / (delE + u[1, 2]))) * u[1, 2]) + (nueThreshPow / (currThreshPow + (u[1, 1] ^ h)))
	@inbounds for i in 2:sIndex
		# du[i, 1] = (((dm * (delM + u[i + 1, 2])) - (dm * (delM + u[i - 1, 2]))) * (u[i + 1, 1] - u[i - 1, 1])) + ((dm * (delM + u[i, 2])) * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1]) + num
		du[i, 1] = ((diffX[i + 1] - diffX[i - 1]) * (u[i + 1, 1] - u[i - 1, 1])) + (diffX[i] * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1]) + num
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (currThreshPow + (u[i, 1] ^ h)))
	end
	@inbounds for i in sIndex1:N1
		# du[i, 1] = (((dm * (delM + u[i + 1, 2])) - (dm * (delM + u[i - 1, 2]))) * (u[i + 1, 1] - u[i - 1, 1])) + ((dm * (delM + u[i, 2])) * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1])
		du[i, 1] = ((diffX[i + 1] - diffX[i - 1]) * (u[i + 1, 1] - u[i - 1, 1])) + (diffX[i] * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1])
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (currThreshPow + (u[i, 1] ^ h)))
	end
	# du[N, 1] = (2 * dm * (delM + u[N, 2]) * u[N1, 1]) - (((2 * dm * (delM + u[N, 2])) + k1) * u[N, 1])
	du[N, 1] = (2 * diffX[N] * u[N1, 1]) - (((2 * diffX[N]) + k1) * u[N, 1])
	du[N, 2] = (2 * de * u[N1, 2]) - (((2 * de) + (mu / (delE + u[N, 2]))) * u[N, 2]) + (nueThreshPow / (currThreshPow + (u[N, 1] ^ h)))

end
function fDFeedbackWConst!(du, u, p, t)
	dm, de, num, nue, nueThreshPow, k1, mu = p

	diffX = [dm * (delM + u[i, 2]) for i in 1:N]

	# du[1, 1] = (2 * dm * (delM + u[1, 2]) * u[2, 1]) - (((2 * dm * (delM + u[1, 2])) + k1) * u[1, 1]) + num
	du[1, 1] = (2 * diffX[1] * u[2, 1]) - (((2 * diffX[1]) + k1) * u[1, 1]) + num
	du[1, 2] = (2 * de * u[2, 2]) - (((2 * de) + (mu / (delE + u[1, 2]))) * u[1, 2]) + (nueThreshPow / (threshPow + (u[1, 1] ^ h)))
	@inbounds for i in 2:sWConstIndex
		# du[i, 1] = (((dm * (delM + u[i + 1, 2])) - (dm * (delM + u[i - 1, 2]))) * (u[i + 1, 1] - u[i - 1, 1])) + ((dm * (delM + u[i, 2])) * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1]) + num
		du[i, 1] = ((diffX[i + 1] - diffX[i - 1]) * (u[i + 1, 1] - u[i - 1, 1])) + (diffX[i] * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1]) + num
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	@inbounds for i in sWConstIndex1:NWConst1
		# du[i, 1] = (((dm * (delM + u[i + 1, 2])) - (dm * (delM + u[i - 1, 2]))) * (u[i + 1, 1] - u[i - 1, 1])) + ((dm * (delM + u[i, 2])) * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1])
		du[i, 1] = ((diffX[i + 1] - diffX[i - 1]) * (u[i + 1, 1] - u[i - 1, 1])) + (diffX[i] * (u[i + 1, 1] + u[i - 1, 1] - (2 * u[i, 1]))) - (k1 * u[i, 1])
		du[i, 2] = (de * (u[i - 1, 2] + u[i + 1, 2])) - (((2 * de) + (mu / (delE + u[i, 2]))) * u[i, 2]) + (nueThreshPow / (threshPow + (u[i, 1] ^ h)))
	end
	# du[NWConst, 1] = (2 * dm * (delM + u[NWConst, 2]) * u[NWConst1, 1]) - (((2 * dm * (delM + u[NWConst, 2])) + k1) * u[NWConst, 1])
	du[NWConst, 1] = (2 * diffX[NWConst] * u[NWConst1, 1]) - (((2 * diffX[NWConst]) + k1) * u[NWConst, 1])
	du[NWConst, 2] = (2 * de * u[NWConst1, 2]) - (((2 * de) + (mu / (delE + u[NWConst, 2]))) * u[NWConst, 2]) + (nueThreshPow / (threshPow + (u[NWConst, 1] ^ h)))

end

# Define callback condition function that checks for when the derivative of the total morphogen concentration is zero
@views oscCheckMCondition(u, t, integrator) = trapezium(integrator(t, Val{1})[:, 1])

# Define affect of successful callback condition that stores maxima and minima parameters and checks if system has reached NESS
function oscCheckMAffect!(integrator)

	# Extract the total morphogen concentration
	@views mCurr::Float64 = trapezium(integrator.u[:, 1])

	# Check if enough time has passed since last maximum/minimum detection and if the observed concentration has changed sufficiently
	# (Prevents errors since adaptive timestep -> 0 near morphogen maximum/minimum due to maximum rate of change of expander concentration)
	if((integrator.t > tPeakMPrev) && ((max(mCurr, uPeakMPrev) / min(mCurr, uPeakMPrev)) - 1. > relTol))
		# println("HERE: $(integrator.t): $(mCurr), $(uPeakMPrev), $((max(mCurr, uPeakMPrev) / min(mCurr, uPeakMPrev)) - 1.)")

        # Julia does not have functionality to use second order derivative to identify maxima/minima due to adaptive timestepping
        # Instead alternate the arrays that the peaks are stored in, assuming inflexion points are not possible (but check for these manually)

        # Store current time and peak amplitude
        global tPeakMPrev = Float64(integrator.t + 0.1)
        global uPeakMPrev = mCurr

        # Store peak amplitude in alternating arrays
        if(typePeakMPrev == 0)
            push!(uPeakM1, uPeakMPrev)
            # println("Peak 1 = ($(tPeakMPrev), $(uPeakMPrev))")
        else
            push!(uPeakM2, uPeakMPrev)
            # println("Peak 2 = ($(tPeakMPrev), $(uPeakMPrev))")
        end

        # Flip flag indicating which peak array to store data in
        global typePeakMPrev = 1 - typePeakMPrev

        # Check if equal (non-zero) numbers of maxima and minima have been detected so far
        nPeakM1::Int = length(uPeakM1)
        if((nPeakM1 > 1) && (nPeakM1 == length(uPeakM2)))

            # Check if net translation occurred
            nPeakM11::Int = nPeakM1 - 1
            #= if(((uPeakM1[nPeakM1] - uPeakM1[nPeakM11] > 0.1) && (uPeakM2[nPeakM2] - uPeakM2[nPeakM21] > 0.1)) || 
                ((uPeakM1[nPeakM1] - uPeakM1[nPeakM11] < -0.1) && (uPeakM2[nPeakM2] - uPeakM2[nPeakM21] < -0.1)))
                # println("Net translation detected\n")

                # Clear arrays containing stored maxima and minima times and magnitudes
                # empty!(uPeakM1)
                # empty!(uPeakM2)

            else =#

                # Compare ratios between peak heights in time to check if system has reached NESS
                magCurr::Vector{Float64} = [abs(uPeakM1[i] - uPeakM2[i]) for i in 1:nPeakM1]
                for i in 1:nPeakM11
        
                    if(abs(min(magCurr[end] / magCurr[i], magCurr[i] / magCurr[end]) - 1) < 0.005)
                        # println("\nSystem reached morphogen NESS\n")

                        # Terminate solver
                        global nessFlag = 1

                        # Moved termination of integrator to newTerminateSTeadyStateCondition() function

                    end
        
                end

            # end

        else

            # Check for inflexion points: number of detected maxima should always be within number of minima +/- 1
            if(abs(nPeakM1 - length(uPeakM2)) > 1)
                # println("Inflexion point detected\n")

                # Clear arrays containing stored maxima and minima times and magnitudes
                empty!(uPeakM1)
                empty!(uPeakM2)

            end

        end

	end

end

# Define updated function that checks for when simulation reaches steady state
function newTerminateSteadyStateCondition(u, t, integrator)

	# Extract distribution and derivative information
	intFromCache = first(get_tmp_cache(integrator))
    DiffEqBase.get_du!(intFromCache, integrator)
	# println("$(t), $(trapezium(integrator.u[:, 1])), $(trapezium(integrator.u[:, 2])) ($(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 1], u[:, 1]))), $(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 2], u[:, 2]))))")

	# First check if distributions have reached steady state
	@views if(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 1], u[:, 1])) ||
			  any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 2], u[:, 2])))

        # Check if system is in an NESS state
        if(nessFlag == 1)
            # println("\nSystem reached morphogen NESS at check 1 at t = $(t)\n")
            terminate!(integrator)
        end
		
        # Else reset steady state flags and times
        #= if(flagSS > 0)
			println("\nSimulated $(flagSS) points before deviation at t = $(t), Mtot = $(trapezium(integrator.u[:, 1])) ($(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 1], u[:, 1]))), $(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 2], u[:, 2]))))\n")
		end =#
		# println("Check 1 at t = $(t) ($(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 1], u[:, 1]))), $(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 2], u[:, 2]))))")
        global flagSS = 0
		global tAnySS = 0.
		global tPrevSS = 0.
		global mTotPrevSS = 0. # These two definitions are technically not required here
		global mTotDevMax = 0.
		return false
		
	# If system reached steady state check if this is the first time
	elseif(flagSS == 0)

		# Change flag to indicate steady state has been reached and store current time
		# println("\nCheck 2 at t = $(t), Mtot = $(trapezium(integrator.u[:, 1])) ($(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 1], u[:, 1]))), $(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 2], u[:, 2]))))\n")
		global flagSS = 1
		global tAnySS = t
		global tPrevSS = t
		@views global mTotPrevSS = trapezium(integrator.u[:, 1])
		global mTotDevMax = copy(mTotPrevSS)
		return false

	# If system has reached steady state but not for long enough yet
	elseif((t < tPrevSS + ssCheckTime) || (flagSS < nSSCheck))

		# println("\nCheck 3 at t = $(t), Mtot = $(trapezium(integrator.u[:, 1])) ($(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 1], u[:, 1]))), $(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 2], u[:, 2]))))\n")
		global flagSS += 1
		@views if(abs(trapezium(integrator.u[:, 1]) - mTotPrevSS) > abs(mTotDevMax - mTotPrevSS))
			@views global mTotDevMax = trapezium(integrator.u[:, 1])
		end
		return false

	# If system has reached steady state for a prolonged period but has varied slowly during this time so as not to be detected
	elseif((t >= tPrevSS + ssCheckTime) && (flagSS >= nSSCheck) && (abs(mTotDevMax - mTotPrevSS) / mTotPrevSS > 0.01))
        
        # Check if system is in an NESS state
        if(nessFlag == 1)
            # println("\nSystem reached morphogen NESS at check 4 at t = $(t)\n")
            terminate!(integrator)
        end

        # Else restart steady-state checks
		# println("\nRestarting steady-state checks at t = $(t), Mtot = $(trapezium(integrator.u[:, 1])) ($(abs(mTotDevMax - mTotPrevSS) / mTotPrevSS) > 0.01)\n")
		# println("\nCheck 4 at t = $(t), Mtot = $(trapezium(integrator.u[:, 1])) ($(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 1], u[:, 1]))), $(any((abs(deriv) > absTol) && (abs(deriv) > relTol * abs(val)) for (deriv, val) in zip(intFromCache[:, 2], u[:, 2]))))\n")
		global flagSS = 1
		global tPrevSS = t
		@views global mTotPrevSS = trapezium(integrator.u[:, 1])
		global mTotDevMax = copy(mTotPrevSS)
		return false
	
	# If system has reached steady state for a prolonged period
	elseif((t >= tPrevSS + ssCheckTime) && (flagSS >= nSSCheck))

		# println("\nCheck 5 at t = $(t) > $(tPrevSS + ssCheckTime), Mtot = $(trapezium(integrator.u[:, 1])) ($(abs(mTotDevMax - mTotPrevSS) / mTotPrevSS) < 0.01)\n")
		return true

	end

end

# Define updated function to terminate the solver when steady state has been reached
newTerminateSteadyStateAffect!(integrator) = terminate!(integrator)

# Define function to calculate an integral using the trapezium rule
@views trapezium(u) = u[1] + u[end] + (2 * sum(u[2:end - 1])) # Only used for in-code, dimensionless functions, to calculate actual integrals need factor of dx/2

# Define function to calculate the average scaling around specific points using the same measure as in BenZvi2010Scaling
function calculateLocalScalingOriginal(arrL1::Vector{Float64}, arrL2::Vector{Float64}, sInd1::Int, NInd1::Int)

    # Initialise array to store the average local scaling around each pre-specified point
    delPos = Vector{Float64}(undef, nBIndexes)

    # Extract local concentrations pre-growth
    @views valL1::Vector{Float64} = arrL1[bIndexes]

    # Loop over pre-specified points
    @inbounds for i in 1:nBIndexes

        # Locate boundary position in new array
        indCurr::Int = 0
        @inbounds for j in sInd1:NInd1
            if(((arrL2[j] >= valL1[i]) && (arrL2[j + 1] <= valL1[i])) || ((arrL2[j] <= valL1[i]) && (arrL2[j + 1] >= valL1[i])))
                indCurr = j
                break
            end
        end
        # println(indCurr)

        # Perform linear interpolation
        if(indCurr != 0)

            m = (arrL2[indCurr + 1] - arrL2[indCurr]) / dxL2
            c = ((arrL2[indCurr] + arrL2[indCurr + 1]) - (m * (xL2[indCurr] + xL2[indCurr + 1]))) / 2.

            # Calculate absolute difference between positions before and after change
            delPos[i] = (valL1[i] - c) / m
            # println(delPos[i])
            # delPos[i] = abs((delPos[i] / xL2[end]) - posL1[i])
			delPos[i] = abs(((delPos[i] - xL2SWidth) / xL2End_SWidth) - posL1[i])
            # println(delPos[i])

        else
            
            # System cannot scale
            delPos[i] = -1.

        end

    end

	# Return the local changes in position
	return delPos

end;

# Define function to calculate the average robustness around specific points using the same measure as in BenZvi2010Scaling
function calculateLocalRobustnessOriginal(arrNu1::Vector{Float64}, arrNu2::Vector{Float64}, arrNu3::Vector{Float64})

    # Initialise array to store the average local scaling around each pre-specified point
    delPos = Vector{Float64}(undef, nBIndexes)

    # Extract local concentrations pre-growth
    @views valNu1::Vector{Float64} = arrNu1[bIndexes]

    # Loop over pre-specified points
    @inbounds for i in 1:nBIndexes

        # Locate boundary position for decreased nu (2)
        indCurr::Int = 0
        @inbounds for j in sIndex1:N1
            if(((arrNu2[j] >= valNu1[i]) && (arrNu2[j + 1] <= valNu1[i])) || ((arrNu2[j] <= valNu1[i]) && (arrNu2[j + 1] >= valNu1[i])))
                indCurr = j
                break
            end
        end
        # println(indCurr)

        # Perform linear interpolation
        posNu2::Float64 = 0.
        if(indCurr != 0)

            m = (arrNu2[indCurr + 1] - arrNu2[indCurr]) / dxL1
            c = ((arrNu2[indCurr] + arrNu2[indCurr + 1]) - (m * (xL1[indCurr] + xL1[indCurr + 1]))) / 2.

            # Calculate absolute difference between positions before and after change
            posNu2 = (valNu1[i] - c) / m
            # println(posNu2)
            # posNu2 = abs((posNu2 / xL1End) - posL1[i])
			posNu2 = abs(((posNu2 - xL1SWidth) / xL1End_SWidth) - posL1[i])
            # println(posNu2)

        else
            
            # System cannot scale
            posNu2 = -1.

        end

        # Repeat calculations for increased nu (3)
        indCurr = 0
        @inbounds for j in 1:N1
            if(((arrNu3[j] >= valNu1[i]) && (arrNu3[j + 1] <= valNu1[i])) || ((arrNu3[j] <= valNu1[i]) && (arrNu3[j + 1] >= valNu1[i])))
                indCurr = j
                break
            end
        end
        # println(indCurr)

        # Perform linear interpolation
        posNu3::Float64 = 0.
        if(indCurr != 0)

            m = (arrNu3[indCurr + 1] - arrNu3[indCurr]) / dxL1
            c = ((arrNu3[indCurr] + arrNu3[indCurr + 1]) - (m * (xL1[indCurr] + xL1[indCurr + 1]))) / 2.

            # Calculate absolute difference between positions before and after change
            posNu3 = (valNu1[i] - c) / m
            # println(posNu3)
            # posNu3 = abs((posNu3 / xL1End) - posL1[i])
			posNu3 = abs(((posNu3 - xL1SWidth) / xL1End_SWidth) - posL1[i])
            # println(posNu3)

        else
            
            # System cannot scale
            posNu3 = -1.

        end

        # Store average change in boundary position
        delPos[i] = (posNu2 .== -1.) || (posNu3 .== -1.) ? -1. : (posNu2 + posNu3) / 2.

    end

	# Return the local changes in position
	return delPos

end;

# Define function to calculate the average scaling around specific points using the same measure as in BenZvi2010Scaling
function calculateLocalPrecision(arr::Vector{Float64})
    
    # Initialise array to store the precision at each pre-specified point
    prec = Vector{Float64}(undef, nBIndexes)

    # Loop over pre-specified points
	grad::Float64 = 0.
    @inbounds for i in 1:nBIndexes

        # Calculate the local morphogen gradient
        if(bIndexes[i] < N)
            grad = (arr[bIndexes[i] + 1] - arr[bIndexes[i] - 1]) / dr2
        else
            grad = (arr[bIndexes[i]] - arr[bIndexes[i] - 1]) / dr
        end

        # Calculate the local precision
        prec[i] = 1. - abs((sqrt(arr[bIndexes[i]] / LCell) / grad))

    end

    # Return the precision at each position
    return prec

end;

# Define function to calculate expander source width
function calculateESourceWidth(arr::Vector{Float64}, currThreshPow::Float64, LBuff::Float64, dxBuff::Float64, NInd::Int)

	# Loop over target tissue
	sourceWidthIndex::Int = 0;
	@inbounds for i in 2:NInd

		# Find where morphogen concentration crosses ETHRESH
		if((arr[i - 1] > thresh) && (arr[i] < thresh))
			sourceWidthIndex = i
			break
		end

	end

	# Morphogen concentration crossed ETHRESH within target tissue
	width1::Float64 = 0.
	if(sourceWidthIndex > 0)

		# Use linear interpolation to find more accurate position
		x2 = Float64(sourceWidthIndex - 1) * dxBuff
		y2 = arr[sourceWidthIndex]
		m = (y2 - arr[sourceWidthIndex - 1]) / dxBuff
		c = y2 - (m * x2)
		pos = (thresh - c) / m

		# Calculate expander source width
		width1 = LBuff - pos

	end # Else morphogen concentration did not cross ETHRESH within tissue, and width is already set to 0

	# Also calculate integrated expander production term
	width2::Float64 = (dxBuff / 2.) * trapezium(currThreshPow ./ (currThreshPow .+ (arr .^ h)))

	return [width1, width2]

end

# Define function to calculate morphogen half-decay length (within the target tissue)
function calculateMHDL(arr::Vector{Float64}, dxBuff::Float64, sInd1::Int)

	# Find index closest to position of half-decay length
	target::Float64 = arr[sInd1] / 2.
	@views hdlIndex::Vector{Int} = findall(x -> (x < target), arr[sInd1:end])
	
	# Check if half-decay length exists
	if(!isempty(hdlIndex))

		# Use linear interpolation to find more accurate position
		currIndex::Int = sInd1 - 1 + hdlIndex[1]
		x2 = currIndex * dxBuff
		y2 = arr[currIndex]
		m = (y2 - arr[currIndex - 1]) / dxBuff
		c = y2 - (m * x2)
		pos = (target - c) / m
	
		return pos

	else
		return Inf # Half-decay length is not measurable in target tissue
	end

end

# Define function to calculate expander half-decay length (within the target tissue)
function calculateEHDL(arr::Vector{Float64}, dxBuff::Float64, NInd::Int)

	# Find index closest to position of half-decay length
	target::Float64 = arr[NInd] / 2.
	@views hdlIndex::Vector{Int} = findall(x -> (x > target), arr)
	
	# Check if half-decay length exists
	if(!isempty(hdlIndex) && (hdlIndex[1] > 1))

		# Use linear interpolation to find more accurate position
		currIndex::Int = hdlIndex[1]
		x2 = currIndex * dxBuff
		y2 = arr[currIndex]
		m = (y2 - arr[currIndex - 1]) / dxBuff
		c = y2 - (m * x2)
		pos = (target - c) / m
	
		return pos

	else
		return Inf # Half-decay length is not measurable in target tissue
	end

end
