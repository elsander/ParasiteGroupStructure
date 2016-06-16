#### functions to calculate the imbalance in group representation betwen classes
#### of nodes and the probability of obtaining that imbalance, given a grouping

## Install package if not already installed (slow)
#try Pkg.installed("Iterators") catch Pkg.add("Iterators") end
using Iterators

## Calculate the imbalance of a given class distribution
Imbalance(presvects...) = prod([
        maximum([ vect[ii] for vect in presvects ]) /
        sum([ vect[ii] for vect in presvects ])
    for ii in 1:length(presvects[1]) ])

## Calculate the pvalue of a given imbalance
Pvalue(presvects...) = Float64(prod([
        prod([ binomial(BigInt(sum(vect) - sum(vect[1:(ii - 1)])),
                        BigInt(vect[ii]))
               for vect in presvects ]) /
        binomial(BigInt(sum(sum(presvects)) -
                            sum(sum([ vect[1:(ii - 1)]
                                for vect in presvects ]))),
                 BigInt(sum([ vect[ii] for vect in presvects ])))
    for ii in 1:length(presvects[1]) ]))

## Go through all possible distributions for a given class of species across the
## groups, only keeping those which have the right number of species
function GetPotentialPresences(groupsizes, numspp)
    vects = [0:min(x, numspp) for x in groupsizes]
    presvect = []
    for p in product(vects...)
        if sum(p) != numspp
            continue
        end
        if length(presvect) == 0
            presvect = p
        else
            presvect = [presvect p]
        end
    end
    return presvect
end

## Go through all possible combinations of class distributions, only keeping
## those which do not violate the group sizes (This is the slowest step...)
function GetCompatiblePresences(groupsizes, presvects...)
    ### HARD-CODED ALLOCATION STEP SIZE ###
    STEP = 10000
    ### # # # # # # # # # # # # # # # # ###
    fullpres = Array(Any, STEP)
    counter = 0
    for p in product(presvects...)
        if +([collect(q) for q in p]...) != groupsizes
            continue
        end
#        if length(fullpres) == 0
#            fullpres = [collect(q) for q in p]
#        else
#            fullpres = [fullpres [collect(q) for q in p]]
#        end
        counter += 1
        fullpres[counter] = [collect(q) for q in p]
        if counter == length(fullpres)
            push!(fullpres, Array(Any, STEP))
        end
    end
    return fullpres[1:counter]
end

## Construct a table indicating the pvalue of each possible imbalance value
function GetAllImbalances(presvects)
    imbalances = [Imbalance(presvects[ii]...) for ii in 1:size(presvects, 1)]
    pvalues = [Pvalue(presvects[ii]...) for ii in 1:size(presvects, 1)]
    ## Sanity check
    if sum(pvalues) < .99999999
        error("P-values do not sum to unity: $(round(sum(pvalues), 10))")
    elseif sum(pvalues) > 1.00000001
        error("P-values sum to value greater than unity: $(round(sum(pvalues), 10))")
    end
    return hcat(imbalances, pvalues)
end

## Calculate the pvalue of the empirically observed distribution by summing the
## pvalues of all imbalance values >= to that observed
Emppvalue(empirical, imbatable) = sum(imbatable[Imbalance(empirical...) .<= imbatable[:,1],2])

## Sample several distributions imbalance values to infer the empirical pvalue
function SampleDists(cc, psi, groupvect, classvect, numsamps, outfile)
    imbavals = Vector{Float64}(numsamps)
    uniqClasses = unique(classvect)
    uniqGroups = sort(unique(groupvect))
    randdist = [zeros(Int16, length(uniqGroups)) for xx in 1:length(uniqClasses)]
    for nn in 1:numsamps
        ## equivalent to shuffle groups or classes; groups are integers, so less memory to move
        tmpgroup = shuffle(groupvect)
        for ii in 1:length(uniqClasses)
            for jj in uniqGroups
                ## add one to jj because groups start at 0
                randdist[ii][jj+1] = length(find(yy -> yy == jj,
                                             tmpgroup[find(xx -> xx == uniqClasses[ii],
                                                           classvect)]))
            end 
        end
        ## record its imbalance
        imbavals[nn] = Imbalance(randdist...)
        if (numsamps < 100) 
            greater = 0
            for imba in imbavals[1:nn]
                if imba >= psi
                    greater = greater + 1
                end
            end
            prob = greater / nn
            write(outfile, "$cc,$psi,$prob,$nn\n")
        elseif (nn % div(numsamps, 100)) == 0
            greater = 0
            for imba in imbavals[1:nn]
                if imba >= psi
                    greater = greater + 1
                end
            end
            prob = greater / nn
            write(outfile, "$cc,$psi,$prob,$nn\n")
        end
    end
    return nothing
end

