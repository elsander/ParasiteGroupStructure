## Import imbalance functions
include("Imbalance.jl")

srand(123456)

## Calculate the imbalance and probability of imbalance for given group/class vectors 
datafile = ARGS[1]
numsamps = parse(Int, ARGS[2])
cc = ARGS[3]
outdir = ARGS[4]

## ensure trailing slash
if outdir[end] != '/'
    outdir = outdir * "/"
end

## read in the class/group vectors
data = readdlm(datafile, Any)
classvect = data[:,1]
groupvect = data[:,2]
## get matrix size
numspp = length(groupvect)
## Sanity check
if numspp != length(classvect)
    error("Group and class vectors differ in length")
end
## for comparison to earlier analysis, reduce to two classes (cc and not-cc)
if cc != "all" 
    classvect[find(xx -> xx != cc, classvect)] = "Not-" * cc
    numclasses = 2
else
    numclasses = length(unique(classvect))
end
## count the number of species in each class
sppcounts = zeros(Int16, numclasses)
for ii in classvect
    sppcounts[findfirst(xx -> xx == ii, unique(classvect))] += 1
end
## get group dimensions/distribution
## (note that the groupings are always read as an integer vector starting at 0)
groupsizes = zeros(Int16, maximum(groupvect) + 1)
for ii in groupvect
    groupsizes[ii+1] += 1
end
## assemble empirical distribution
empirical = [zeros(Int16, maximum(groupvect) + 1) for xx in 1:numclasses]
for ii in unique(classvect)
    for jj in unique(groupvect)
        empirical[findfirst(xx -> xx == ii, unique(classvect))][jj+1] =
            sum(groupvect[classvect .== ii] .== jj)
    end 
end
## Sanity checks
if sum(groupsizes) != sum(sppcounts) != numspp
    error("Group sizes do not sum to total number of species")
end
if size(hcat(empirical...)) != (length(groupsizes), length(sppcounts))
    error("Empirical observation does not conform with number of groups, classes")
end
if sum(empirical) != groupsizes
    error("Impossible empirical distribution given group sizes")
end
if [sum(xx) for xx in empirical] != sppcounts
    error("Impossible empirical distribution given species counts per class")
end
## Find pvalue of empirical observation
if numsamps > 0
    psi = Imbalance(empirical...)
    ## Repeatedly, uniformly sample potential distributions
    outfile = open(outdir * "$(rsplit(basename(datafile), ".", limit=2)[1])-$cc", "w")
    SampleDists(cc, psi, groupvect, classvect, numsamps, outfile)
    println(outdir * "$(rsplit(basename(datafile), ".", limit=2)[1])-$cc")
    close(outfile)
else
    presencevects = [GetPotentialPresences(groupsizes, pp) for pp in sppcounts]
    fullpresences = GetCompatiblePresences(groupsizes, presencevects...)
    imbatable = GetAllImbalances(fullpresences)
    psi = Imbalance(empirical...)
    prob = Emppvalue(empirical, imbatable)
    println("ψ = $psi, p(ψ) = $prob")
    outfile = open(outdir * "$(rsplit(basename(datafile), ".", limit=2)[1])", "w")
    write(outfile, "$cc,$psi,$prob\n")
    close(outfile)
end
