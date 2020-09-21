using DataStructures
#40.76.67.99
#20.55.3.62
using JLD2
using FileIO
include("phylowrap.jl")
const D1000M1 = Dataset("datasets/1000M1")
const D1000M4 = Dataset("datasets/1000M4")
store = OrderedDict{ExpKey, Union{Missing, Tuple{Float64, PhyTree}}}()

function adddataset!(s, ds :: Dataset, method :: TreeMethod)
    for i in allinputs(ds)
        s[ExpKey(i, method)] = missing
    end
end

for ds in [D1000M1, D1000M4]
    adddataset!(store, ds, UPGMA(PDis))
end

# ds = D1000M1
# adddataset!(store, ds, IQTree(GTR))

if isfile("store.jld2")
    old_store = load("store.jld2")["store"]

    # merge the old store onto the existing one
    for (k, v) in old_store
        if (!haskey(store, k)) || ismissing(store[k])
            # store[k] = old_store[k]
        end
    end
end

execute(expkey :: ExpKey) = buildtree(expkey.method, expkey.path)

for (k, v) in store
    if ismissing(v)
        @show (k, v)
        try
            r = @timed execute(k)
            store[k] = (r.time, r.value)
            @show (r.time, r.value)
            @save "store.jld2" store
        catch e
            store[k] = missing
            @show e
        end
    end
end