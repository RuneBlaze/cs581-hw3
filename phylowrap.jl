include("wrappers.jl")
using Glob
using Serialization
@enum DistanceType LogDet JC PDis
@enum Model GTR JCM

abstract type TreeMethod end
struct NJ <: TreeMethod
    distance :: DistanceType
end

struct BME <: TreeMethod
    distance :: DistanceType
end

struct UPGMA <: TreeMethod
    distance :: DistanceType
end

struct RAxML <: TreeMethod
    model :: Model
end

struct IQTree <: TreeMethod
    model :: Model
end

struct FastTree2 <: TreeMethod
    model :: Model
end

struct MPBoot <: TreeMethod
end

const TreeMethod2Speed = Dict(
    NJ => 0,
    BME => 1,
    UPGMA => 2,
    RAxML => 3,
    IQTree => 4,
    FastTree2 => 5,
)

function lessthan(x :: TreeMethod, y :: TreeMethod)
    lessthan(TreeMethod2Speed(x), TreeMethod2Speed(y))
end

const FastaPath = AbstractString
const TTPath = AbstractString

# buildtree blindly executes the command and returns the file handle
# to the output tree.
function buildtree(method :: TreeMethod, input :: FastaPath)
    error("no concrete method specified")
end

struct Dataset
    rootpath :: String
end

const PhyTree = String

struct ExpKey
    path :: String
    method :: TreeMethod
end

function lessthan(x :: ExpKey, y :: ExpKey)
    lessthan(x.method, y.method)
end

function allinputs(ds :: Dataset)
    glob("*/rose.aln.true.phy", ds.rootpath)
end

const Method2FastMEDist = Dict(LogDet => _LogDet, JC => _JC69, PDis => _PDis)

function buildtree(method :: NJ, input :: FastaPath)
    d = Method2FastMEDist[method.distance]
    fastme(_NJ, d, input)
end

function buildtree(method :: BME, input :: FastaPath)
    d = Method2FastMEDist[method.distance]
    fastme(_BalME, d, input)
end

function buildtree(method :: FastTree2, input :: FastaPath)
    cmd = method.model == GTR ? `FastTree -gtr -nt` : `FastTree -nt`
    run(pipeline(cmd, stdin=input, stdout="outtree"))
    read("outtree", String)
end

function buildtree(method :: IQTree, input :: FastaPath)
    specify = Dict(GTR => "GTR", JCM => "JC")[method.model]
    cmd = `iqtree -s $input -m $specify`
    run(cmd)
    read(input * ".treefile", String)
end

function buildtree(method :: UPGMA, input :: FastaPath)
    fasta = replace(input, ".phy" => ".fasta")
    cmd = `megacc -a UPGMA.mao -d $fasta -o outtree`
    run(cmd)
    read("outtree.nwk", String)
end

function buildtree(method :: MPBoot, input :: FastaPath)
    cmd = `./mpboot -s $input`
    run(cmd)
    read(input * ".treefile", String)
end