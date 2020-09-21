@enum FastMEMethod begin
    _BalME
    _OLSME
    _BIONJ
    _NJ
    _UNJ
end

@enum FastMEModel begin
    _PDis
    _RYSym
    _RY
    _JC69
    _K2P
    _F81
    _F84
    _TN93
    _LogDet
end

HumanizeFastMEMethod = Dict(
    _BalME => "B",
    _OLSME => "O",
    _BIONJ => "I",
    _NJ => "N",
    _UNJ => "U",
)

HumanizeFastMEModel = Dict(
    _PDis => "P",
    _RYSym => "Y",
    _RY => "R",
    _JC69 => "J",
    _K2P => "K",
    _F81 => "1",
    _F84 => "4",
    _TN93 => "T",
    _LogDet => "L"
)

function fastme(method :: FastMEMethod, model :: FastMEModel, input :: AbstractString)
    m = HumanizeFastMEMethod[method]
    d = HumanizeFastMEModel[model]
    cmd = `fastme --dna=$d -m $m -i "$input" -o outtree -O mat.txt`
    run(cmd)
    read("outtree", String)
end