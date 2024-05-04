module Misc

using DelimitedFiles

function load_data(dirs, nbuf, lev, tmax; column=11, thorn="testsubcyclingmc2")
    dats = []
    for dir in dirs
        tsvs = filter(file -> endswith(file, ".x.tsv"), readdir(dir))
        tsvs = filter(file -> startswith(file, thorn*"-error"), tsvs)
        x_in = []
        f_in = []
        for fname in tsvs
            data = readdlm(joinpath(dir, fname), Float64, comments = true)
            if isapprox(data[1, 2], tmax; rtol = 1e-12)
                x = data[findall(x -> x == lev, data[:, 4]), 8]
                f = data[findall(x -> x == lev, data[:, 4]), column]
                x_in = x[1+nbuf:length(x)-nbuf]
                f_in = f[1+nbuf:length(x)-nbuf]
            end
        end
        push!(dats, [x_in, f_in])
    end
    return dats
end

end
