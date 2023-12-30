module Misc

using DelimitedFiles

function load_data(dirs, nbuf, lev, tmax; column=6)
    dats = []
    for dir in dirs
        tsvs = filter(file -> endswith(file, ".tsv"), readdir(dir))
        x_in = []
        f_in = []
        for fname in tsvs
            data = readdlm(joinpath(dir, fname), Float64, comments = true)
            if isapprox(data[1, 2], tmax; rtol = 1e-12)
                x = data[findall(x -> x == lev, data[:, 3]), 5]
                f = data[findall(x -> x == lev, data[:, 3]), column]
                x_in = x[1+nbuf:length(x)-nbuf]
                f_in = f[1+nbuf:length(x)-nbuf]
            end
        end
        push!(dats, [x_in, f_in])
    end
    return dats
end

end
