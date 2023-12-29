module Misc

using DelimitedFiles
using Plots

function plot_error(
    f_analytical::Function,
    dirs;
    nbuf = 2,
    lev = 1,
    tmax = 2.0,
    cord = 0,
    parent_dir = expanduser("~/docker-workspace/simulations/"),
)
    plt = plot()
    for i = 1:length(dirs)
        fulldir = parent_dir * dirs[i]
        tsvs = filter(file -> endswith(file, ".tsv"), readdir(fulldir))
        for fname in tsvs
            data = readdlm(joinpath(fulldir, fname), Float64, comments = true)
            if isapprox(data[1, 2], tmax; rtol = 1e-12)
                x = data[findall(x -> x == lev, data[:, 3]), 5]
                f = data[findall(x -> x == lev, data[:, 3]), 6]
                x_in = x[1+nbuf:length(x)-nbuf]
                f_in = f[1+nbuf:length(x)-nbuf]
                error = (f_in - f_analytical.(tmax, x_in))
                plt = plot!(x_in, error * (2^cord)^(i-1))
            end
        end
    end
    return (plt)
end

end
