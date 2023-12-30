module Misc

using DelimitedFiles
using Plots

function plot_error(
    parent_dir,
    dirs,
    f_analytical::Function,
    nbuf,
    lev,
    tmax;
    f_error::Function = (err, i) -> err
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
                plt = plot!(x_in, f_error(error, i))
            end
        end
    end
    return (plt)
end

end
