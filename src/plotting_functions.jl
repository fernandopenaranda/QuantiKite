function dos_plot(filename)
    # Read data from a .dat file
    data = readdlm(filename)
    # Extract columns for x and y
    x = data[:, 1]
    y = data[:, 2]
    f = Figure()
    ax = Axis(f[1,1], xlabel = "E (eV)", ylabel = "DOS (a.u.)")
    # Create a scatter plot
    lines!(ax, x, y)
    return f
end

function linear_optical_conductivity_plot(filename; ylims = [-0.5,3])
    # Read data from a .dat file
    data = readdlm(filename)
    # Extract columns for x and y
    x = data[:, 1]
    y = data[:, 2]
    z = data[:,3]
    f = Figure()
    ax = Axis(f[1,1], xlabel = "ħω (eV)", ylabel = "σ_ab (ω) [e^2/ħ]")
    # Create a scatter plot
    lines!(ax, x, y, color = :purple)
    lines!(ax, x, z, color = :orange)
    ylims!(ax, ylims)
    terms = ["Re[σ_ab]",  "Im[σ_ab]"]
    ccolors = [:orange, :purple]
    elems = [[MarkerElement(color = col, marker=:circle, markersize = 15, strokecolor = :black)] for col in ccolors]
    axislegend(ax, elems, terms; position=:rt)
    return f
end