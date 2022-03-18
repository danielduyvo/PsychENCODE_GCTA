using Pkg
Pkg.activate(homedir())

using CairoMakie, CSV, DataFrames, Statistics

tri1_path =
"data/Trimester_1_EUR_25_HPC_gene_filtered_1Mbase_window/output/results/all_variance.txt"
tri2_path = 
"data/Trimester_2_EUR_15_HPC_gene_filtered_1Mbase_window/output/results/all_variance.txt"
tri1_sqtl_path =
"data/tri1/hg19_sQTL_15HPC/results/prelim_1Mbase_window.csv"
adult_path =
"data/EUR_SPC_gene_HRC_1Mbase_window/output/results/all_variance.txt"

readh2 = function(path::String)
    return CSV.read(path, DataFrame,
                    delim = ' ',
                    header = 1,
                    types = [String, Float64, Float64, Float64, Float64, Float64, 
                             Float64, Float64, Float64, Float64, Float64, Float64, 
                             Float64, Float64, Float64, Float64, Int64 ])
end
tri1 = readh2(tri1_path)
tri1sqtl = readh2(tri1_sqtl_path)
tri2 = readh2(tri2_path)
adult = readh2(adult_path)

violinsplit = function(gp::GridPosition, arr1::Vector{<:Number}, arr2::Vector{<:Number})
    xs = fill(0, length(arr1) + length(arr2))
    ys = vcat(arr1, arr2)
    side = vcat(fill(:left, length(arr1)), fill(:right, length(arr2)))
    color = @. ifelse(side == :left, :teal, :orange)
    violin(gp, xs, ys, side = side, color = color)
end

overlapdensity2 = function(ax::Axis,
                        arr1::AbstractArray{<:Number,1},
                        arr2::AbstractArray{<:Number,1})
    fillcolor1 = (:teal, 0.2)
    fillcolor2 = (:orange, 0.2)
    density!(ax, arr1, color = fillcolor1,
            strokecolor = :teal, strokewidth = 3, strokearound = true)
    density!(ax, arr2, color = fillcolor2,
            strokecolor = :orange, strokewidth = 3, strokearound = true)
end

overlapdensity3 = function(ax::Axis,
                        arr1::AbstractArray{<:Number,1},
                        arr2::AbstractArray{<:Number,1},
                        arr3::AbstractArray{<:Number,1})
    fillcolor1 = (:teal, 0.2)
    fillcolor2 = (:orange, 0.2)
    fillcolor3 = (:purple, 0.2)
    density!(ax, arr1, color = fillcolor1,
            strokecolor = :teal, strokewidth = 3, strokearound = true)
    density!(ax, arr2, color = fillcolor2,
            strokecolor = :orange, strokewidth = 3, strokearound = true)
    density!(ax, arr3, color = fillcolor3,
            strokecolor = :purple, strokewidth = 3, strokearound = true)
end


boxplotsplit2 = function(ax::Axis,
                        arr1::AbstractArray{<:Number,1},
                        arr2::AbstractArray{<:Number,1})
    fillcolor1 = (:teal, 0.2)
    fillcolor2 = (:orange, 0.2)
    xs = fill(0, length(arr1) + length(arr2))
    ys = vcat(arr1, arr2)
    dodge = vcat(fill(1, length(arr1)), fill(2, length(arr2)))
    color = Tuple[]
    for el in dodge
        if (el == 1)
            push!(color, fillcolor1)
        else
            push!(color, fillcolor2)
        end
    end
    boxplot!(xs, ys, dodge = dodge, color = color, orientation = :horizontal)
end

boxplotsplit3 = function(ax::Axis,
                        arr1::AbstractArray{<:Number,1},
                        arr2::AbstractArray{<:Number,1},
                        arr3::AbstractArray{<:Number,1})
    fillcolor1 = (:teal, 0.2)
    fillcolor2 = (:orange, 0.2)
    fillcolor3 = (:purple, 0.2)
    xs = fill(0, length(arr1) + length(arr2) + length(arr3))
    ys = vcat(arr1, arr2, arr3)
    dodge = vcat(fill(1, length(arr1)),
                 fill(2, length(arr2)),
                 fill(3, length(arr3)))
    color = Tuple[]
    for el in dodge
        if (el == 1)
            push!(color, fillcolor1)
        elseif (el == 2)
            push!(color, fillcolor2)
        else
            push!(color, fillcolor3)
        end
    end
    boxplot!(ax, xs, ys, dodge = dodge, color = color, orientation = :horizontal)
end

# Cis h2 for adults and both trimesters
f = Figure(resolution = (1000, 1000))
axtop = Axis(f[1,1])
axbot = Axis(f[2,1])
linkxaxes!(axtop, axbot)
boxplotsplit3(axbot, tri1.Cis_h2, tri2.Cis_h2, adult.Cis_h2)
overlapdensity3(axtop, tri1.Cis_h2, tri2.Cis_h2, adult.Cis_h2)
Label(f[1, 1, Top()], "Tri 1 vs. Tri 2 vs. Adult Cis-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
axbot.xlabel = "Cis-h\u00b2"
axbot.yreversed = true
Legend(f[1:2,2],
       [
        [LineElement(color = :teal, linestyle = nothing)],
        [LineElement(color = :orange, linestyle = nothing)],
        [LineElement(color = :purple, linestyle = nothing)]
       ],
       ["Tri 1", "Tri 2", "Adult"])
save("tri1_vs_tri2_vs_adult_cis_h2.pdf", f)

# Cis h2 for both trimesters
f = Figure(resolution = (1000, 1000))
axtop = Axis(f[1,1])
axbot = Axis(f[2,1])
linkxaxes!(axtop, axbot)
boxplotsplit2(axbot, tri1.Cis_h2, tri2.Cis_h2)
overlapdensity2(axtop, tri1.Cis_h2, tri2.Cis_h2)
Label(f[1, 1, Top()], "Tri 1 vs. Tri 2 Cis-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
axbot.xlabel = "Cis-h\u00b2"
axbot.yreversed = true
Legend(f[1:2,2],
       [
        [LineElement(color = :teal, linestyle = nothing)],
        [LineElement(color = :orange, linestyle = nothing)]
       ],
       ["Tri 1", "Tri 2"])
save("tri1_vs_tri2_cis_h2.pdf", f)

# Cis h2 for RNA-seq and sQTL for Trimester 1
f = Figure(resolution = (1000, 1000))
axtop = Axis(f[1,1])
axbot = Axis(f[2,1])
linkxaxes!(axtop, axbot)
boxplotsplit2(axbot, tri1.Cis_h2, tri1sqtl.Cis_h2)
overlapdensity2(axtop, tri1.Cis_h2, tri1sqtl.Cis_h2)
Label(f[1, 1, Top()], "Tri 1 RNA-seq vs. sQTL Cis-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
axbot.xlabel = "Cis-h\u00b2"
axbot.yreversed = true
Legend(f[1:2,2],
       [
        [LineElement(color = :teal, linestyle = nothing)],
        [LineElement(color = :orange, linestyle = nothing)]
       ],
       ["RNA-seq", "sQTL"])
save("tri1_rnaseq_vs_sqtl_cis_h2.pdf", f)

# Trans h2 for RNA-seq and sQTL for Trimester 1
f = Figure(resolution = (1000, 1000))
axtop = Axis(f[1,1])
axbot = Axis(f[2,1])
linkxaxes!(axtop, axbot)
boxplotsplit2(axbot, tri1.Trans_h2, tri1sqtl.Trans_h2)
overlapdensity2(axtop, tri1.Trans_h2, tri1sqtl.Trans_h2)
Label(f[1, 1, Top()], "Tri 1 RNA-seq vs. sQTL Trans-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
axbot.xlabel = "Trans-h\u00b2"
axbot.yreversed = true
Legend(f[1:2,2],
       [
        [LineElement(color = :teal, linestyle = nothing)],
        [LineElement(color = :orange, linestyle = nothing)]
       ],
       ["RNA-seq", "sQTL"])
save("tri1_rnaseq_vs_sqtl_trans_h2.pdf", f)

# Cis h2 and SE for adults and both trimesters
f = Figure(resolution = (1000, 1000))
axtop = Axis(f[1,1])
axbot = Axis(f[2,1])
axbotright = Axis(f[2,2])
linkxaxes!(axtop, axbot)
linkyaxes!(axbot, axbotright)
leg = Legend(f[1,2],
       [
        [LineElement(color = :teal, linestyle = nothing)],
        [LineElement(color = :orange, linestyle = nothing)],
        [LineElement(color = :purple, linestyle = nothing)]
       ],
       ["Tri 1", "Tri 2", "Adult"])
leg.tellheight = false
leg.tellwidth = false
boxplotsplit3(axbot, tri1.Cis_h2, tri2.Cis_h2, adult.Cis_h2)
boxplotsplit3(axbotright, tri1.SE_Cis_h2, tri2.SE_Cis_h2, adult.SE_Cis_h2)
overlapdensity3(axtop, tri1.Cis_h2, tri2.Cis_h2, adult.Cis_h2)
Label(f[1, 1, Top()], "Tri 1 vs. Tri 2 vs. Adult Cis-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
Label(f[2, 2, Top()], "Tri 1 vs. Tri 2 vs. Adult\nStandard Error of Cis-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
axbot.xlabel = "Cis-h\u00b2"
axbotright.xlabel = "SE of Cis-h\u00b2"
axbot.yreversed = true
axbotright.yreversed = true
save("tri1_vs_tri2_vs_adult_cis_h2.pdf", f)

# Trans h2 and SE for adults and both trimesters
f = Figure(resolution = (1000, 1000))
axtop = Axis(f[1,1])
axbot = Axis(f[2,1])
axbotright = Axis(f[2,2])
linkxaxes!(axtop, axbot)
linkyaxes!(axbot, axbotright)
leg = Legend(f[1,2],
       [
        [LineElement(color = :teal, linestyle = nothing)],
        [LineElement(color = :orange, linestyle = nothing)],
        [LineElement(color = :purple, linestyle = nothing)]
       ],
       ["Tri 1", "Tri 2", "Adult"])
leg.tellheight = false
leg.tellwidth = false
boxplotsplit3(axbot, tri1.Trans_h2, tri2.Trans_h2, adult.Trans_h2)
boxplotsplit3(axbotright, tri1.SE_Trans_h2, tri2.SE_Trans_h2, adult.SE_Trans_h2)
overlapdensity3(axtop, tri1.Trans_h2, tri2.Trans_h2, adult.Trans_h2)
Label(f[1, 1, Top()], "Tri 1 vs. Tri 2 vs. Adult Trans-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
Label(f[2, 2, Top()], "Tri 1 vs. Tri 2 vs. Adult\nStandard Error of Trans-h\u00b2", valign = :bottom,
    padding = (0, 0, 5, 0))
axbot.xlabel = "Trans-h\u00b2"
axbotright.xlabel = "SE of Trans-h\u00b2"
axbot.yreversed = true
axbotright.yreversed = true
save("tri1_vs_tri2_vs_adult_trans_h2.pdf", f)

temp = innerjoin(tri1, tri2, on = :ID, makeunique=true)
temp1heritable = temp.Cis_h2 - (temp.SE_Cis_h2 .* 1.96) .> 0 # length: 2065
temp2heritable = temp.Cis_h2_1 - (temp.SE_Cis_h2_1 .* 1.96) .> 0 # length: 1069
temp1and2heritable = temp1heritable .&& temp2heritable # length: 519
temp1r = sqrt(sum(temp1heritable)/pi) # 25.6
temp2r = sqrt(sum(temp2heritable)/pi) # 18.4
dist1to2 = 23.5 # stolen from some circle intersection calculator

# Venn diagram of heritable genes
f = Figure(resolution = (1000, 1000))
ax = Axis(f[1,1], aspect = DataAspect())
leg = Legend(f[1,2],
       [
        [LineElement(color = :teal, linestyle = nothing)],
        [LineElement(color = :orange, linestyle = nothing)]
       ],
       ["Tri 1", "Tri 2"])
poly!(ax, Circle(Point2f(-dist1to2/2, 0), temp1r), color = (:teal, 0.5), strokecolor = :teal, strokewidth = 1)
poly!(ax, Circle(Point2f(dist1to2/2, 0), temp2r), color = (:orange, 0.5), strokecolor = :orange, strokewidth = 1)
hidedecorations!(ax)
text!(ax, "$(sum(temp1heritable) - sum(temp1and2heritable))", position = (-dist1to2/2 - 10, 0), textsize = 25)
text!(ax, "$(sum(temp1and2heritable))", position = (0 + 5, 0), textsize = 25)
text!(ax, "$(sum(temp2heritable) - sum(temp1and2heritable))", position = (dist1to2/2 + 10, 0), textsize = 25)
Label(f[1, 1, Top()], "Heritable genes in Tri 1 and Tri 2", valign = :bottom,
    padding = (0, 0, 5, 0))
save("tri1_vs_tri2_venn.pdf", f)

