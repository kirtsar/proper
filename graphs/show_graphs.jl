using ImageMagick
using GraphPlot
using Compose
using LightGraphs
using Cairo
using Fontconfig
using Images
using ImageView

using Colors



function gshow(g :: DiGraph; labels = Vector{String}([]))
    n = nv(g)
    m = Int(log2(n))
    if isempty(labels)
        for i in  0 : n - 1
            push!(labels, bitstring(i)[end - m + 1 : end])
        end
    end
    nodescolor = [colorant"white" for i in 1 : n]
    n = length(readdir("pics")) + 1
    draw(PNG("pics/graph$n.png", 14cm, 14cm), 
    gplot(g, nodelabel = labels, nodefillc=nodescolor))
    img = load("pics/graph$n.png")
    imshow(img)
end


function gshow(adjmat :: Matrix)
    gshow(DiGraph(adjmat))
end


using TikzGraphs
function latexshow(adjmat :: Matrix, labels = Vector{String}([]))
    n = size(adjmat)[1]
    g = DiGraph(adjmat)
    if isempty(labels)
        for i in  0 : n - 1
            push!(labels, bitstring(i)[end - n + 1 : end])
        end
    end
    TikzGraphs.plot(g, labels = labels)
end