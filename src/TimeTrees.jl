"""
TimeTrees
---------

Basic rooted phylogenetic time tree manipulation module which includes TimeTree
    and Node classes.  The TimeTree class can be initialized from a Newick
    string and provides methods for producing visualizations using ASCII art.

TimeTrees provides thw following two types:

* `TimeTree`: a rooted phylogenetic time tree.
* `Node`: the basic building block of a phylogenetic tree.

Read the documentation for these types for further information.
"""
module TimeTrees

export Node, isRoot, isLeaf, addChild!, edgeLength,
    getDecendentCount, getSorted,
    TimeTree, getLeaves, getNodes, getNewick, plot

import Base.show

"""
The `Node` type is the basic element of any `TimeTree`. It
defines the following attributes:

* `parent`: the parent node. If the node is a root, then this
attribute refers to the current node object.

* `children`: the list of child nodes. If empty, this node is
a leaf.

* `height`: the height or age of the node (usually measured
from the most recent sample/leaf node.

* `label`: a (possibly empty) string labeling this node. Ofen
only used for leaves.

An empty node can be constructed using the `Node()` method.

Several methods on nodes are defined: `isRoot`, `isLeaf`, `addChild!`,
`edgeLength`, `getLeaves`, `getNodes`, `getSorted`, `getNewick`.  Read the
documentation for these methods for further information.
"""
type Node
    parent::Node
    children::Array{Node, 1}
    height::AbstractFloat
    label::AbstractString

    Node() = begin
        n = new()
        n.parent = n
        n.children = []
        n.height = 0.0
        n.label = ""
        return n
    end
end

function show(io::IO, n::Node)
    if isRoot(n)
        print(io, "Root")
    else
        print(io, "Non-root")
    end

    print(io, " node (age: $(n.height), children: $(length(n.children))")

    if length(n.label)>0
        print(io,", label: $(n.label))")
    else
        print(io,", no label)")
    end

end

"""
`getNewick(n::Node)` retrieves the Newick representation of the subtree below `n`.
"""
function getNewick(n::Node)
    res = ""
    if length(n.children)>0
        res = string(res, "(")
        for i in 1:length(n.children)
            if i > 1
                res = string(res, ",")
            end
            res = string(res, getNewick(n.children[i]))
        end
        res = string(res, ")")
    end
    res = string(res, n.label, ":", edgeLength(n))
end

"""
`isRoot(n::Node)` returns `true` if `n` is a root node.
"""
function isRoot(n::Node)
    return n.parent == n
end

"""
`isLeaf(n::Node)` returns `true` if `n` is a leaf node.
"""
function isLeaf(n::Node)
    return length(n.children) == 0
end

"""
`addChild!(n::Node, c::Node)` adds node `c` as a child of node `n`.
"""
function addChild!(n::Node, c::Node)
    c.parent = n
    push!(n.children, c)
end

"""
`edgeLength(n::Node)` returns the length of the edge above node `n`.  Always
returns 0 if `n` is a root node.
"""
function edgeLength(n::Node)
    if isRoot(n)
        return 0.0
    else
        return n.parent.height - n.height
    end
end

"""
`getLeaves(n::Node)` returns array of leaf nodes below `n`.
"""
function getLeaves(n::Node)
    if isLeaf(n)
        return [n]
    else
        res = []
        for c in n.children
            res = [res; getLeaves(c)]
        end

        return res
    end
end

"""
`getNodes(n::Node)` returns array of nodes below (and including) node `n`.
"""
function getNodes(n::Node)
    res = [n]
    for c in n.children
        res = [res; getNodes(c)]
    end

    return res
end

"""
`getDecendentCount(n::Node)` returns number of decendents below (and including) `n`.
"""
function getDecendentCount(n::Node)
    res = 1
    for c in n.children
        res += getDecendentCount(c)
    end

    return res
end

"""
`getSorted(n::Node, [rev=false])` produces a copy of the clade below `n`,
sorted so that children appear in order of the number of decendents they have.
Setting `rev=true` reverses the sort.
"""
function getSorted(n::Node; rev = false)

    copy = Node()
    copy.label = n.label
    copy.height = n.height

    for c in n.children
        addChild!(copy, getSorted(c))
    end

    perm = sortperm([getDecendentCount(c) for c in copy.children], rev=rev)
    copy.children = copy.children[perm]

    return copy
end

"""
The `TimeTree` type represents a rooted phylogenetic tree.  It defines a
single attribute `root::Node` which specifies the root node of the tree.

`TimeTree` objects can be constructed in two ways:

1. Using the default constructor `TimeTree(root::Node)` to create a tree
with the given root node, and

2. Using the constructor `TimeTree(newick::String)` which creates a
tree by parsing a Newick string given as an argument.  For instance,
`TimeTree("((A:1,B:1):1,C:2):0;")` constructs a tree with 3 leaves and
a total height of two.

There are several methods which act on trees: `getLeaves`,
`getNodes`, `getSorted`, `getNewick` and `plot`.
"""
type TimeTree
    root::Node
end

"""
`getLeaves(t::TimeTree)` returns array of leaf nodes `t` contains.
"""
function getLeaves(t::TimeTree)
    return getLeaves(t.root)
end

"""
`getNodes(t::TimeTree)` returns array of nodes that `t` contains.
"""
function getNodes(t::TimeTree)
    return getNodes(t.root)
end

"""
`getSorted(t::TimeTree, [rev=false])` returns a copy of `t` in which children are
sorted in order of the number of of their respective decendents.
"""
function getSorted(t::TimeTree; rev = false)
    return TimeTree(getSorted(t.root))
end

function show(io::IO, t::TimeTree)
    print(io, string("A phylogenetic tree with ",
        length(getLeaves(t)), " leaves (",
        length(getNodes(t)), " nodes in total)"))
end

"""
`getNewick(t::TimeTree)` retrieves the Newick representation of `t`.
"""
function getNewick(t::TimeTree)
    return string(getNewick(t.root), ";")
end

"""
Construct a `TimeTree` from a Newick string.
"""
function TimeTree(newick::AbstractString)

    ### Parser ###

    patterns = Dict{AbstractString,Regex}(
        "open_paren" => r"^\(",
        "close_paren" => r"^\)",
        "colon" => r"^:",
        "comma" => r"^,",
        "number" => r"^\d+(\.\d*)?([eE]-?\d+)?",
        "string" => r"^((\w+)|(\"[^\"]*\")|(^'[^']*'))"
    )

    function matchToken(token::AbstractString; mustMatch = false)

        # Skip whitespace
        while i<=length(newick) && (newick[i] == ' ' || newick[i] == '\t')
            i += 1
        end

        m = i<=length(newick) ? match(patterns[token], newick[i:end]) : nothing
        if m == nothing
            if mustMatch
                c = newick[i]
                error("Expected token $token at index $i but got $c instead.")
            else
                return nothing
            end
        end

        i += length(m.match)
        return m.match
    end

    function ruleT()
       return TimeTree(ruleN())
    end

    function ruleN()
        node = Node()
        if matchToken("open_paren") != nothing
            while true
                addChild!(node, ruleN())
                if matchToken("comma") == nothing
                    break
                end
            end
            matchToken("close_paren", mustMatch = true)
        end

        node.label = ruleL()
        node.height = ruleH()

        return node
    end

    function ruleL()
        res = matchToken("string")
        if res == nothing
            return ""
        end

        if startswith(res, "\"") || startswith(res, "'")
            return res[2:end-1]
        else
            return res
        end
    end

    function ruleH()
        res = matchToken("colon")
        if res == nothing
            return 0.0
        end

        return float(matchToken("number", mustMatch = true))
    end

    i = 1
    tree =  ruleT()

    ### Post-processing ###

    function getTime(node::Node, currentTime::AbstractFloat)
        currentTime += node.height

        maxTime = currentTime
        for child in node.children
            maxTime = max(maxTime, getTime(child, currentTime))
        end

        return maxTime
    end

    function branchLengthToHeight(node::Node, treeHeight::AbstractFloat, currentTime::AbstractFloat)
        currentTime += node.height
        node.height = treeHeight - currentTime

        for child in node.children
            branchLengthToHeight(child, treeHeight, currentTime)
        end
    end

    height = getTime(tree.root, 0.0)
    branchLengthToHeight(tree.root, height, 0.0)

    return tree
end

"""
`plot(t::TimeTree)` produces an ASCII representation of `t`.
"""
function plot(t::TimeTree; width = 70, labelLeaves = true, dots = true)

    leaves = getLeaves(t)
    nodes = getNodes(t)
    pos = zeros(length(nodes))

    grid = fill(' ', length(leaves), width)

    function computePos(n::Node)
        idx = findfirst(nodes, n)
        if isLeaf(n)
            pos[idx] = findfirst(leaves, n)
        else
            for c in n.children
                pos[idx] += computePos(c)
            end
            pos[idx] /= length(n.children)
        end
        return pos[idx]
    end
    computePos(t.root)

    function heightToIdx(h::Float64)
        return round(Int, h/t.root.height*(width-1) + 1)
    end

    # Edges
    for i in 1:length(nodes)
        n = nodes[i]

        if isRoot(n)
            continue
        else
            pi = findfirst(nodes, n.parent)
            x1 = heightToIdx(nodes[i].height)
            y1 = round(Int, pos[i])
            x2 = heightToIdx(nodes[pi].height)
            y2 = round(Int, pos[pi])
            ymin = min(y1, y2)
            ymax = max(y1, y2)

            grid[y1,x1:x2] = '-'
            grid[ymin:ymax, x2] = '|'
            grid[ymin,x2] = '/'
            grid[ymax,x2] = '\\'
        end
    end

    # Nodes
    for i in 1:length(nodes)
        x = heightToIdx(nodes[i].height)
        y = round(Int, pos[i])
        if isLeaf(nodes[i])
            grid[y, x] = '*'
            if x>1 && labelLeaves && dots
                grid[y,x-1:-1:1] = '⋅'
            end
        else
            grid[y, x] = '+'
        end
    end

    # Display
    for i in 1:length(leaves)
        if labelLeaves
            println(string(join(reverse(vec(grid[i,:]))), " ", leaves[i].label))
        else
            println(string(join(reverse(vec(grid[i,:])))))
        end
    end

end

end
