"""
TimeTrees
---------
Basic rooted phylogenetic time tree manipulation module which includes Tree
and Node classes.  The Tree class can be initialized from a Newick string
and provides methods for producing visualizations using ASCII art.
"""
module TimeTrees

export Node, isRoot, addChild, edgeLength, Tree,
    getLeaves, getNodes

import Base.show

"""
The `Node` type is the basic element of any `TimeTree`. It
defines attributes which include links to a parent and
any children, the height of the node, and its label.
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
    if length(n.children)>0
        print(io, "(")
        for i in 1:length(n.children)
            if i > 1
                print(io, ",")
            end
            show(io, n.children[i])
        end
        print(io, ")")
    end
    print(io, string(n.label, ":", edgeLength(n)))
end

"""
Returns `true` if `n` is a root node.
"""
function isRoot(n::Node)
    return n.parent == n
end

"""
Returns `true` if `n` is a leaf node.
"""
function isLeaf(n::Node)
    return length(n.children) == 0
end

"""
Add node `c` as a child of node `n`.
"""
function addChild(n::Node, c::Node)
    c.parent = n
    push!(n.children, c)
end

"""
Returns the length of the edge above node `n`.  Always
returns 0 if `n` is a root node.
"""
function edgeLength(n::Node)
    if isRoot(n)
        return 0.0
    else
        return n.parent.height - n.height
    end
end

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

function getNodes(n::Node)
    res = [n]
    for c in n.children
        res = [res; getNodes(c)]
    end

    return res
end

"""
A phylogenetic tree.
"""
type Tree
    root::Node
end

function show(io::IO, t::Tree)
    show(io, t.root)
    print(io, ";")
end

"""
Construct a `Tree` from a Newick string.
"""
function Tree(newick::AbstractString)

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
       return Tree(ruleN())
    end

    function ruleN()
        node = Node()
        if matchToken("open_paren") != nothing
            while true
                addChild(node, ruleN())
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
Produce an ASCII representation of a given Tree object.
"""
function plotASCII(tree::Tree, width = 70, labelLeaves = true)

end

end
