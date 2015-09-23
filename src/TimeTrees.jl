module TimeTrees

export Node, isRoot, addChild, edgeLength, Tree

import Base.show

type Node
    parent::Node
    children::Array{Node, 1}
    height::FloatingPoint
    label::String

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

function isRoot(n::Node)
    return n.parent == n
end

function addChild(n::Node, c::Node)
    c.parent = n
    push!(n.children, c)
end

function edgeLength(n::Node)
    if isRoot(n)
        return 0.0
    else
        return n.parent.height - n.height
    end
end

type Tree
    root::Node
end

function show(io::IO, t::Tree)
    show(io, t.root)
    print(io, ";")
end

# Assemble tree from Newick string
function Tree(newick::String)

    ### Parser ###

    patterns = {
        "open_paren" => r"^\(",
        "close_paren" => r"^\)",
        "colon" => r"^:",
        "comma" => r"^,",
        "number" => r"^\d+(\.\d*)?([eE]-?\d+)?",
        "string" => r"^((\w+)|(\"[^\"]*\")|(^'[^']*'))"
    }

    function matchToken(token::String; mustMatch = false)

        # Skip whitespace
        while isblank(newick[i]) && i<=length(newick)
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

        if beginswith(res, "\"") || beginswith(res, "'")
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

    function getTime(node::Node, currentTime::FloatingPoint)
        currentTime += node.height

        maxTime = currentTime
        for child in node.children
            maxTime = max(maxTime, getTime(child, currentTime))
        end

        return maxTime
    end

    function branchLengthToHeight(node::Node, treeHeight::FloatingPoint, currentTime::FloatingPoint)
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

end
