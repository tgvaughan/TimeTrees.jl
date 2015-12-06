using TimeTrees
using Base.Test

# Test newick parsing:
t = TimeTree("((A:1,B:1):1,C:2):0;")
r = t.root

# Test getNodes and getLeaves
nodes = getNodes(t)
leaves = getLeaves(t)
l = leaves[1]

# Test basic methods
@test isRoot(r) == true
@test isRoot(l) == false
@test isLeaf(r) == false
@test isLeaf(l) == true

# Test tree sorting
@test getDecendentCount(r) == 5
s = getSorted(t)
@test s.root.children[1].label == "C"
