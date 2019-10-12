include("classifiers.jl")
clf3 = load_classification("/DATA/classification3.clf")
fam = clf3[16][1]
pi = Pairing(fst)
as_perm(fam, pi, 1)

