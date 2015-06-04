#meta-analysis 

#simulation of data
# 0: reference marker; 1: alternate marker
mark <- setClass("mark", representation(type = "numeric"))
ref <- new("mark", type = 0)
alt <- new("mark", type = 1)

marker <- setClass("marker", representation(i = "mark", j = "mark"))


marker(i = mark(0), j = mark(0))



genome <- setClass("genome", slots = "matrix", )