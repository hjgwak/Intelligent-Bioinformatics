import sys

__author__ = "hjgwak"
__version__ = "1.0.1"

def printUSAGE() :
	print """
###########################################################
# python convertARFF [options] [input file]
# options :
#     -k	#	:	use k-mer to convert (1~4, Default 1)
#     -sum	0,1	:	Set 1 : sum / Set 0 : bitwise (Default 0)
#     -weight	0,1	:	Set 1 : give weight for amino-acid by order of octamer (12344321)
#                       !! don't use this option when set 0 for -sum (Default 0)
#     -abs	0,1	:	Set 1 : using abs value of difference between left and right substring
#                   !! don't use this option when set 0 for -sum (Default 0)
###########################################################
"""

def printDescription(Sum) :
	if Sum :
		print """@relation octamer_cleave

@attribute Hydrophobic numeric
@attribute Positive numeric
@attribute Negative numeric
@attribute Polar numeric
@attribute Charged numeric
@attribute Small numeric
@attribute Tiny numeric
@attribute Aliphatic numeric
@attribute Aromatic numeric
@attribute Proline numeric
@attribute cleave {cleave, not_cleave}

@data
"""
	else :
		print """@relation octamer_cleave

@attribute Hydrophobic {Same, Different}
@attribute Positive {Same, Different}
@attribute Negative {Same, Different}
@attribute Polar {Same, Different}
@attribute Charged {Same, Different}
@attribute Small {Same, Different}
@attribute Tiny {Same, Different}
@attribute Aliphatic {Same, Different}
@attribute Aromatic {Same, Different}
@attribute Proline {Same, Different}
@attribute cleave {cleave, not_cleave}

@data
"""

def setDefault(Sum) :
	res = 0b0000000000
	if Sum :
		res = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	return res

def convertBit2List(bit, n) :
	divider = 2**(n-1)
	res_list = []
	for i in range(n) :
		div = bool(bit & divider)
		if div :
			res_list.append(1)
		else :
			res_list.append(0)
		divider /= 2

	return res_list

def mergeVector(vec, ten, Sum, weight) :
	res = []
	if Sum :
		ten_list = convertBit2List(ten, 10)
		for i in range(10) :
			res.append(vec[i] + weight * ten_list[i])
	else :
		res = vec | ten

	return res

def diff(x, y, Abs) :
	res = x - y
	if Abs :
		res = abs(res)
	return res

def SameDiff(boolean) :
	if boolean :
		return "Same"
	else :
		return "Different"

# main

TVD_table = {    # Hydrophobic Positive Negative Polar Charged Small Tiny Aliphatic Aromatic Proline
	"A" : 0x218, # 1           0        0        0     0       1     1    0         0        0
	"R" : 0x160, # 0           1        0        1     1       0     0    0         0        0
	"N" : 0x050, # 0           0        0        1     0       1     0    0         0        0
	"D" : 0x0F0, # 0           0        1        1     1       1     0    0         0        0
	"C" : 0x210, # 1           0        0        0     0       1     0    0         0        0
	"Q" : 0x040, # 0           0        0        1     0       0     0    0         0        0
	"E" : 0x0E0, # 0           0        1        1     1       0     0    0         0        0
	"G" : 0x218, # 1           0        0        0     0       1     1    0         0        0
	"H" : 0x362, # 1           1        0        1     1       0     0    0         1        0
	"I" : 0x204, # 1           0        0        0     0       0     0    1         0        0
	"L" : 0x204, # 1           0        0        0     0       0     0    1         0        0
	"K" : 0x360, # 1           1        0        1     1       0     0    0         0        0
	"M" : 0x200, # 1           0        0        0     0       0     0    0         0        0
	"F" : 0x202, # 1           0        0        0     0       0     0    0         1        0
	"P" : 0x011, # 0           0        0        0     0       1     0    0         0        1
	"S" : 0x058, # 0           0        0        1     0       1     1    0         0        0
	"T" : 0x250, # 1           0        0        1     0       1     0    0         0        0
	"W" : 0x242, # 1           0        0        1     0       0     0    0         1        0
	"Y" : 0x242, # 1           0        0        1     0       0     0    0         1        0
	"V" : 0x214  # 1           0        0        0     0       1     0    1         0        0
}

if '-h' in sys.argv or '-help' in sys.argv :
	printUSAGE()
	exit(-1)

# open data file
data = open(sys.argv[-1], 'r')

# getting options from command
k = 1
if '-k' in sys.argv :
	k = int(sys.argv[sys.argv.index('-k') + 1])
Sum = False
if '-sum' in sys.argv :
	Sum = bool(int(sys.argv[sys.argv.index('-sum') + 1]))
weight = False
if '-weight' in sys.argv :
	weight = bool(int(sys.argv[sys.argv.index('-weight') + 1]))
Abs = False
if '-abs' in sys.argv :
	Abs = bool(int(sys.argv[sys.argv.index('-abs') + 1]))

if ((not Sum) and weight) or ((not Sum) and Abs) :
	print "-weight or -abs option only use when set 1 for -sum\n"
	printUSAGE()
	exit(-1)

# print description
printDescription(Sum)

# convert data
for line in data.readlines() :
	line = line.rstrip('\r\n')
	(octamer, cleave) = line.split(',')
	octamer = octamer.upper()
	pos_weight = 1
	# calc left side substring
	left = setDefault(Sum)
	for i in range(k) :
		tenbit = TVD_table[octamer[3 - i]]
		if weight :
			pos_weight = 4 - i
		left = mergeVector(left, tenbit, Sum, pos_weight)
	# calc right side substring
	right = setDefault(Sum)
	for i in range(k) :
		tenbit = TVD_table[octamer[4 + i]]
		if weight :
			pos_weight = 4 - i
		right = mergeVector(right, tenbit, Sum, pos_weight)
	# convert bit 2 string
	if not Sum :
		left = convertBit2List(left, 10)
		right = convertBit2List(right, 10)
	# compare between left and right side
	res = []
	if Sum :
		for i in range(10) :
			res.append(str(diff(left[i], right[i], Abs)))
	else :
		for i in range(10) :
			res.append(SameDiff(not(left[i] ^ right[i])))
	# add cleave information
	if cleave == '1' :
		res.append("cleave")
	else :
		res.append("not_cleave")
	# report resut
	print ','.join(res)

data.close()
