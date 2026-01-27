import Levenshtein as lv


def compareString(str1, str2, method="lv"):
    allowed = {'A', 'T', 'C', 'G'}

    # if the string is not null and it only contains the characters A, T, C, or G
    if (str1 and set(str1).issubset(allowed)) and (str2 and set(str2).issubset(allowed)):
        if method == "lv":
            return lv.distance(str1, str2)
    
    else:
        return None


def NoneToZero(num):
    if num is None:
        return 0
    else:
        return num
    

def NoneToNA(input):
    if input is None:
        return 'NA'
    else:
        return input


def compareGt(gt1, gt2):

    # pad genotypes to length 2 for compatibility (eg. for handling hemizygous regions)
    gt1 += [None] * (2 - len(gt1))
    gt2 += [None] * (2 - len(gt2))

    # get distances for both permutations of comparing alleles
    vert_dist = [compareString(gt1[0], gt2[0]), compareString(gt1[1], gt2[1])]
    cross_dist = [compareString(gt1[1], gt2[0]), compareString(gt1[0], gt2[1])]

    v_sum = sum(NoneToZero(dist) for dist in vert_dist)
    c_sum = sum(NoneToZero(dist) for dist in cross_dist)

    # assume the lesser distance is the correct comparison order
    best_sum, best_dist = (v_sum, vert_dist) if v_sum <= c_sum else (c_sum, cross_dist)

	# returns the sum of the differences between alleles, as well as the individual differences
    return best_sum, NoneToNA(best_dist[0]), NoneToNA(best_dist[1])

'''
def compareGCContent(str1, str2):
    str1_gc = [count for count in str1[]]
    str2_gc = [count for count in str2[]]

    return str1_gc - str2_gc
''' 


def stateCheck(rdr):
    return (not rdr.pause) and (not rdr.end_state)
