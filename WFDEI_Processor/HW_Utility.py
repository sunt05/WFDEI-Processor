import numpy as np


def memoize(fn):
    """returns a memoized version of any function that can be called
    with the same list of arguments.
    Usage: foo = memoize(foo)
    http://stackoverflow.com/a/17268784/920789"""

    def handle_item(x):
        if isinstance(x, dict):
            return make_tuple(sorted(x.items()))
        elif hasattr(x, '__iter__'):
            return make_tuple(x)
        else:
            return x

    def make_tuple(L):
        return tuple(handle_item(x) for x in L)

    def foo(*args, **kwargs):
        items_cache = make_tuple(sorted(kwargs.items()))
        args_cache = make_tuple(args)
        if (args_cache, items_cache) not in foo.past_calls:
            foo.past_calls[(args_cache, items_cache)] = fn(*args, **kwargs)
        return foo.past_calls[(args_cache, items_cache)]
    foo.past_calls = {}
    foo.__name__ = 'memoized_' + fn.__name__
    return foo



# split sequential list into sublists if sequence is broken
def splitSeq2D(seq):
    sub = [seq[0]]
    res = []
    i = 1
    while i < len(seq):
        if int(seq[i, 0]) - int(seq[i - 1, 0]) == 1:
            sub.append(seq[i])
        else:
            res.append(sub)
            sub = [seq[i]]
        i += 1
    res.append(sub)
    return res


def group1st(data):
    keys = list(set([x[0] for x in data]))
    return [[y for y in data if y[0] == x] for x in keys]

# split sequential list into sublists if not overlapped


def splitOlp2D(xseq):
    seq = sorted(xseq)
    sub = [seq[0]]
    res = []
    i = 1
    while i < len(seq):
        if seq[i - 1][0] < seq[i][0] < seq[i - 1][1] or seq[i - 1][0] < seq[i][1] < seq[i - 1][1]:
            sub.append(seq[i])
        else:
            res.append(sub)
            sub = [seq[i]]
        i += 1
    res.append(sub)
    return res

# when overlapped, pick out the longest one
# if several longest ones exist, choose the first occurence


def selectOlpLong(seq):
    resOlp = splitOlp2D(seq)
    resOlpMax = [max([y[1] - y[0] for y in x]) for x in resOlp]
    resOlpLong = [[x1 for x1 in x if x1[1] - x1[0] == y][0]
                  for [x, y] in zip(resOlp, resOlpMax)]
    return resOlpLong


def HWFinder(rawdata, T1, T2):
    seq = np.arange(len(rawdata))

    temp = rawdata[:, 1]
    data = np.transpose([seq, temp])

    # split rawdata into sublists based on T2
    sub = [x for x in data if x[1] > T2]
    subl = [np.array(x) for x in splitSeq2D(np.array(sub)) if len(x) >= 3]

    # search HWs in sublists
    dictHW = {}

    def HWSearch(data, T1, T2):
        # shortest scenario:
        if np.sort(data[:, 1])[-3] <= T1:
            return False
        # longer scenarios:
        else:
            if np.mean(data[:, 1]) > T1:
                dictHW[str(data[0, 0]) + str(data[-1, 0])] = data
            else:
                return HWSearch(data[:-1], T1, T2), HWSearch(data[1:], T1, T2)

    HWSearch = memoize(HWSearch)

    [HWSearch(x, T1, T2) for x in subl]

    indList = [x[[0, -1], 0] for x in dictHW.values()]

    # combine results
    indListGrp = [np.array(x) for x in group1st(indList)]

    # get the longest sublist
    indListMax = [[x[0, 0], max(x[:, 1])] for x in indListGrp]

    # when overlaped, select the longer ones
    indListLong = selectOlpLong(indListMax)

    return np.sort(indListLong)
