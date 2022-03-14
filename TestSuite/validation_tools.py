import os,sys, getopt
import numpy as np

import struct
import math

import colorcodes as cc
_c = cc.Colorcodes()

def compute_ulp(x):
    """Return the value of the least significant bit of a
    float x, such that the first float bigger than x is x+ulp(x).
    Then, given an expected result x and a tolerance of n ulps,
    the result y should be such that abs(y-x) <= n * ulp(x).
    The results from this function will only make sense on platforms
    where native doubles are represented in IEEE 754 binary64 format.

    From official test_math.py. Python 3.9 has now math.ulp
    """

    x = abs(float(x))
    if math.isnan(x) or math.isinf(x):
        return x

    # Find next float up from x.
    n = struct.unpack('<q', struct.pack('<d', x))[0]
    x_next = struct.unpack('<d', struct.pack('<q', n + 1))[0]
    if math.isinf(x_next):
        # Corner case: x was the largest finite float. Then it's
        # not an exact power of two, so we can take the difference
        # between x and the previous float.
        x_prev = struct.unpack('<d', struct.pack('<q', n - 1))[0]
        return x - x_prev
    else:
        return x_next - x


def processArgv(argv):
    ref_dir = None
    data_dir = None

    usage = 'validate_benchmark.py -d <data_dir> -r <ref_dir>'
    try:
        opts, args = getopt.getopt(argv,"hd:r:")
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-d"):
            data_dir = arg + "/"
        elif opt in ("-r"):
            ref_dir = arg + "/"

    return (ref_dir, data_dir)

def printResult(condition, msg, ntabs = 1):
    """Pretty print test results"""

    res = np.zeros(2, dtype='i8')
    if condition:
        status = (_c.green + b'passed' + _c.reset).decode()
    else:
        status = (_c.red + b'failed' + _c.reset).decode()
        res[1] += 1
    res[0] += 1
    tabs = "\t"*ntabs
    print(tabs + msg.ljust(50) + f': {status}')

    return res

def printSummary(results):
    """Pretty print validation test summary"""

    print("")
    if(results[1] == 0):
        print(f'All benchmark validation tests passed!')
    else:
        t = 'test'
        if results[1] > 1:
            t += 's'
        tc = _c.red + (f'{results[1]} benchmark validation {t} failed').encode()  + _c.reset
        msg = tc + (f' out of {results[0]}').encode()
        print(msg.decode())

def tableTest(fname, ref_dir, data_dir, tol = 11, usecols = None, max_rows = None, threshold = -1, percol = False):

    # Validate nusselt number
    extra = ''
    if usecols is not None:
        extra = f' usecols = {usecols}'
    if max_rows is not None:
        extra = f' max_rows = {max_rows}'
    if extra:
        extra = ' (' + extra + ' )'
    print(f'Validating {fname}{extra}')
    res = np.zeros(2, dtype='i8')

    cond = True

    if cond:
        cond = os.path.exists(ref_dir + fname)
        res += printResult(cond, 'Checking if reference exists')

    if cond:
        cond = os.path.exists(data_dir + fname)
        res += printResult(cond, 'Checking if data exists')

    if cond:
        ref = np.genfromtxt(ref_dir + fname, usecols=usecols, max_rows = max_rows)
        data = np.genfromtxt(data_dir + fname, usecols=usecols, max_rows = max_rows)
        cond = (ref.shape == data.shape)
        res += printResult(cond, f'Checking file size match (shape: {ref.shape})')
    
    if cond:
        max_ulp = 0
        if percol:
            col_max = np.max(np.abs(ref), axis = 0)
            for idx, r in np.ndenumerate(ref):
                if r > threshold:
                    d = data[idx]
                    diff = np.abs(r-d)
                    ulp = diff/(compute_ulp(col_max[idx[1]]))
                    if ulp > max_ulp:
                        max_ulp = ulp
                    if ulp > tol:
                        print((r.item(), d.item(), diff, ulp))
        else:
            for r, d in np.nditer([ref, data]):
                if r > threshold:
                    diff = np.abs(r-d)
                    ulp = diff/(compute_ulp(r))
                    if ulp > max_ulp:
                        max_ulp = ulp
                    if ulp > tol:
                        print((r.item(), d.item(), diff, ulp))

        cond = (max_ulp < tol)
        if max_ulp > 1e3:
            msg = f'Checking error tolerance (max ulp: {max_ulp:.3e})'
        else:
            msg = f'Checking error tolerance (max ulp: {max_ulp:.0f})'
        res += printResult(cond, msg)

    return res
