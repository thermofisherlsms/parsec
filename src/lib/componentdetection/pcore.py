from lib.componentdetection.cdl import *
from lib.componentdetection.pdata import *

def kv_scan(s):
    """
    Tuple (file name , scan) where 0 is int for msOrder (e.g ms1 or ms2 or msN)
    """
    scan = get_pscan(s)
    return (s.FileName, scan)



def map_peak_mf(t,w):
    p =  detect_peaks_mf(t[1],w)
    return (t,p)


def ping_detect(x,y):
    return -1


def detect_peaks(t,sn):
    # get trace s/n
    stn = max([i.Intensity for i in t.Trace]) / min([i.Intensity for i in t.Trace])
    if stn > sn:
        #pad the trace with zeros to detect 1 or 2 point peaks
        #if len(t.Trace) < 3:
        pad = MZPlus()
        pad.FileID = t.Trace[0].FileID
        t.Trace.insert(0, pad)
        t.Trace.insert(len(t.Trace) + 1,pad)
        peaks = detect_peaks_w_derivatives(t.Trace)
        if peaks is not None:
                return peaks
        else:
            return []
    else:
        return []

def map_peak(p,t):
    #rounding up to group peaks
    rt = round(p.RT, 1)
    rt = math.ceil(rt * 10) / 10
    upper = str(round(rt + t, 1))
    lower = str(round(rt - t, 1))
    k = p.FileID + "_" + lower + "_" + upper
    return (k,[p])


def map_feature(f):
    """
    Map features where the keys are the A0 and RT so we can group and
    do alignment (RT offsets relative to reference file
    """
    return (f.FileID,f)



def bin_masses(scansPerFile):

    fileId = scansPerFile[0]
    scanList = scansPerFile[1]

    massBins = {}

    for scan in scanList:
        BinMassesByFile(fileId, scan, massBins)

    # create a list of key value tuples for spark
    flattenedList = []
    for k, masses in massBins.items():
            flattenedList.append((k, masses))

    return flattenedList

def create_mass_traces(fileMassBin):

    traceList = []
    binKeyItems = fileMassBin[0].split(",")
    fileId = binKeyItems[0]
    binNumber = int(binKeyItems[1])
    binMassList = fileMassBin[1]

    traces = build_mass_traces(binNumber, binMassList)
    # traceList.append(traces)
    for j in traces:
        traceList.append((fileId, j))

    return traceList

def BinMassesByFile(fileId, scan, binMasses):

    for i in range(len(scan.Masses)):
        binNumber = round(scan.Masses[i])
        binKey = "{},{}".format(fileId, binNumber)

        if(binKey in binMasses):
            massList = binMasses[binKey]
        else:
            massList = []
            binMasses[binKey] = massList

        mz = MZPlus()
        #temporary_hack
        #mz.FileID = id
        mz.FileID = fileId
        mz.Mass = scan.Masses[i]
        mz.Intensity = scan.Intensities[i]
        mz.RetentionTime = scan.RetentionTime
        mz.RT = scan.RetentionTime
        mz.ScanNumber = scan.ScanNumber
        massList.append(mz)

def flatten_features(f):
    for p in f:
        return f[1]

def featureA0_to_formula(f):
    # TODO we must deconvolute with charge state
    # temporary_hack to assign mw of assumed +1 charge states
    f.Formulas = mass_to_formula(f.A0MZ - 1.00727647)
    return f

def create_feature(d):

    '''
    Read text file and create features
    :param f:
    :return:
    '''

    #feature instance
    f = Feature()

    #data is comma separated
    d = d[0].split(',')

    #populate data
    f.RT = d[0]
    f.A0MZ = d[1]
    f.Intensity = d[5]
    f.FileID = d[7]

    return f



