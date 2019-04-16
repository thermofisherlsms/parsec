import json
import math
from lib.componentdetection.pdata import *
from lib.componentdetection.chemcalc_services import *
import numpy as np
from scipy import signal
from numpy import sqrt, pi, exp
import peakutils
from scipy.spatial import distance
from math import *

c13 = 1.003355

def square_rooted(x):
    return round(sqrt(sum([a * a for a in x])), 3)

def cosine_similarity(x, y):
    numerator = sum(a * b for a, b in zip(x, y))
    denominator = square_rooted(x) * square_rooted(y)
    return round(numerator / float(denominator), 3)

def read_file(file):
    with open(file) as data_file:
        data = json.load(data_file)
        return (data)

def get_pscan(spectrum):
    scan = PScan()
    scan.FileName = spectrum.FileName
    scan.BasePeakIntensity = spectrum.BasePeakIntensity
    scan.BasePeakMass = spectrum.BasePeakMass
    scan.LowMass = spectrum.LowMass
    scan.HighMass = spectrum.HighMass
    scan.Masses = spectrum.Masses
    scan.PacketType = spectrum.PacketType
    scan.RetentionTime = spectrum.RetentionTime
    scan.ScanNumber = spectrum.ScanNumber
    scan.ScanType = spectrum.ScanType
    scan.TIC = spectrum.TIC
    scan.Intensities = spectrum.Intensities
    scan.Baselines = spectrum.Baselines
    scan.Noises = spectrum.Noises
    scan.Resolutions = spectrum.Resolutions
    scan.Charges = spectrum.Charges
    scan.MsOrder = spectrum.MsOrder
    if spectrum.MsOrder >= 2:
        scan.PrecursorMass = spectrum.PrecursorMass
    return scan

def scan_from_json(spectrum):
    scan = PScan()
    for i in spectrum:
        if i[0] == 'FileName':
            scan.FileName = i[1]
        if i[0] == 'BasePeakIntensity':
            scan.BasePeakIntensity = i[1]
        if i[0] == 'BasePeakMass':
            scan.BasePeakMass = i[1]
        if i[0] == 'LowMass':
            scan.LowMass = i[1]
        if i[0] == 'HighMass':
            scan.HighMass = i[1]
        if i[0] == 'Masses':
            scan.Masses = i[1]
        if i[0] == 'PacketType':
            scan.PacketType = i[1]
        if i[0] == 'RetentionTime':
            scan.RetentionTime = i[1]
        if i[0] == 'ScanNumber':
            scan.ScanNumber = i[1]
        if i[0] == 'ScanType':
            scan.ScanType = i[1]
        if i[0] == 'TIC':
            scan.TIC = i[1]
        if i[0] == 'Intensities':
            scan.Intensities = i[1]
        if i[0] == 'Baselines':
            scan.Baselines = i[1]
        if i[0] == 'Noises':
            scan.Noises = i[1]
        if i[0] == 'Resolutions':
            scan.Resolutions = i[1]
        if i[0] == 'Charges':
            scan.Charges = i[1]
        if i[0] == 'PrecursorMass':
                scan.PrecursorMass = i[1]
    return scan


def get_scan_list(file):
    data = read_file(file)
    scanList = []
    for s in data:
        spectrum = s.items()
        scan = scan_from_json(spectrum)
        scanList.append(scan)
    return scanList

def get_ppm_error(t, m):
    """
    calculate theoretical (t) and measured (m) ppm
    """
    return (((t - m) / t) * 1e6)


def get_mass_given_ppm(t, p):
    """
    calculate theoretical mass (t) given a ppm tolerance (p)
    """
    return (t - ((p * t) / 1e6))


def find_nearest_mass(m, t, p):
    """
    find nearest mass in the mass list dictionary (m), target mass (t) and ppm tolerance (p)
    """
    for k, mass in m.items():
        if abs(get_ppm_error(mass, t)) < p:
            return k
    return None


def detect_peaks_mf(t, w):
    """
    Detect peaks from mass trace (t) where w is the expected peak width at base (choose odd value)    """
    s = []
    rt = []
    mz = []
    for i in t.Trace:
        s.append(i.Intensity)
        rt.append(i.RT)
        mz.append(i.Mass)
    m = signal.medfilt(s, w)
    r = s - m
    p = Peak()
    p.Intensity = max(r)
    p.RT = round(rt[s.index(max(s))], 1)
    p.MZ = mz[s.index(max(s))]
    p.FileID = t.FileID
    return p

def get_median(a):
    return np.median(a)

def get_max(a):
    return max(a)

def remove_zeros(x):
    """
    remove zeros from a list
    """
    return list(filter(lambda a: a != 0, x))

def get_mz_list(scan):
    """
    Create mz objects from a scan
    """
    l = []
    for i in range(len(scan.Masses)):
        m = MZ()
        m.FileID = scan.FileName
        m.ScanNumber = scan.ScanNumber
        m.RT = scan.RetentionTime
        m.Intensity = scan.Intensities[i]
        m.Mass = scan.Masses[i]
        l.append(m)
    return l

def get_mass_traces_from_bins(scansPerFile):

    fileId = scansPerFile[0]
    scanList = scansPerFile[1]

    binMasses = {}

    for scan in scanList:
        bin_mzs(scan, binMasses)

    traceList = []

    for k,bin in binMasses.items():
        traces = build_mass_traces(k, bin)
        #traceList.append(traces)
        for j in traces:
            traceList.append((fileId, j))

    return traceList


def bin_mzs(scan, binMasses):

    for i in range(len(scan.Masses)):
        binNumber = round(scan.Masses[i])

        if(binNumber in binMasses):
            massList = binMasses[binNumber]
        else:
            massList = []
            binMasses[binNumber] = massList

        mz = MZPlus()
        mz.FileID = id
        mz.Mass = scan.Masses[i]
        mz.Intensity = scan.Intensities[i]
        mz.RetentionTime = scan.RetentionTime
        mz.RT = scan.RetentionTime
        mz.ScanNumber = scan.ScanNumber
        massList.append(mz)

def build_mass_traces(binIndex, massBin):

    # sort by intenstity, get most intense mz and anything within 5ppm of that. Mark processed so not picked up again

    MaxPPMDiff = 10
    MinTracePoints = 2
    ScalingFactor = 0

    if (binIndex <= 500):  # scaling factor for smaller masses to increase PPM tolerance, bin = rounded mz value
        f = divmod(binIndex, 100)
        f = MaxPPMDiff - f[0]
        ScalingFactor = f

    intensitySort = sorted(massBin, key=lambda  mz: mz.Intensity, reverse=True)
    massSort = sorted(massBin, key=lambda  mz: mz.Mass)
    massSortLength = len(massSort)
    massTraceList = []

    for sortedMass in intensitySort:
        if(not sortedMass.Processed):
            startIndex = massSort.index(sortedMass)
            massTrace = MassTrace()
            massTrace.Trace = []
            massTrace.Trace.append(sortedMass)
            massTrace.HighestIntensity = sortedMass.Intensity
            massTrace.HighestIntensityMass = sortedMass.Mass
            massTrace.HighestIntensityRt = sortedMass.RT
            sortedMass.Processed = True
            ppmDiff = 0
            counter = 1;

            while((ppmDiff <= (MaxPPMDiff + ScalingFactor)) and (startIndex + counter < massSortLength)):
                nextMass = massSort[startIndex + counter]
                massDiff = get_absolute_mass_error(sortedMass.Mass, nextMass.Mass)

                if(massDiff <= (MaxPPMDiff + ScalingFactor)):
                    massTrace.Trace.append(nextMass)
                    nextMass.Processed = True
                    counter += 1
                else:
                    break

            ppmDiff = 0
            counter = 1;
            while ((ppmDiff <= (MaxPPMDiff + ScalingFactor)) and (startIndex - counter >= 0)):
                previousMass = massSort[startIndex - counter]
                massDiff = get_absolute_mass_error(sortedMass.Mass, previousMass.Mass)
                if (massDiff <= (MaxPPMDiff + ScalingFactor)):
                    massTrace.Trace.append(previousMass)
                    previousMass.Processed = True
                    counter += 1
                else:
                    break

            if(len(massTrace.Trace) >= MinTracePoints):
                massTrace.Trace = sorted(massTrace.Trace, key=lambda m: m.RetentionTime)
                massTraceList.append(massTrace)

    return massTraceList

def get_absolute_mass_error(t, m):
    """
    calculate theoretical (t) and measured (m) ppm
    """
    return abs(((t - m) / t) * 1e6)

def PeakUtilsDetect2(trace):
    peaks = []
    x = []
    y = []
    for i in range(len(trace)):
        x.append(trace[i].RT)
        y.append(trace[i].Intensity)
    x = np.asarray(x)
    y = np.asarray(y)

    indexes = peakutils.indexes(y, thres=0.005, min_dist=.02)

    for i in range(len(indexes)):
        p = Peak()
        p.MZ = trace[indexes[i]].Mass
        p.FileID = trace[0].FileID
        p.RT = x[indexes[i]]
        p.Intensity = y[indexes[i]]
        p.Charge = 1 #todo figure out charge value, assume 1 for now
        #need this in case the trace is zero padded
        m = [x for x in trace if x.Mass !=0]
        p.SN = stn = max([i.Intensity for i in m]) \
                     / min([i.Intensity for i in m])
        peaks.append(p)
    return peaks

def detect_peaks_w_derivatives(trace):
    peaks = []
    x = []
    y = []
    for i in range(len(trace)):
        x.append(trace[i].RT)
        y.append(trace[i].Intensity)
    x = np.asarray(x)
    y = np.asarray(y)
    vector = y
    indexes = signal.argrelextrema(np.array(vector), comparator=np.greater, order=3)
    for i in range(len(indexes[0])):
        p = Peak()
        rt_index = indexes[0][i]
        p.RT = x[rt_index]
        p.Intensity = y[rt_index]
        p.MZ = trace[rt_index].Mass
        p.FileID = trace[0].FileID
        p.Charge = 1  # todo figure out charge value, assume 1 for now
        peaks.append(p)
    return peaks


def indices_satisfying(target, peakList,err):

  return [index for index,p in enumerate(peakList) if abs(get_ppm_error(p.MZ,target + c13)) < err or
          abs(get_ppm_error(p.MZ,target + 2*c13)) < err  and len(peakList) != 0
          or abs(get_ppm_error(p.MZ,target + 3*c13)) < err and len(peakList) >= 2]

def detect_smallmolecule_features(keyedPeaks,err, assignFormula = False):

    features = []
    retentionTimeKey = keyedPeaks[0]
    peaks = keyedPeaks[1]
    intensitySortedPeaks = sorted(peaks, key=lambda p: p.Intensity, reverse=True)
    massSortedPeaks = sorted(peaks, key=lambda p: p.MZ)
    massPeaksLength = len(massSortedPeaks)
    # will tracked processed peaks via index
    peakIndexProcessList = [0] * len(massSortedPeaks)

    for intensityIndex in range(len(intensitySortedPeaks)):

        currentPeak = intensitySortedPeaks[intensityIndex]

        # get index in the mass sorted list of current peak, then increment index for checking next peak values
        peakList = []
        massIndex = massSortedPeaks.index(intensitySortedPeaks[intensityIndex])
        if (peakIndexProcessList[massIndex] == 1):
            continue
        massIndexCounter = massIndex + 1
        maxMass =  currentPeak.MZ + 3 * c13

        while((massIndexCounter < massPeaksLength) and  (massSortedPeaks[massIndexCounter ].MZ <= maxMass)):
            isotopePeak = massSortedPeaks[massIndexCounter]

            # check mass differences between current peak and potential isotopes, verify intensity also
            if abs(get_ppm_error(isotopePeak.MZ, currentPeak.MZ + c13)) < err and (len(peakList) == 0) and (isotopePeak.Intensity <= currentPeak.Intensity):
                peakList.append(isotopePeak)
                peakIndexProcessList[massIndexCounter] = 1
                massIndexCounter = massIndexCounter + 1
                continue
            if abs(get_ppm_error(isotopePeak.MZ, currentPeak.MZ + 2*c13)) < err and (len(peakList) == 1) and (isotopePeak.Intensity <= peakList[0].Intensity):
                peakList.append(isotopePeak)
                peakIndexProcessList[massIndexCounter] = 1
                massIndexCounter = massIndexCounter + 1
                continue
            if abs(get_ppm_error(isotopePeak.MZ, currentPeak.MZ + 3 * c13)) < err and (len(peakList) == 2) and (isotopePeak.Intensity <= peakList[1].Intensity):
                peakList.append(isotopePeak)
                peakIndexProcessList[massIndexCounter] = 1
                break
            massIndexCounter = massIndexCounter + 1


        # results = [p for p in intensitySortedPeaks if abs(GetPPMError(p.MZ,peak.MZ + c13)) < err and peak.Intensity > p.Intensity
        #            or abs(GetPPMError(p.MZ,peak.MZ + 2*c13)) < err and peak.Intensity > p.Intensity or
        #            abs(GetPPMError(p.MZ,peak.MZ + 3*c13)) < err and peak.Intensity > p.Intensity]

        # matt removed
        #inds = indices_satisfying(peakMz,intensitySortedPeaks,err)

        #sort by intensity
        isotopePeaks = sorted(peakList, key= lambda p: p.Intensity,reverse=True)

        #if the peak has no isotopes discard it
        if len(isotopePeaks) == 0:
            peakIndexProcessList[massIndex] = 1
        else:
            f = Feature()
            f.Peaks.append(intensitySortedPeaks[intensityIndex])
            f.FileID = intensitySortedPeaks[intensityIndex].FileID
            f.A0MZ = intensitySortedPeaks[intensityIndex].MZ
            f.RT = intensitySortedPeaks[intensityIndex].RT
            f.Intensity = intensitySortedPeaks[intensityIndex].Intensity
            f.Charge = intensitySortedPeaks[intensityIndex].Charge
            for p in isotopePeaks:
                f.Peaks.append(p)
            #prune feature. It is possible that A3 is added before A2
            # matt removed
            # ppmError = get_ppm_error((f.Peaks[0].MZ + c13), f.Peaks[1].MZ)
            # if abs(ppmError) < err:
            features.append(f)

            #TODO we will miss overlapping peaks and peaks that are within
            #the window are not handled yet
            #remove isotopes from original list
            # matt removed
            # for j in inds:
            #     #TODO strange index bug here to fix
            #     try:
            #         if j == len(intensitySortedPeaks):
            #             peakIndexProcessList[j-1] = 1
            #         else:
            #             if initialPass ==False:
            #                 peakIndexProcessList[j] = 1
            #                 initialPass = True
            #             else:
            #                 peakIndexProcessList[j-1] = 1
            #
            #     except:
            #         #catch the index out of range exception
            #         pass
            #remove presumed A0 from the original list
            peakIndexProcessList[massIndex] = 1

        # TODO we must deconvolute with charge state
        #  temporary_hack to assign mw of assumed +1 charge states
        if assignFormula:
            for f in features:
                f.formulas = mass_to_formula(f.A0 - 1.00727647)

    return (retentionTimeKey, features)


