import urllib.request, urllib.parse, urllib.error, urllib.request,urllib.error, urllib.parse
import json


def mass_to_formula(m):

    chemcalcURL = 'http://www.chemcalc.org/chemcalc/em'

    # Define a molecular formula string
    mfRange = 'C0-100H0-100N0-10O0-10'
    # target mass
    mass = m

    # Define the parameters and send them to Chemcalc
    # other options (mass tolerance, unsaturation, etc.
    params = {'mfRange': mfRange, 'monoisotopicMass': mass}

    response = urllib.request.urlopen(chemcalcURL, bytes(urllib.parse.urlencode(params), encoding="utf-8"))

    # Read the output and convert it from JSON into a Python dictionary
    jsondata = response.read()
    formulas = json.loads(jsondata)
    return formulas['results'][0:5]


def formula_info(f):

    ccurl = 'http://www.chemcalc.org/chemcalc/mf'

    # Define a molecular formula string
    mf = f

    # Define the parameters and send them to Chemcalc
    params = {'mf': mf, 'isotopomers': 'jcamp,xy'}
    response = urllib.request.urlopen(ccurl, urllib.parse.urlencode(params).encode('utf-8'))

    # Read the output and convert it from JSON into a Python dictionary
    jsondata = response.read()
    return json.loads(jsondata.decode('utf-8'))







