class Peak:
    def __init__(self):
        self.RT = 0
        self.PeakArea = 0
        self.Trace = MassTrace()

class PScan:
    def __init__(self):
        self.FileName = ""
        self.PrecursorMass = 0
        self.MsOrder = 0
        self.BasePeakIntensity = 0
        self.BasePeakMass = 0
        self.HighMass = 0
        self.LowMass = 0
        self.PacketType = 21
        self.RetentionTime = 0
        self.ScanNumber = 0
        self.ScanType = ""
        self.TIC = 0
        self.Baselines = []
        self.Intensities = []
        self.Masses = []
        self.Noises = []
        self.Resolutions = []
        self.Charges = []

class MassTrace:
    def __init__(self):
        self.FileID = ""
        self.TraceID = ""
        self.Trace = []
        self.HighestIntensity = 0
        self.HighestIntensityMass = 0
        self.HighestIntensityRt = 0


class MZ:
    def __init__(self):
        self.FileID = ""
        self.Mass = 0
        self.Intensity = 0
        self.RT = 0
        self.ScanNumber = 0

class MZPlus(MZ):
    def __init__(self):
        MZ.__init__(self)
        self.Processed = 0

class Peak:
    def __init__(self):
        self.FileID = ""
        self.RT = 0
        self.Intensity = 0
        self.MZ = 0
        self.ChiSq = 0
        self.SN = 0
        self.Charge = 0

    def as_dict(self):
        return {'PeakID': self.FileID, 'RT': self.RT, 'Intensity': self.Intensity, 'QualityScore': self.ChiSq}

class Feature:
    def __init__(self):
        self.FeatureID = ""
        self.FileID = ""
        self.Key = ""
        self.RT = 0
        self.Intensity = 0
        self.A0MZ = 0
        self.Charge = 0
        self.RtOffset= 0
        self.Distance = 0
        self.Peaks = []
        self.Formulas = []

class AdductFeature(Feature):
    def __init__(self):
        Feature.__init__(self)
        self.AdductName = ""
        self.Formulas = []
        self.FormulaScores = []

class Component:
    def __init__(self):
        self.ComponentID = ""
        self.ComponentMass = 0
        self.ComponentRt = 0
        self.ComponentIntensity = 0
        self.BaseFeature = ()
        self.AdductFeatures = []
        self.Formulas = []
        self.FormulaScores = []