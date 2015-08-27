from lmfit import Parameters

######## DEFINE PSF FUNCTION #######################

class PsfFunction:

    name = str
    gaussianFWHM = float
    moffatAlpha = float
    moffatBeta = float

######## DEFINE COMPONENT OBJECT ###################

class Component:

    number = int
    name = str
    parameters = Parameters()
