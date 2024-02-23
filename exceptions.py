class AutodeException(Exception):
    """Base autodE exception"""

class CouldNotPlotSmoothProfile(AutodeException):
    """A smooth reaction profile cannot be plotted"""

    