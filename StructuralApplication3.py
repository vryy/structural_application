# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosStructuralApplication import *
application = KratosStructuralApplication()
application_name = "KratosStructuralApplication"

_ImportApplication(application, application_name)
