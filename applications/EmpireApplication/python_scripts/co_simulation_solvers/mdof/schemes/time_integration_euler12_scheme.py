from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from time_integration_base_scheme import TimeIntegrationBaseScheme

# Importing tools
from co_simulation_tools import ValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateScheme(scheme_settings):
    return TimeIntegrationEuler12Scheme(scheme_settings)

# PMT to be checked, seems to be written acceleration based, might need correction

class TimeIntegrationEuler12Scheme(TimeIntegrationBaseScheme):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """
    def __init__(self, scheme_settings):

        default_settings = {
                "type"      : "euler12",
                "settings"  : {}
            }

        # add buffer size - this is not user-specified
        # each derived scheme specifies it
        scheme_settings["settings"].update({"buffer_size":3})

        ValidateAndAssignDefaults(default_settings, scheme_settings)

        # base scheme settings
        super(TimeIntegrationEuler12Scheme, self).__init__(scheme_settings["settings"])

    def _AssembleLHS(self, model):
        """
        """
        return model.m

    def _AssembleRHS(self, model):
        """
        """
        self.buffer[3,0,:] = self.force

        RHS = self.buffer[3,1,:] * self.dt**2
        RHS += np.dot(2 * model.m - model.k * self.dt**2, self.buffer[0,1,:])
        RHS += np.dot(- model.m + 0.5 * model.b * self.dt, self.buffer[0,2,:])

        return RHS

    def UpdateDerivedValues(self):
        """
        """
        self.buffer[1,0,:] = (self.buffer[0,0,:] - self.buffer[0,2,:]) * 0.5 / self.dt
        self.buffer[2,0,:] = (self.buffer[0,0,:] - 2 * self.buffer[0,1,:] + self.buffer[0,2,:]) / self.dt**2