"""Static definitions for the compartments and model parameters.

This code is automatically generated by generator/make_constants.py

Any manual changes will be lost.
"""
from collections import (
    namedtuple,
)
import os

import numpy as np
from numba.core import (
    types,
)
from numba.typed import (
    Dict,
)


#######################
# Primitive variables #
#######################

_BaseCompartment = namedtuple('BaseCompartment', [
    'S', 
    'E', 
    'I', 
    'R', 
])

_Variant = namedtuple('Variant', [
    'ancestral', 
    'alpha',     
    'beta',      
    'gamma',     
    'delta',     
    'omega',     
])

_BaseParameter = namedtuple('BaseParameter', [
    'alpha', 
    'sigma', 
    'gamma', 
    'pi',    
    'new_e', 
])

_VariantParameter = namedtuple('VariantParameter', [
    'beta',  
    'kappa', 
    'rho',   
])

_RiskGroup = namedtuple('RiskGroup', [
    'lr', 
    'hr', 
])

_ProtectionStatus = namedtuple('ProtectionStatus', [
    'unprotected',          
    'non_escape_protected', 
    'escape_protected',     
    'omega_protected',      
])

_ImmuneStatus = namedtuple('ImmuneStatus', [
    'non_escape_immune', 
    'escape_immune',     
    'omega_immune',      
])

_VaccinationStatus = namedtuple('VaccinationStatus', [
    'unvaccinated', 
    'vaccinated',   
])

_RemovedVaccinationStatus = namedtuple('RemovedVaccinationStatus', [
    'unvaccinated',     
    'vaccinated',       
    'newly_vaccinated', 
])

_VaccineType = namedtuple('VaccineType', [
    'unprotected',          
    'non_escape_protected', 
    'escape_protected',     
    'omega_protected',      
    'non_escape_immune',    
    'escape_immune',        
    'omega_immune',         
])

BASE_COMPARTMENT = _BaseCompartment(*_BaseCompartment._fields)
VARIANT = _Variant(*_Variant._fields)
BASE_PARAMETER = _BaseParameter(*_BaseParameter._fields)
VARIANT_PARAMETER = _VariantParameter(*_VariantParameter._fields)
RISK_GROUP = _RiskGroup(*_RiskGroup._fields)
PROTECTION_STATUS = _ProtectionStatus(*_ProtectionStatus._fields)
IMMUNE_STATUS = _ImmuneStatus(*_ImmuneStatus._fields)
VACCINATION_STATUS = _VaccinationStatus(*_VaccinationStatus._fields)
REMOVED_VACCINATION_STATUS = _RemovedVaccinationStatus(*_RemovedVaccinationStatus._fields)
VACCINE_TYPE = _VaccineType(*_VaccineType._fields)

PARAMETERS = Dict.empty(
    types.UniTuple(types.unicode_type, 2),
    types.int8,
)
PARAMETERS[('alpha', 'all')] = np.int8(0)
PARAMETERS[('sigma', 'all')] = np.int8(1)
PARAMETERS[('gamma', 'all')] = np.int8(2)
PARAMETERS[('pi', 'all')] = np.int8(3)
PARAMETERS[('new_e', 'all')] = np.int8(4)
PARAMETERS[('beta', 'ancestral')] = np.int8(5)
PARAMETERS[('beta', 'alpha')] = np.int8(6)
PARAMETERS[('beta', 'beta')] = np.int8(7)
PARAMETERS[('beta', 'gamma')] = np.int8(8)
PARAMETERS[('beta', 'delta')] = np.int8(9)
PARAMETERS[('beta', 'omega')] = np.int8(10)
PARAMETERS[('kappa', 'ancestral')] = np.int8(11)
PARAMETERS[('kappa', 'alpha')] = np.int8(12)
PARAMETERS[('kappa', 'beta')] = np.int8(13)
PARAMETERS[('kappa', 'gamma')] = np.int8(14)
PARAMETERS[('kappa', 'delta')] = np.int8(15)
PARAMETERS[('kappa', 'omega')] = np.int8(16)
PARAMETERS[('rho', 'ancestral')] = np.int8(17)
PARAMETERS[('rho', 'alpha')] = np.int8(18)
PARAMETERS[('rho', 'beta')] = np.int8(19)
PARAMETERS[('rho', 'gamma')] = np.int8(20)
PARAMETERS[('rho', 'delta')] = np.int8(21)
PARAMETERS[('rho', 'omega')] = np.int8(22)
PARAMETERS_NAMES = ['_'.join(k) for k in PARAMETERS]

VACCINE_TYPES = Dict.empty(
    types.UniTuple(types.unicode_type, 1),
    types.int8,
)
VACCINE_TYPES[('unprotected',)] = np.int8(0)
VACCINE_TYPES[('non_escape_protected',)] = np.int8(1)
VACCINE_TYPES[('escape_protected',)] = np.int8(2)
VACCINE_TYPES[('omega_protected',)] = np.int8(3)
VACCINE_TYPES[('non_escape_immune',)] = np.int8(4)
VACCINE_TYPES[('escape_immune',)] = np.int8(5)
VACCINE_TYPES[('omega_immune',)] = np.int8(6)
VACCINE_TYPES_NAMES = ['_'.join(k) for k in VACCINE_TYPES]

COMPARTMENTS = Dict.empty(
    types.UniTuple(types.unicode_type, 3),
    types.int8,
)
COMPARTMENTS[('S', 'unprotected', 'unvaccinated')] = np.int8(0)
COMPARTMENTS[('S', 'unprotected', 'vaccinated')] = np.int8(1)
COMPARTMENTS[('S', 'non_escape_protected', 'unvaccinated')] = np.int8(2)
COMPARTMENTS[('S', 'non_escape_protected', 'vaccinated')] = np.int8(3)
COMPARTMENTS[('S', 'escape_protected', 'unvaccinated')] = np.int8(4)
COMPARTMENTS[('S', 'escape_protected', 'vaccinated')] = np.int8(5)
COMPARTMENTS[('S', 'omega_protected', 'unvaccinated')] = np.int8(6)
COMPARTMENTS[('S', 'omega_protected', 'vaccinated')] = np.int8(7)
COMPARTMENTS[('S', 'non_escape_immune', 'vaccinated')] = np.int8(8)
COMPARTMENTS[('S', 'escape_immune', 'vaccinated')] = np.int8(9)
COMPARTMENTS[('S', 'omega_immune', 'vaccinated')] = np.int8(10)
COMPARTMENTS[('E', 'ancestral', 'unvaccinated')] = np.int8(11)
COMPARTMENTS[('E', 'ancestral', 'vaccinated')] = np.int8(12)
COMPARTMENTS[('E', 'alpha', 'unvaccinated')] = np.int8(13)
COMPARTMENTS[('E', 'alpha', 'vaccinated')] = np.int8(14)
COMPARTMENTS[('E', 'beta', 'unvaccinated')] = np.int8(15)
COMPARTMENTS[('E', 'beta', 'vaccinated')] = np.int8(16)
COMPARTMENTS[('E', 'gamma', 'unvaccinated')] = np.int8(17)
COMPARTMENTS[('E', 'gamma', 'vaccinated')] = np.int8(18)
COMPARTMENTS[('E', 'delta', 'unvaccinated')] = np.int8(19)
COMPARTMENTS[('E', 'delta', 'vaccinated')] = np.int8(20)
COMPARTMENTS[('E', 'omega', 'unvaccinated')] = np.int8(21)
COMPARTMENTS[('E', 'omega', 'vaccinated')] = np.int8(22)
COMPARTMENTS[('I', 'ancestral', 'unvaccinated')] = np.int8(23)
COMPARTMENTS[('I', 'ancestral', 'vaccinated')] = np.int8(24)
COMPARTMENTS[('I', 'alpha', 'unvaccinated')] = np.int8(25)
COMPARTMENTS[('I', 'alpha', 'vaccinated')] = np.int8(26)
COMPARTMENTS[('I', 'beta', 'unvaccinated')] = np.int8(27)
COMPARTMENTS[('I', 'beta', 'vaccinated')] = np.int8(28)
COMPARTMENTS[('I', 'gamma', 'unvaccinated')] = np.int8(29)
COMPARTMENTS[('I', 'gamma', 'vaccinated')] = np.int8(30)
COMPARTMENTS[('I', 'delta', 'unvaccinated')] = np.int8(31)
COMPARTMENTS[('I', 'delta', 'vaccinated')] = np.int8(32)
COMPARTMENTS[('I', 'omega', 'unvaccinated')] = np.int8(33)
COMPARTMENTS[('I', 'omega', 'vaccinated')] = np.int8(34)
COMPARTMENTS[('R', 'ancestral', 'unvaccinated')] = np.int8(35)
COMPARTMENTS[('R', 'ancestral', 'vaccinated')] = np.int8(36)
COMPARTMENTS[('R', 'ancestral', 'newly_vaccinated')] = np.int8(37)
COMPARTMENTS[('R', 'alpha', 'unvaccinated')] = np.int8(38)
COMPARTMENTS[('R', 'alpha', 'vaccinated')] = np.int8(39)
COMPARTMENTS[('R', 'alpha', 'newly_vaccinated')] = np.int8(40)
COMPARTMENTS[('R', 'beta', 'unvaccinated')] = np.int8(41)
COMPARTMENTS[('R', 'beta', 'vaccinated')] = np.int8(42)
COMPARTMENTS[('R', 'beta', 'newly_vaccinated')] = np.int8(43)
COMPARTMENTS[('R', 'gamma', 'unvaccinated')] = np.int8(44)
COMPARTMENTS[('R', 'gamma', 'vaccinated')] = np.int8(45)
COMPARTMENTS[('R', 'gamma', 'newly_vaccinated')] = np.int8(46)
COMPARTMENTS[('R', 'delta', 'unvaccinated')] = np.int8(47)
COMPARTMENTS[('R', 'delta', 'vaccinated')] = np.int8(48)
COMPARTMENTS[('R', 'delta', 'newly_vaccinated')] = np.int8(49)
COMPARTMENTS[('R', 'omega', 'unvaccinated')] = np.int8(50)
COMPARTMENTS[('R', 'omega', 'vaccinated')] = np.int8(51)
COMPARTMENTS[('R', 'omega', 'newly_vaccinated')] = np.int8(52)
COMPARTMENTS_NAMES = ['_'.join(k) for k in COMPARTMENTS]

AGGREGATES = Dict.empty(
    types.UniTuple(types.unicode_type, 2),
    types.int8,
)
AGGREGATES[('S', 'ancestral')] = np.int8(0)
AGGREGATES[('S', 'alpha')] = np.int8(1)
AGGREGATES[('S', 'beta')] = np.int8(2)
AGGREGATES[('S', 'gamma')] = np.int8(3)
AGGREGATES[('S', 'delta')] = np.int8(4)
AGGREGATES[('S', 'omega')] = np.int8(5)
AGGREGATES[('E', 'ancestral')] = np.int8(6)
AGGREGATES[('E', 'alpha')] = np.int8(7)
AGGREGATES[('E', 'beta')] = np.int8(8)
AGGREGATES[('E', 'gamma')] = np.int8(9)
AGGREGATES[('E', 'delta')] = np.int8(10)
AGGREGATES[('E', 'omega')] = np.int8(11)
AGGREGATES[('I', 'ancestral')] = np.int8(12)
AGGREGATES[('I', 'alpha')] = np.int8(13)
AGGREGATES[('I', 'beta')] = np.int8(14)
AGGREGATES[('I', 'gamma')] = np.int8(15)
AGGREGATES[('I', 'delta')] = np.int8(16)
AGGREGATES[('I', 'omega')] = np.int8(17)
AGGREGATES[('R', 'ancestral')] = np.int8(18)
AGGREGATES[('R', 'alpha')] = np.int8(19)
AGGREGATES[('R', 'beta')] = np.int8(20)
AGGREGATES[('R', 'gamma')] = np.int8(21)
AGGREGATES[('R', 'delta')] = np.int8(22)
AGGREGATES[('R', 'omega')] = np.int8(23)
AGGREGATES[('N', 'unvaccinated')] = np.int8(24)
AGGREGATES[('N', 'vaccinated')] = np.int8(25)
AGGREGATES_NAMES = ['_'.join(k) for k in AGGREGATES]

NEW_E = Dict.empty(
    types.UniTuple(types.unicode_type, 1),
    types.int8,
)
NEW_E[('ancestral',)] = np.int8(0)
NEW_E[('alpha',)] = np.int8(1)
NEW_E[('beta',)] = np.int8(2)
NEW_E[('gamma',)] = np.int8(3)
NEW_E[('delta',)] = np.int8(4)
NEW_E[('omega',)] = np.int8(5)
NEW_E_NAMES = ['_'.join(k) for k in NEW_E]

FORCE_OF_INFECTION = Dict.empty(
    types.UniTuple(types.unicode_type, 1),
    types.int8,
)
FORCE_OF_INFECTION[('ancestral',)] = np.int8(0)
FORCE_OF_INFECTION[('alpha',)] = np.int8(1)
FORCE_OF_INFECTION[('beta',)] = np.int8(2)
FORCE_OF_INFECTION[('gamma',)] = np.int8(3)
FORCE_OF_INFECTION[('delta',)] = np.int8(4)
FORCE_OF_INFECTION[('omega',)] = np.int8(5)
FORCE_OF_INFECTION_NAMES = ['_'.join(k) for k in FORCE_OF_INFECTION]

NATURAL_IMMUNITY_WANED = Dict.empty(
    types.UniTuple(types.unicode_type, 1),
    types.int8,
)
NATURAL_IMMUNITY_WANED[('ancestral',)] = np.int8(0)
NATURAL_IMMUNITY_WANED[('alpha',)] = np.int8(1)
NATURAL_IMMUNITY_WANED[('beta',)] = np.int8(2)
NATURAL_IMMUNITY_WANED[('gamma',)] = np.int8(3)
NATURAL_IMMUNITY_WANED[('delta',)] = np.int8(4)
NATURAL_IMMUNITY_WANED[('omega',)] = np.int8(5)
NATURAL_IMMUNITY_WANED_NAMES = ['_'.join(k) for k in NATURAL_IMMUNITY_WANED]

VACCINE_IMMUNITY_WANED = Dict.empty(
    types.UniTuple(types.unicode_type, 1),
    types.int8,
)
VACCINE_IMMUNITY_WANED[('non_escape_immune',)] = np.int8(0)
VACCINE_IMMUNITY_WANED[('escape_immune',)] = np.int8(1)
VACCINE_IMMUNITY_WANED[('omega_immune',)] = np.int8(2)
VACCINE_IMMUNITY_WANED_NAMES = ['_'.join(k) for k in VACCINE_IMMUNITY_WANED]

COMPARTMENT_GROUPS = Dict.empty(
    types.UniTuple(types.unicode_type, 2),
    types.int8[:],
)
COMPARTMENT_GROUPS[('S', 'ancestral')] = np.array([
    COMPARTMENTS[('S', 'unprotected', 'unvaccinated')],
    COMPARTMENTS[('S', 'unprotected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('S', 'alpha')] = np.array([
    COMPARTMENTS[('S', 'unprotected', 'unvaccinated')],
    COMPARTMENTS[('S', 'unprotected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('S', 'beta')] = np.array([
    COMPARTMENTS[('S', 'unprotected', 'unvaccinated')],
    COMPARTMENTS[('S', 'unprotected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('S', 'gamma')] = np.array([
    COMPARTMENTS[('S', 'unprotected', 'unvaccinated')],
    COMPARTMENTS[('S', 'unprotected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('S', 'delta')] = np.array([
    COMPARTMENTS[('S', 'unprotected', 'unvaccinated')],
    COMPARTMENTS[('S', 'unprotected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('S', 'omega')] = np.array([
    COMPARTMENTS[('S', 'unprotected', 'unvaccinated')],
    COMPARTMENTS[('S', 'unprotected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'escape_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'unvaccinated')],
    COMPARTMENTS[('S', 'omega_protected', 'vaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'unvaccinated')],
    COMPARTMENTS[('S', 'non_escape_immune', 'vaccinated')],
    COMPARTMENTS[('S', 'escape_immune', 'unvaccinated')],
    COMPARTMENTS[('S', 'escape_immune', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('E', 'ancestral')] = np.array([
    COMPARTMENTS[('E', 'ancestral', 'unvaccinated')],
    COMPARTMENTS[('E', 'ancestral', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('E', 'alpha')] = np.array([
    COMPARTMENTS[('E', 'alpha', 'unvaccinated')],
    COMPARTMENTS[('E', 'alpha', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('E', 'beta')] = np.array([
    COMPARTMENTS[('E', 'beta', 'unvaccinated')],
    COMPARTMENTS[('E', 'beta', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('E', 'gamma')] = np.array([
    COMPARTMENTS[('E', 'gamma', 'unvaccinated')],
    COMPARTMENTS[('E', 'gamma', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('E', 'delta')] = np.array([
    COMPARTMENTS[('E', 'delta', 'unvaccinated')],
    COMPARTMENTS[('E', 'delta', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('E', 'omega')] = np.array([
    COMPARTMENTS[('E', 'omega', 'unvaccinated')],
    COMPARTMENTS[('E', 'omega', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('I', 'ancestral')] = np.array([
    COMPARTMENTS[('I', 'ancestral', 'unvaccinated')],
    COMPARTMENTS[('I', 'ancestral', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('I', 'alpha')] = np.array([
    COMPARTMENTS[('I', 'alpha', 'unvaccinated')],
    COMPARTMENTS[('I', 'alpha', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('I', 'beta')] = np.array([
    COMPARTMENTS[('I', 'beta', 'unvaccinated')],
    COMPARTMENTS[('I', 'beta', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('I', 'gamma')] = np.array([
    COMPARTMENTS[('I', 'gamma', 'unvaccinated')],
    COMPARTMENTS[('I', 'gamma', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('I', 'delta')] = np.array([
    COMPARTMENTS[('I', 'delta', 'unvaccinated')],
    COMPARTMENTS[('I', 'delta', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('I', 'omega')] = np.array([
    COMPARTMENTS[('I', 'omega', 'unvaccinated')],
    COMPARTMENTS[('I', 'omega', 'vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('R', 'ancestral')] = np.array([
    COMPARTMENTS[('R', 'ancestral', 'unvaccinated')],
    COMPARTMENTS[('R', 'ancestral', 'vaccinated')],
    COMPARTMENTS[('R', 'ancestral', 'newly_vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('R', 'alpha')] = np.array([
    COMPARTMENTS[('R', 'alpha', 'unvaccinated')],
    COMPARTMENTS[('R', 'alpha', 'vaccinated')],
    COMPARTMENTS[('R', 'alpha', 'newly_vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('R', 'beta')] = np.array([
    COMPARTMENTS[('R', 'beta', 'unvaccinated')],
    COMPARTMENTS[('R', 'beta', 'vaccinated')],
    COMPARTMENTS[('R', 'beta', 'newly_vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('R', 'gamma')] = np.array([
    COMPARTMENTS[('R', 'gamma', 'unvaccinated')],
    COMPARTMENTS[('R', 'gamma', 'vaccinated')],
    COMPARTMENTS[('R', 'gamma', 'newly_vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('R', 'delta')] = np.array([
    COMPARTMENTS[('R', 'delta', 'unvaccinated')],
    COMPARTMENTS[('R', 'delta', 'vaccinated')],
    COMPARTMENTS[('R', 'delta', 'newly_vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('R', 'omega')] = np.array([
    COMPARTMENTS[('R', 'omega', 'unvaccinated')],
    COMPARTMENTS[('R', 'omega', 'vaccinated')],
    COMPARTMENTS[('R', 'omega', 'newly_vaccinated')],
], dtype=np.int8)
COMPARTMENT_GROUPS[('N', 'unvaccinated')] = np.array([
    v for k, v in COMPARTMENTS.items() if k[2] == 'unvaccinated'
])
COMPARTMENT_GROUPS[('N', 'vaccinated')] = np.array([
    v for k, v in COMPARTMENTS.items() if k[2] == 'vaccinated'
])
COMPARTMENT_GROUPS[('N', 'vaccine_eligible')] = np.array([
    v for k, v in COMPARTMENTS.items() if k[2] == 'unvaccinated' and k[0] not in ['E', 'I']
])
COMPARTMENT_GROUPS[('N', 'total')] = np.array([
    v for k, v in COMPARTMENTS.items()
])

# Turning off the JIT is operationally 1-to-1 with
# saying something is broken in the ODE code and
# I need to figure it out.
DEBUG = int(os.getenv('NUMBA_DISABLE_JIT', 0))
