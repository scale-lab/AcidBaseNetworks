# -*- coding: utf-8 -*-


from . import analysis
from . import automation
from . import chemtheory
from . import coding
from . import optimize
from . import simulation
from . import synthesis
from . import tools

from .tools.CompoundList import CompoundList
from .tools.MixtureList import MixtureList
from .tools.Containers import Container,WellPlate384,\
							 WellPlate96,WellPlate24,\
							 microtube20,microtube15,\
							 tube50,tube15,\
                                 WellPlate1536LDV,WellPlate384LDV,\
                                 WellPlate384PP,MaldiPlate384,MaldiPlate1536
from .tools.Task import Task,TransferTask,TaskList,\
						rowsum_TaskList,blocksum_TaskList,\
						maldispot_TaskList,dilution_TaskList
from .tools import Containers
from .tools import Plating
from .tools import Path
from .synthesis.UgiLibrary import UgiLibrary
from .synthesis.PasseriniLibrary import PasseriniLibrary
from .analysis import MassAnalysis
from .synthesis import SynthesisUtilities
from .synthesis.PlateSheet import PlateSheet
from .automation.Andrew import Andrew
from .automation.Echo import Echo
from .automation.Bruker import Bruker

from .coding.MixtureCoder import MixtureCoder


from sys import platform
if platform == "win32":
    from .tools import BrukerMS
