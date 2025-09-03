# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import Dict

import yaml
from dephon.enum import CorrectionType
from monty.json import MSONable
from vise.util.mix_in import ToYamlFileMixIn


@dataclass
class DephonCorrection(MSONable, ToYamlFileMixIn):
    corrections: Dict[CorrectionType, float]

    @property
    def total_correction_energy(self):
        return sum([v for v in self.corrections.values()])

    def to_yaml(self):
        return yaml.dump({str(k): v for k, v in self.corrections.items()})

    @classmethod
    def from_yaml(cls, filename: str):
        with open(filename) as file:
            d = yaml.safe_load(file)
        dd = {}
        for k, v in d.items():
            dd[CorrectionType.from_string(k)] = v
        return cls(dd)

