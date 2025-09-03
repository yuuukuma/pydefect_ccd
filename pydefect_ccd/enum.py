# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from monty.json import MSONable
from vise.util.enum import ExtendedEnum


class Carrier(MSONable, ExtendedEnum):
    e, h = "e", "h"

    @property
    def pn(self):
        return "n" if self is Carrier.e else "p"

    @property
    def charge(self):
        if self == self.e:
            return -1
        elif self == self.h:
            return 1
        else:
            raise ValueError

    @classmethod
    def from_carrier_charge(cls, carrier_charge):
        if carrier_charge == 1:
            return cls.h
        elif carrier_charge == -1:
            return cls.e
        raise ValueError

    def is_occupied(self, occupation):
        return occupation > 0.1 if self is Carrier.e \
            else occupation < 0.9


class CorrectionType(MSONable, ExtendedEnum):
    extended_FNV = "extended FNV"
    kumagai2023 = "kumagai2023"


class BandEdge(MSONable, ExtendedEnum):
    vbm, cbm = "vbm", "cbm"

