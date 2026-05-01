# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from monty.json import MSONable
from vise.util.enum import ExtendedEnum


class Carrier(MSONable, ExtendedEnum):
    """Carrier type used to select band edges and occupation criteria."""
    e, h = "e", "h"

    @property
    def carrier_type(self) -> str:
        """Return the semiconductor convention, ``n`` for e and ``p`` for h."""
        return "n" if self is Carrier.e else "p"

    @property
    def charge(self):
        """Return the carrier charge sign."""
        if self == self.e:
            return -1
        elif self == self.h:
            return 1
        else:
            raise ValueError

    @classmethod
    def from_carrier_charge(cls, carrier_charge):
        """Return a carrier from charge sign."""
        if carrier_charge == 1:
            return cls.h
        elif carrier_charge == -1:
            return cls.e
        raise ValueError

    def is_occupied(self, occupation):
        """Return whether a state occupation is consistent with this carrier."""
        return occupation > 0.1 if self is Carrier.e else occupation < 0.9
