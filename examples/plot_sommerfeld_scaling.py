# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
"""Example plot of Sommerfeld scaling factors."""

from matplotlib import pyplot as plt

from pydefect_ccd.local_enum import Carrier
from pydefect_ccd.sommerfeld_scaling import SommerfeldScaling


def main():
    """Plot electron Sommerfeld scaling for several masses and charges."""
    fig, ax = plt.subplots()
    temperatures = [100, 200, 300, 400, 500]

    for mass in [0.1, 1.0]:
        for dielectric_constant, line_style \
                in zip([10], ["-", "--", ":"]):
            for defect_charge in [-1]:
                scaling = SommerfeldScaling(
                    epsilon0=dielectric_constant,
                    electron_effective_mass=mass,
                    hole_effective_mass=1.0,
                    Ts=temperatures,
                )
                scaling.add_to_ax(
                    ax,
                    carrier_type=Carrier.e,
                    defect_charge=defect_charge,
                    ls=line_style,
                )

    SommerfeldScaling.set_label(ax)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
