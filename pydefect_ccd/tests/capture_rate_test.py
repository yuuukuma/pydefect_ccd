# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import pytest

from dephon.capture_rate import CaptureRate


@pytest.fixture
def capture_rate():
    return CaptureRate(Wif=1.23456,
                       summed_phonon_overlaps=[0.1234, 1.1234, 2.1234],
                       velocities=[10.001, 10.0, 10.0],
                       temperatures=[100.01234, 200, 300],
                       site_degeneracy=4.0,
                       spin_selection_factor=0.5,
                       volume=1000.1234)


def test_capture_rate_str(capture_rate):
    print(capture_rate)

