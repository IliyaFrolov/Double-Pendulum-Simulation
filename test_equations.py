from equations import *
from numpy import sqrt
import pytest

@pytest.mark.parametrize('length, angle, coordinate, x_expected, y_expected', [
    (5, pi/4, 1, 5*sqrt(2)/2, -5*sqrt(2)/2),
    (3, 3*pi, 'X', 0, 3),
    (7, -19/4*pi, 'w', -7*sqrt(2)/2, 7*sqrt(2)/2)
])
def test_polar_to_cartesian(length, angle, coordinate, x_expected, y_expected):
    convert = Polar_to_Cartesian(length)
    
    assert round(convert(angle, 'x'), 3) == round(x_expected , 3)
    assert round(convert(angle, 'y'), 3) == round(y_expected, 3)
    with pytest.raises(Exception):
        convert(angle, 1)

@pytest.mark.parametrize('length_1, mass_1, length_2, mass_2, theta_1, omega_1, theta_2, omega_2', [
    (2, 3, 4, 5, pi, pi/2, pi/4, 2*pi), 
    (4, 4, 3, 1, 5*pi, 0, pi/2, 6),
    (0, 1, 10, 4, 2, 3*pi, 0, 2)
])
class TestEnergy():

    def test_kineticenergy(self, length_1, mass_1, length_2, mass_2, theta_1, omega_1, theta_2, omega_2):
        kinetic_1 = KineticEnergy(length_1, mass_1, length_2, mass_2, 1)
        kinetic_2 = KineticEnergy(length_1, mass_1, length_2, mass_2, 2)

        assert kinetic_1.L1 == length_1
        assert kinetic_1.m1 == mass_1
        assert kinetic_1.L2 == length_2
        assert kinetic_1.m2 == mass_2
        assert kinetic_1(theta_1, omega_1, theta_2, omega_2) == mass_1/2*length_1**2*omega_1**2
        assert kinetic_2(theta_1, omega_1, theta_2, omega_2) == mass_1/2*(length_1**2*omega_1**2+length_2**2*omega_2**2+2*length_1*length_2*omega_1*omega_2*cos(theta_1-theta_2))
        with pytest.raises(Exception):
            KineticEnergy(length_1, mass_1, length_2, mass_2, 3)

    def test_potentialenergy(self, length_1, mass_1, length_2, mass_2, theta_1, omega_1, theta_2, omega_2):
        potential_1 = PotentialEnergy(length_1, mass_1, length_2, mass_2, 1)
        potential_2 = PotentialEnergy(length_1, mass_1, length_2, mass_2, 2)

        assert potential_1.L1 == length_1
        assert potential_1.m1 == mass_1
        assert potential_1.L2 == length_2
        assert potential_1.m2 == mass_2
        assert potential_1(theta_1, theta_2) == mass_1*9.8*length_1*(1-cos(theta_1))
        assert potential_2(theta_1, theta_2) == mass_1*9.8*(length_1*(1-cos(theta_1))+length_2*(1-cos(theta_2)))
        with pytest.raises(Exception):
            PotentialEnergy(length_1, mass_1, length_2, mass_2, -5)