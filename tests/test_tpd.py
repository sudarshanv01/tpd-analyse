
import pytest
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from tpd_analyse.tpd import PlotTPD
from ase.io import read
import numpy as np

@pytest.fixture
def generate_tpd_results():
    """ Generate a results dictionary for the tests. """

    files = ['tests/input/exposure_0p25.csv']
    atoms = read('tests/input/co.traj')
    thermo_gas = 0.00012 * np.array([2121.52, 39.91, 39.45])
    vibration_energies_gas = IdealGasThermo(thermo_gas, atoms = atoms,
            geometry='linear', symmetrynumber=1, spin=0)
    vibration_energies_ads = HarmonicThermo(0.00012 * np.array([2044.1, 282.2, 201.5, 188.5, 38.3, 11.5]))
    order = 1

    TPDClass = PlotTPD(exp_data=files,
                        order=order,
                        thermo_ads=vibration_energies_ads,
                        thermo_gas=vibration_energies_gas,
                        plot_temperature=np.linspace(100, 500, 50), 
                        T_switch=[170],
                        T_max=250,
                        T_rate_min=[250,300],
                        beta=2)
    # TPDClass.get_results()

    return TPDClass


def test_first_order(generate_tpd_results):
    """ Run a simple first order calculation and check if outputs are present. """
    # assert hasattr(generate_tpd_results, 'results')
    assert 'results' in generate_tpd_results.__dict__


