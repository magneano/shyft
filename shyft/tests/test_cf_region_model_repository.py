import os
import unittest
import shutil, tempfile
from shyft.repository.netcdf.cf_region_model_repository import CFRegionModelRepository
from shyft import shyftdata_dir
from shyft.api.pt_gs_k import PTGSKModel



class CFRegionModelRepositoryTestCase(unittest.TestCase):


    # Create a temporary directory
    test_dir = tempfile.mkdtemp()
    region = {'region_model_id': 'test',  # a unique name identifier of the simulation
              'domain': {'EPSG': 32633,
                         'nx': 400,
                         'ny': 80,
                         'step_x': 1000,
                         'step_y': 1000,
                         'lower_left_x': 100000,
                         'lower_left_y': 6960000},
              'repository': {'class': CFRegionModelRepository,
                             'params': {
                                 'data_file': os.path.join(shyftdata_dir,
                                                        'netcdf/orchestration-testdata/cell_data.nc')}},
              }

    model = {'model_t': PTGSKModel,  # model to construct
             'model_parameters': {
                 'ae': {
                     'ae_scale_factor': 1.5},
                 'gs': {
                     'calculate_iso_pot_energy': False,
                     'fast_albedo_decay_rate': 6.752787747748934,
                     'glacier_albedo': 0.4,
                     'initial_bare_ground_fraction': 0.04,
                     'max_albedo': 0.9,
                     'max_water': 0.1,
                     'min_albedo': 0.6,
                     'slow_albedo_decay_rate': 37.17325702015658,
                     'snow_cv': 0.4,
                     'snow_cv_altitude_factor': 0.0,
                     'snow_cv_forest_factor': 0.0,
                     'tx': -0.5752881492890207,
                     'snowfall_reset_depth': 5.0,
                     'surface_magnitude': 30.0,
                     'wind_const': 1.0,
                     'wind_scale': 1.8959672005350063,
                     'winter_end_day_of_year': 100},
                 'kirchner': {
                     'c1': -3.336197322290274,
                     'c2': 0.33433661533385695,
                     'c3': -0.12503959620315988},
                 'p_corr': {
                     'scale_factor': 1.0},
                 'pt': {'albedo': 0.2,
                        'alpha': 1.26},
                 'routing': {
                     'alpha': 0.9,
                     'beta': 3.0,
                     'velocity': 0.0}
             }
             }

    region_model_repo = CFRegionModelRepository(region, model)


    def test_get_region_model(self):

        region_model = self.region_model_repo.get_region_model('test')
        self.assertIsInstance(region_model, PTGSKModel, 'Correct model type not returned from CFRegionModelRepository')

    def test_cell_data_to_netcdf(self):
        region_model = self.region_model_repo.get_region_model('test')
        self.region_model_repo.cell_data_to_netcdf(region_model, os.path.join(self.test_dir,'test'))

        # open the file and be sure it works
        output_nc = os.path.join(self.test_dir,'test_cell_data.nc')
        self.region['repository']['params']['data_file'] = output_nc

        tmp_rm = CFRegionModelRepository(self.region, self.model).get_region_model('test')
        self.assertIsInstance(tmp_rm, PTGSKModel, 'Error with {}'.format(output_nc))

        shutil.rmtree(self.test_dir)


if __name__ == '__main__':
    unittest.main()
