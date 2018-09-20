from os import path
import numpy as np

from . import interfaces
from .netcdf.cf_region_model_repository import BoundingBoxRegion
from shyft import api

from shyft import shyftdata_dir
from shyft.orchestration.configuration.config_interfaces import RegionConfig, ModelConfig, RegionConfigError
from shyft.orchestration.configuration.dict_configs import DictModelConfig, DictRegionConfig
from six import iteritems # This replaces dictionary.iteritems() on Python 2 and dictionary.items() on Python 3

class DeserializeRegionRipositoryError(Exception):
    pass

class DeserializeRegionRipository(interfaces.RegionModelRepository):
    """
    Repository that delivers fully specified shyft api region_models
    based on data found in netcdf files.
    """

    def __init__(self, region, model):
        """
        Parameters
        ----------
        region: either a dictionary suitable to be instantiated as a
            RegionConfig Object or a sublcass of the interface RegionConfig
            containing regional information, like
            catchment overrides, and which netcdf file to read
        model: either a dictionary suitable to be instantiated as a
            ModelConfig Object or a subclass of interface ModelConfig
            Object containing model information, i.e.
            information concerning interpolation and model
            parameters
        """

        if not isinstance(region, RegionConfig):
            region_config = DictRegionConfig(region)
        if not isinstance(model, ModelConfig):
            model_config = DictModelConfig(model)
        else:
            region_config = region
            model_config = model

        if not isinstance(region_config, RegionConfig) or \
           not isinstance(model_config, ModelConfig):
            raise interfaces.InterfaceError()
        self._rconf = region_config
        self._mconf = model_config
        self._region_model = model_config.model_type() # region_model
        self._mask = None
        self._epsg = self._rconf.domain()["EPSG"] # epsg
        filename = self._rconf.repository()["params"]["data_file"]
        filename = path.expandvars(filename)
        if not path.isabs(filename):
            # Relative paths will be prepended the data_dir
            filename = path.join(shyftdata_dir, filename)
        if not path.isfile(filename):
            raise DeserializeRegionRipositoryError("No such file '{}'".format(filename))
        self._data_file = filename
        self._catch_ids = self._rconf.catchments()
        self.bounding_box = None

    def get_region_model(self, region_id, catchments=None):
        """
        Return a fully specified shyft api region_model for region_id, based on data found
        in netcdf dataset.

        Parameters
        -----------
        region_id: string
            unique identifier of region in data

        catchments: list of unique integers
            catchment indices when extracting a region consisting of a subset
            of the catchments has attribs to construct params and cells etc.

        Returns
        -------
        region_model: shyft.api type
        """

        byte_vector = api.byte_vector_from_file(self._data_file)
        gcdv = api.GeoCellDataVector.deserialize(byte_vector) # geo cell data vector
        c_ids = [gcd.catchment_id() for gcd in gcdv]
        xcoord = np.array([gcd.mid_point().x for gcd in gcdv])
        ycoord = np.array([gcd.mid_point().y for gcd in gcdv])
        xcoord_m = xcoord
        ycoord_m = ycoord

        #m_catch = np.ones(len(c_ids), dtype=bool)
        #mask out cells that are not listed in config
        if self._catch_ids is not None:
            raise DeserializeRegionRipositoryError("Catchment id masking not yet available.")
        #    m_catch = np.in1d(c_ids, self._catch_ids)
        #    xcoord_m = xcoord[m_catch]
        #    ycoord_m = ycoord[m_catch]

        #dataset_epsg = None
        #if 'crs' in Vars.keys():
        #    dataset_epsg = Vars['crs'].epsg_code.split(':')[1] # TODO set dataset epsg to cell epsg
        #if not dataset_epsg:
        #    raise interfaces.InterfaceError("netcdf: epsg attr not found in group elevation")

        target_cs = "+init=EPSG:{}".format(self._epsg)
        #source_cs = "+init=EPSG:{}".format(dataset_epsg)
        source_cs = "+init=EPSG:{}".format(self._epsg)

        if not source_cs==target_cs:
            raise DeserializeRegionRipositoryError("Source crs ({}) must be equal to target crs ({})".format(source_cs,
                                                                                                             target_cs))

        # Construct bounding region
        # TODO:Throuw out cells now in bounding region
        box_fields = set(("lower_left_x", "lower_left_y", "step_x", "step_y", "nx", "ny", "EPSG"))
        if box_fields.issubset(self._rconf.domain()):
            raise DeserializeRegionRipositoryError("Bounding box masking not yet available.")
            #tmp = self._rconf.domain()
            #epsg = tmp["EPSG"]
            #x_min = tmp["lower_left_x"]
            #x_max = x_min + tmp["nx"]*tmp["step_x"]
            #y_min = tmp["lower_left_y"]
            #y_max = y_min + tmp["ny"]*tmp["step_y"]
            #bounding_region = BoundingBoxRegion(np.array([x_min, x_max]),
            #                                    np.array([y_min, y_max]), epsg, self._epsg)
        else:
            dataset_epsg = self._epsg # TODO: get epsg from cells and set to dataset epsg
            bounding_region = BoundingBoxRegion(xcoord_m, ycoord_m, dataset_epsg, self._epsg)
        self.bounding_box = bounding_region.bounding_box(self._epsg)

        # Construct region parameter:
        region_parameter = self._region_model.parameter_t()
        for p_type_name, value_ in iteritems(self._mconf.model_parameters()):
            if hasattr(region_parameter, p_type_name):
                sub_param = getattr(region_parameter, p_type_name)
                for p, v in iteritems(value_):
                    if hasattr(sub_param, p):
                        setattr(sub_param, p, v)
                    else:
                        raise RegionConfigError("Invalid parameter '{}' for parameter set '{}'".format(p, p_type_name))
            else:
                raise RegionConfigError("Invalid parameter set '{}' for selected model '{}'".format(p_type_name, self._region_model.__name__))

        # Construct cells
        cell_vector = self._region_model.cell_t.vector_t()
        for gcd in gcdv:
            cell = self._region_model.cell_t()
            cell.geo = gcd
            cell_vector.append(cell)

        # TODO: make sure target crs equals source crs. Find Crs in cell-vector!
        #target_cs = "+init=EPSG:{}".format(self._epsg)
        #source_cs = "+init=EPSG:{}".format(dataset_epsg)

        c_ids_unique = list(np.unique(c_ids))

        # Construct catchment overrides
        catchment_parameters = self._region_model.parameter_t.map_t()
        for cid, catch_param in iteritems(self._rconf.parameter_overrides()):
            if cid in c_ids_unique:
                param = self._region_model.parameter_t(region_parameter)
                for p_type_name, value_ in iteritems(catch_param):
                    if hasattr(param, p_type_name):
                        sub_param = getattr(param, p_type_name)
                        for p, v in iteritems(value_):
                            if hasattr(sub_param, p):
                                setattr(sub_param, p, v)
                            else:
                                raise RegionConfigError("Invalid parameter '{}' for catchment parameter set '{}'".format(p, p_type_name))
                    else:
                        raise RegionConfigError("Invalid catchment parameter set '{}' for selected model '{}'".format(p_type_name, self._region_model.__name__))

                catchment_parameters[cid] = param
        region_model = self._region_model(cell_vector, region_parameter, catchment_parameters)
        region_model.bounding_region = bounding_region
        region_model.catchment_id_map = c_ids_unique

        def do_clone(x):
            clone = x.__class__(x)
            clone.bounding_region = x.bounding_region
            clone.catchment_id_map = x.catchment_id_map
            # clone.gis_info = polygons  # cell shapes not included yet
            return clone

        region_model.clone = do_clone
        return region_model