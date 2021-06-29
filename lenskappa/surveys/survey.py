import os

from numpy.core.numeric import full
import lenskappa
import logging
import toml
import logging
import argparse
from pydoc import locate
from copy import copy
import hashlib
from abc import ABCMeta, abstractmethod


class Survey(metaclass=ABCMeta):
    """
    Individual surveys are essentially plugins. They must have
    certain attributes and methods, to ensure they behave correctly
    in the core code, but how they do those things is totally up to them

    A survey should consist of several things:
        A Region, defining where the survey is on the sky
        A Catalog, containing the objects
        Optionally, a bright star mask
        A SurveyDataManger, defined as a class attribute

        It's up to the individual class to decide when and how to load these.
        It must define a setup method for this

        It should also have several methods for use during the weighted number counting
            mask_external_catalog: Mask a catalog that is NOT part of the survey,
                based on some defined region WITHIN the survey
            get_objects: Get objects from the survey's catalog, based on a region inside the survey area.
                            Should apply the bright star mask if it exists.
    """
    def __init__(self, *args, **kwargs):
        if not hasattr(self, "datamanager"):
            logging.critical("No data manager found for the survey!")
            return None
        self.setup(*args, **kwargs)
        self._validate()

    def _validate(self):
        try:
            region = self._region
        except:
            logging.error("No region found for the survey")
        
        try:
            catalog = self._catalog
        except:
            logging.error("Now catalog found for this survey")

    @abstractmethod
    def setup(self, *args, **kwargs):
        pass

    def generate_circular_tile(self, radius):
        """
        This should probably be overridden for some kinds of 
        """

        return self._region.generate_circular_tile(radius)

    @abstractmethod
    def mask_external_catalog(self, external_catalog, external_catalog_region, internal_region, *args, **kwargs):
        """
        Apply the bright star mask for a region inside the survey to a catalog from outside the survey region.
        
        Parameters:
            external_catalog: <catalog.Catalog> The catalog for the external objects
            external_region: <region.SkyRegion> A region defining the location of the catalog catalog
            internal_region: <region.SkyRegion> A region defining the location inside the survey to get the masks from
        """
    
    @abstractmethod
    def get_objects(self, internal_region, mask = True, *args, **kwargs):
        """
        Get objects within in a particular region of the survey.
        Either with or without masking objects near brigh stars.

        Parameters:
            internal_region <region.SkyRegion> Region inside the survey area to get objects for
            mask: <bool> Whether or not to mask out objects based on the bright star masks
        
        """
        pass
    



class SurveyDataManager:

    def __init__(self, survey):
        """
        Basic configuration manager for keeping track of survey data
        Should probably never be explicitly instantiated by a user
        Keeps track of things with '.toml' files (found in /config)

        It could be fun to spin this off into its own library for
        managing datasets on a local machine
        """
        self._cmds = ['add', 'remove']
        self._get_survey_config(survey)
    def __del__(self):
        pass
    
    def _get_survey_config(self, survey):
        """
        Gets the survey config file.

        Parameters:
            survey: <str> Name of the survey. It's expected that the
                    config file will be named survey.toml        
        """
        self._survey = survey
        base = os.path.dirname(lenskappa.__file__)
        self._basepath = os.path.join(base, 'surveys')
        base_configpath = os.path.join(base, 'surveys', 'config')
        fname = '.'.join([survey, 'toml'])
        self._survey_config_location = os.path.join(base_configpath, fname)
        try:
            survey_config = toml.load(self._survey_config_location)
            self._validate_survey_config(survey_config)

        except:
            logging.critical("Unable to find config file for suvey {} ".format(survey))
            exit(1) #PW: Using exit here means we should always try to load survey data early
    
    def _validate_survey_config(self, configdata):
        """
        Validate the configuration file for the survey. Make sure it behaves
        According to a few constraints
        """
        try:
            self.data_inputs = configdata['data_config']
            self._configdata = configdata

        except:
            logging.error("No data types were found for survey {}".format(self._survey))
            exit(1)
        
        try:
            metadata = self.data_inputs['meta']
            self._setup_metadata(metadata)
        except:
            logging.info("No metadata types found for this survey")
        
        try:
            support_data = configdata['support_data']
            self._setup_supportdata(support_data)
        except Exception as e:
            logging.info("No support data found for this survey")

        try:
            self.data = configdata['data']
        except:
            self.data = {}
            self._configdata.update({'data': self.data})
    
    def _setup_metadata(self, metadata):
        """
        Save metadata, and setup an argparser that can handle
        inputs from the various scripts
        """
        self._metadata = metadata


        self.argparser = argparse.ArgumentParser(description="hsc")
        for cmd, cmdata in self._metadata.items():
            try:
                required = cmdata['required']
                cmd_input = copy(cmdata)
                cmd_input.pop('required')
            except:
                cmd_input = cmdata

            self.argparser.add_argument(dest=cmd, **cmd_input)
        self.argparser.add_argument(dest="path", nargs=1, help="Path to the datafile")
    
    def _setup_supportdata(self, support_data):
        self._support_data = {}
        refs = {}
        for dataname, datavalue in support_data.items():
            try:
                path = datavalue['path']
                format = datavalue['format']
            except:
                continue

            for dataid, datainfo in datavalue.items():
                if datainfo.startswith('ref:'):
                    ref_key = '.'.join(['support_data', dataname, dataid])
                    self._parse_config_ref(ref_key, datainfo)
            self._support_data.update({dataname: datavalue})



    def _parse_config_ref(self, ref_key, ref_value):
        dict_vals = ref_value.split(':')[1]
        try:
            refs = self._config_refs
        except:
            self._config_refs = {}

        self._config_refs.update({ref_key: ref_value})
        ref_path = dict_vals.split('.')
        return_val = self._configdata
        for path in ref_path:
            return_val = return_val[path]

        ref_path = ref_key.split('.')
        item = self._configdata
        for key in ref_path[:-1]:
            item = item[key]
        item.update({ref_path[-1]: return_val})
        



    def _write_config(self):
        """
        Write the config in its current state to the config file
        """
        try:
            config_ref = self._config_refs
        except:
            config_ref = {}
        for config_key, ref_value in config_ref.items():
            keys = config_key.split('.')
            loc = self._configdata
            for key in keys[:-1]:
                loc = loc[key]
            
            loc.update({keys[-1]: ref_value})

        with open(self._survey_config_location, 'w') as f:
            toml.dump(self._configdata, f)

    def add(self, data, *args, **kwargs):
        """
        Add an datafile to a given survey. Uses a hash table to avoid
        having to name the files.

        Parameters:
        data {<str>: <str>}: Key, value pairs, where the key corresponds to one of the types of metadata
            defined in the suvey's config file. All metadata keys must be present
        """
        cwd = os.getcwd()
        datapath = os.path.join(cwd, data['path'])
        if not os.path.exists(datapath):
            logging.error("File {} does not exist".format(datapath))
            return
        keys = [data[key] for key in self._metadata.keys()]
        key = '_'.join(keys)
        hash_key = hashlib.md5(key.encode('utf-8')).hexdigest()
        if hash_key in self.data:
            logging.warning("Already have a path to this file. Overwriting...")
        
        self.data.update({hash_key: datapath})
        self._write_config()


    def remove(self, data):
        """
        Remove a data file from a given survey. Uses a hash table to avoid
        having to name the files.

        Parameters:
        data {<str>: <str>}: Key, value pairs, where the key corresponds to one of the types of metadata
            defined in the suvey's config file. All metadata keys must be present
        """

        keys = [data[key] for key in self._metadata.keys()]
        key = '_'.join(keys)
        hash_key = hashlib.md5(key.encode('utf-8')).hexdigest()
        try:
            self.data.pop(hash_key)
        except:
            logging.warning("No entry in the database for this file type")
        self._write_config()



    def get_file_location(self, data):
        """
        Gets the location of a particular file for survey.

        Parameters:
        data {<str>: <str>}: Key, value pairs, where the key corresponds to one of the types of metadata
            defined in the suvey's config file. All metadata keys must be present
        """
        keys = [data[key] for key in self._metadata.keys()]
        key = '_'.join(keys)
        hash_key = hashlib.md5(key.encode('utf-8')).hexdigest()
        try:
            full_path = self.data[hash_key]
            if full_path.endswith('.'):
                full_path = full_path[:-1]

            return full_path
        except:
            logging.error("File not found for parameters: {}".format(data))
    
    def get_support_file_location(self, data, *args, **kwargs):
        try:
            support_datatype = data['type']
        except:
            logging.error("No support data type found")
            return
        
        try:
            support_data_id = data['id']
        except:
            logging.error("No ID passed for support data retrieval")
            return

        if support_data_id in self._support_data[support_datatype]['ids']:
            fname = '.'.join([support_data_id, self._support_data[support_datatype]['format']])
            full_path =  os.path.join(self._basepath, self._support_data[support_datatype]['path'], fname)
            if full_path.endswith('.'):
                full_path = full_path[:-1]
            if os.path.exists(full_path):
                return full_path
            else:
                logging.warning("Unable to find support file {}".format(data))
                return

    def parse_cmd_input(self, args):

        """
        Parse command line inputs
        """
        cmd = args[0]
        if cmd == '-h' or cmd == '--help':
            output = self.argparser.parse_args([cmd])
            return

        input_args = args[1:]
        try:
            fn = getattr(self, cmd)
        except:
            logging.warning("Command {} not defined".format(cmd))
            return

        try:
            parsed_args = self.argparser.parse_args(input_args)
        except Exception as e:
            logging.warning("{} can't parse arguments".format(self._survey))
            logging.warning(e)

        arg_dict = {}
        for name, value in vars(parsed_args).items():
            try:
                if len(value) == 1:
                    arg_dict.update({name: value[0]})
            except:
                pass
                #Handle multi-input arguments
                
        return fn(arg_dict)
    
if __name__ == "__main__":
    manager = SurveyDataManager('hsc')
