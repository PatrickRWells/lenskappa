import os
from lenskappa import surveys
import logging
import toml
import logging
import argparse
from pydoc import locate
import lenskappa
from copy import copy
import hashlib

class Survey:
    pass

class SurveyDataManager:

    def __init__(self, survey):
        """
        Basic configuration manager for keeping track of survey data
        Should probably never be explicitly instantiated by a user
        Keeps track of things with '.toml' files (found in /config)
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
        base = os.path.dirname(surveys.survey.__file__)
        fname = '.'.join([survey, 'toml'])
        self._survey_config_location = os.path.join(base, 'config', fname)
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
        except:
            logging.error("No data types were found for survey {}".format(self._survey))
            exit(1)
        
        try:
            metadata = self.data_inputs['meta']
            self._setup_metadata(metadata)
        except:
            logging.info("No metadata types found for this survey")
        
        try:
            self.data = configdata['data']
        except:
            self.data = {}
            
        self._configdata = configdata
    
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
    
    def _write_config(self):
        """
        Write the config in its current state to the config file
        """

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
            return self.data[hash_key]
        except:
            logging.error("File not found for parameters: {}".format(data))            
    
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
    