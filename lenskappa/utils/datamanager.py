import logging
import toml
import os
import argparse
from copy import copy
import hashlib
import appdirs
from sys import exit
import lenskappa

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

    def _get_survey_config(self, survey: str):
        """
        Gets the survey config file.

        Parameters:
            survey: <str> Name of the survey. It's expected that the
                    config file will be named survey.toml
        """
        self._survey = survey
        base = os.path.dirname(lenskappa.__file__)
        self._basepath = os.path.join(base, 'datasets/surveys', survey)
        fname = '.'.join([survey, 'toml'])
        self._survey_config_location = os.path.join(self._basepath, fname)
        try:
            print(self._survey_config_location)
            survey_config = toml.load(self._survey_config_location)
            self._validate_survey_config(survey_config)

        except Exception as e:
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

        self._load_data()

    def _load_data(self):
        """
        Loads any available data for the survey.
        Paths to the data are stored outside the main lenskappa library
        """
        top_dir = appdirs.user_data_dir()
        if not os.path.exists(top_dir):
            logging.CRITICAL("Unable to load data, the main application support path doesn't exist!")
            exit()
        lenskapa_data_dir = appdirs.user_data_dir("lenskappa")
        if not os.path.exists(lenskapa_data_dir):
            os.makedirs(lenskapa_data_dir)

        survey_fname = '.'.join([self._survey, 'toml'])
        self._survey_data_path = os.path.join(lenskapa_data_dir, survey_fname)
        if not os.path.exists(self._survey_data_path):
            self._survey_data = {}
        else:
            self._survey_data = toml.load(self._survey_data_path)
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

    def _write_data(self):
        with open(self._survey_data_path, 'w') as f:
            toml.dump(self._survey_data, f)

    def add(self, data, *args, **kwargs):
        """
        Add an datafile to a given survey. Uses a hash table to avoid
        having to name the files.
        This should probably be updated to use a database, rather than toml

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
        key = '_'.join([self._survey, key])
        hash_key = hashlib.md5(key.encode('utf-8')).hexdigest()
        if hash_key in self._survey_data:
            logging.warning("Already have a path to this file. Overwriting...")

        self._survey_data.update({hash_key: datapath})

        self._write_data()


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
        key = '_'.join([self._survey, key])
        hash_key = hashlib.md5(key.encode('utf-8')).hexdigest()
        try:
            self._survey_data.pop(hash_key)
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
        key = '_'.join([self._survey, key])
        hash_key = hashlib.md5(key.encode('utf-8')).hexdigest()
        try:
            full_path = self._survey_data[hash_key]
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
