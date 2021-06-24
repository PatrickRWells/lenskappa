import os
from lenskappa.surveys import survey
import toml
import logging

class Survey:
    pass

class SurveyDataManager:
    def __init__(self):
        survey_location = (survey.__file__)
        survey_dir = os.path.dirname(survey_location)
        self._config_location = os.path.join(survey_dir, 'config')
    
    def _get_survey_config_(self ,survey_name):
        survey_filename = '.'.join([survey_name, 'toml'])
        survey_config_location = os.path.join(self._config_location, survey_filename)
        try:
            survey_data = toml.load(survey_config_location)
        except:
            logging.warning("Unable to find config file for suvey {} "\
                            "I will create one.".format(survey_name))
            survey_data = {}
        return survey_config_location, survey_data

        
    def add_surveydata(self, survey_name, data):
        survey_config_location, survey_data = self._get_survey_config_(survey_name)
        for key, val in data.items():
            self._recursively_add_data(survey_data, key, val)
        with open(survey_config_location, 'w') as f:
            toml.dump(survey_data, f)


    def _recursively_add_data(self, data, name, newdata):
        if name not in data.keys():
            if type(newdata) != dict:
                data.update({name: newdata})
            else:
                data.update({name: {}})
                for key, val in newdata.items():
                    self._recursively_add_data(data[name], key, val)
        else:
            if type(newdata) != dict:
                logging.warning("Got a value that already exists in config! It will be overwritten")
                data[name] = newdata
            else:
                for key, val in newdata.items():
                    self._recursively_add_data(data[name], key, val)



    def get_surveydata(self, surveyname, data):
        path, survey_config = self._get_survey_config_(surveyname)
        return self._recursively_get_data(survey_config, data)
            
    
    def _recursively_get_data(self, data, retrieve):
        return_data = {}
        if type(retrieve) == dict:
            for key, val in retrieve.items():
                if key not in data.keys():
                    logging.warning("Error: Config item {} not found in the survey config")
                else: 
                    ret_data = self._recursively_get_data(data[key], val)
                return_data.update({key: ret_data})

        elif type(retrieve) == list:
            for name in retrieve:
                try:
                    data_val = data[name]
                    return_data.update({name: data_val})
                except:
                    logging.warning("Error: Config item {} not found in the survey config".format(name))
        else:
            try:
                data_val = data[retrieve]
                return_data.update({retrieve: data_val})
            except:
                logging.warning("Error: Config item {} not found in the survey config".format(retrieve))
        return return_data
