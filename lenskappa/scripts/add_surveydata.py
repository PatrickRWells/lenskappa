import os
import lenskappa
from lenskappa.surveys.survey import Survey, SurveyDataManager
import argparse
import toml

def add_surveydata(args):
    parser = argparse.ArgumentParser(description='ArgParser for add_surveydata')
    parser.add_argument('survey_name', type=str, nargs = 1, help = 'Name of the survey')
    parser.add_argument('data_type', type=str, nargs=1, help='Type of data being added (should be mask or catalog')
    parser.add_argument('path', type=str, nargs=1, help = 'Path to the data')
    argdata = parser.parse_args(args)
    
    payload = {'data': {argdata.data_type[0]: {'type': argdata.data_type[0], 'path': argdata.path[0]} } }
    manager = SurveyDataManager()
    manager.add_surveydata(argdata.survey_name[0], payload)