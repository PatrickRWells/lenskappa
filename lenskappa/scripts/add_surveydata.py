import os
import lenskappa
from lenskappa.surveys.survey import Survey
from lenskappa import SurveyDataManager
import argparse
import logging
import toml

def add_surveydata(args):
    parser = argparse.ArgumentParser(description='ArgParser for add_surveydata')
    parser.add_argument('survey_name', type=str, nargs = 1, help = 'Name of the survey')
    parser.add_argument('-s', '--show-input', action='store_true', help = 'Show the expected inputs for the given survey')
    parser.add_argument('data', type=str, nargs='*')
    argdata = parser.parse_args(args)
    try:
        data = SurveyDataManager(argdata.survey_name[0])
    except Exception as e:
        logging.error("Couldn't find survey {}".format(argdata.survey_name[0]))
        return
    if argdata.show_input:
        input_args = ['-h']
    else:
        args[0] = 'add'
        input_args = args

    data.parse_cmd_input(input_args)
