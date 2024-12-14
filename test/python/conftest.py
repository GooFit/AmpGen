import pytest 
from ampgen import make_ampgen_model 
from pathlib import Path 
import json 
import os 

def pytest_addoption(parser):
    parser.addoption("--model", action="store", default=None, help="my option: type1 or type2")

def make_testname(case): 
    return (f"{case[2]}/{case[3]}")

def pytest_generate_tests(metafunc):
    option_value = metafunc.config.option.model
    if 'model' in metafunc.fixturenames: 
        model = make_ampgen_model(option_value, src=Path(option_value).stem +".cpp") 
        json_fname = os.getenv("AMPGENROOT") + "/test/ref/" + Path(option_value).stem + ".json"
        values = [] 
        for l in model.labels : 
            for evt in range(0,16) : 
                values.append( (model, json_fname, l, evt ) )
        metafunc.parametrize("model", values, ids = make_testname ) 
