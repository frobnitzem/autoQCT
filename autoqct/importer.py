import importlib.util
#import sys

def import_py(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    #sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module
