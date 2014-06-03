# run appconfig by calling:
# $ python run.py appconfig/configs/configfile.cfg

import subprocess
import os
import os.path
from argparse import ArgumentParser

import time
import datetime

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')


parser = ArgumentParser()
parser.add_argument("input_file")
#parser.add_argument("--id", nargs='?', default="temp")
args = parser.parse_args()

input_file = os.path.abspath(args.input_file)
input_file_name = os.path.basename(input_file)
output_dir = os.path.abspath(os.path.join("..","Results",os.path.splitext(input_file_name)[0],st))
#output_dir = os.path.join(output_dir, args.id)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

current_dir = os.path.dirname(os.path.realpath(__file__))

build_path = os.path.abspath(os.path.join("..", "build-release-script"))

if not os.path.exists(build_path):
    os.makedirs(build_path)

print "Cleaning"
subprocess.call(["make", "clean"], cwd=build_path)
print "Running qmake"
#subprocess.call(["qmake", current_dir, "CONFIG+=nompi"], cwd=build_path)
subprocess.call(["qmake", current_dir], cwd=build_path)
print "Running make"
subprocess.call(["make"], cwd=build_path)

lib_path = os.path.join("..","src")
app_path = os.path.join(build_path,"appconfig")
env = dict(os.environ)
env['LD_LIBRARY_PATH'] = lib_path

#subprocess.call(["./appconfig", input_file, output_dir], cwd=app_path, env=env)
subprocess.call(["mpirun", "-n", "4", "./appconfig", input_file, output_dir], cwd=app_path, env=env)
