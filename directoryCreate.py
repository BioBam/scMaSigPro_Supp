import os
import json
import sys

# Check if the input file name has been provided
if len(sys.argv) != 2:
    print("Please provide the JSON file as an argument.")
    sys.exit()

# Load the JSON file
with open(sys.argv[1], 'r') as f:
    directories = json.load(f)

def create_directories(base_path, directories):
    for key, value in directories.items():
        new_path = os.path.join(base_path, key)
        os.makedirs(new_path, exist_ok=True)
        print(f"Created directory: {new_path}")
        if isinstance(value, dict):
            create_directories(new_path, value)

# Start creating directories
create_directories(".", directories)
