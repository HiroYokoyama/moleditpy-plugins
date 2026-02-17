import json
import sys
import os

def validate_json(file_path):
    print(f"Validating {file_path}...")
    
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return False
        
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            
        if not isinstance(data, list):
            print("Error: Root element must be a list (Array).")
            return False
            
        print(f"Success: {file_path} is valid JSON and has a valid structure.")
        return True
        
    except json.JSONDecodeError as e:
        print(f"Error: JSON syntax error in {file_path}:")
        print(f"  Line {e.lineno}, Column {e.colno}: {e.msg}")
        return False
    except Exception as e:
        print(f"Error: An unexpected error occurred while validating {file_path}:")
        print(f"  {str(e)}")
        return False

if __name__ == "__main__":
    # Path relative to the script location or current working directory
    # The script is expected to be in 'utils/' and 'plugins.json' in 'explorer/'
    target_file = os.path.join("explorer", "plugins.json")
    
    # If not found (e.g. running from utils/), try parent directory's explorer/
    if not os.path.exists(target_file):
        target_file = os.path.join("..", "explorer", "plugins.json")

    success = validate_json(target_file)
    if not success:
        sys.exit(1)
    sys.exit(0)
