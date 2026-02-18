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
            
        # Conflict tracking
        seen_ids = {}
        seen_names = {}
        seen_shas = {}
        errors = 0

        for index, plugin in enumerate(data):
            plugin_id = plugin.get("id")
            plugin_name = plugin.get("name")
            plugin_sha = plugin.get("sha256")

            # Basic field presence check
            if not plugin_id:
                print(f"Error at index {index}: Missing 'id'")
                errors += 1
            if not plugin_name:
                print(f"Error at index {index} ({plugin_id}): Missing 'name'")
                errors += 1
            if not plugin_sha:
                print(f"Error at index {index} ({plugin_id}): Missing 'sha256'")
                errors += 1

            # Conflict checks
            if plugin_id:
                if plugin_id in seen_ids:
                    print(f"Error: Duplicate ID found: '{plugin_id}' (indices {seen_ids[plugin_id]} and {index})")
                    errors += 1
                seen_ids[plugin_id] = index

            if plugin_name:
                if plugin_name in seen_names:
                    print(f"Error: Duplicate Name found: '{plugin_name}' (indices {seen_names[plugin_name]} and {index})")
                    errors += 1
                seen_names[plugin_name] = index

            if plugin_sha:
                if plugin_sha in seen_shas:
                    print(f"Error: Duplicate SHA256 found: '{plugin_sha}' (indices {seen_shas[plugin_sha]} and {index})")
                    print(f"  ID 1: {data[seen_shas[plugin_sha]].get('id')}")
                    print(f"  ID 2: {plugin_id}")
                    errors += 1
                seen_shas[plugin_sha] = index

        if errors > 0:
            print(f"Validation failed with {errors} error(s).")
            return False

        print(f"Success: {file_path} is valid JSON and passed all conflict checks.")
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
    import sys
    
    if len(sys.argv) > 1:
        target_file = sys.argv[1]
    else:
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
