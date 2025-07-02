import argparse
import os
import json
import jsonschema

def validate_json_files(directory, schema_file):
    """
    Validates all JSON files in the specified directory against the given schema.

    Args:
        directory: The directory containing the JSON files.
        schema_file: The path to the JSON schema file.
    """

    with open(schema_file) as f:
        try:
            schema = json.load(f)
        except json.JSONDecodeError as e:
            print(f"schema file {schema_file} did not load")
            print( e )
            quit()


    import os
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".json"):
            filepath = os.path.join(directory, filename)
            #print( "Opening: ", filepath, flush = True )
            with open(filepath) as f:
                try:
                    data = json.load(f)
                except:
                    print(f"{filepath} did not load")
                    quit()

            try:
                jsonschema.validate(instance=data, schema=schema)
                print(f"{filepath} is valid.")
            except jsonschema.exceptions.ValidationError as e:
                print(f"{filepath} is invalid: {e}")

def main():
    parser = argparse.ArgumentParser(description="Process some files.")
    parser.add_argument( "-d", "--directory",
        type=str,
        default="TestRdesJson",  # Set the default directory here
        help="Path to directory containing files (default: TestRdesJson)",
    )
    args = parser.parse_args()

    # Get the directory path from the parsed arguments
    directory = args.directory

    # Ensure the directory exists
    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        exit(1)
    schema_file = "rdesigneurSchema.json"
    validate_json_files(directory, schema_file)

if __name__ == "__main__":
    main()
