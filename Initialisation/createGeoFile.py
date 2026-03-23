import math

def read_parameters_from_lua(lua_file, keys):
    """
    Reads specified keys (e.g., GL, L) from a Lua file with equations.
    Returns a dictionary {key: value}.
    """
    variables = {}
    results = {}

    with open(lua_file, "r") as f:
        for line in f:
            # Remove comments
            line = line.split("--")[0].strip()
            if not line or "=" not in line:
                continue

            key, expr = line.split("=", 1)
            key = key.strip()
            expr = expr.strip().rstrip(";")  # remove Lua semicolon
            expr = expr.replace("^", "**")  # Lua power -> Python

            # Evaluate
            value = eval(expr, {"__builtins__": None, "math": math}, variables)
            variables[key] = value

            # Capture only requested keys
            if key in keys:
                results[key] = value

    # Ensure all keys were found
    for k in keys:
        if k not in results:
            raise ValueError(f"Key '{k}' not found in Lua file")

    return results

def populate_geo_file(template_file, output_file, params):
    """
    Replaces placeholders in a .geo file with values from params dict.
    Placeholders are like {GL}, {L}, etc.
    """
    with open(template_file, "r") as f:
        content = f.read()

    for key, value in params.items():
        content = content.replace(f"{{{key}}}", str(value))

    with open(output_file, "w") as f:
        f.write(content)


lua_file = "parameters.lua"
geo_template = "Initialisation/template.geo"
#geo2_template = "Initialisation/template2D.geo"

geo_output = "footprint.geo"
#geo2_output = "footprint2D.geo"

# Parameters you want to extract from Lua
keys = ["GL", "L", "GLRes","BoundaryRes"]




params = read_parameters_from_lua(lua_file, keys)
populate_geo_file(geo_template, geo_output, params)
#populate_geo_file(geo2_template, geo2_output, params)
