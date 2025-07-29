import re

# Sanitize column names to remove or replace slashes and parentheses
def sanitize(col):
    col = col.replace("/", "_")        # Replace slash with underscore
    col = col.replace("(", "").replace(")", "")  # Remove parentheses
    return re.sub(r"\s+", "_", col)    # Replace spaces with underscores
