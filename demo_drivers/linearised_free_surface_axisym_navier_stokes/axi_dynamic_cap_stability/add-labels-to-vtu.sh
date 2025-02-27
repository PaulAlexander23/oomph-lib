#!/bin/zsh

# Ensure correct usage
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 replacements.csv 'file_pattern'"
    echo "Example: $0 replacements.csv 'folder/*.txt'"
    exit 1
fi

CSV_FILE="$1"
FILE_PATTERN="$2"

# Ensure the CSV file exists
if [[ ! -f "$CSV_FILE" ]]; then
    echo "Error: CSV file '$CSV_FILE' not found!"
    exit 1
fi

# Expand the file pattern and check if any files match
FILES=(${~FILE_PATTERN})
if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "Error: No files matched the pattern '$FILE_PATTERN'."
    exit 1
fi

# Read CSV and apply replacements to each file
while IFS=, read -r old new; do
    # Escape special characters in 'old' term
    old_escaped=$(echo "$old" | sed 's/[&/]/\\&/g')

    # Loop through matched files and apply replacements
    for FILE in "${FILES[@]}"; do
        TMP_FILE=$(mktemp)
        sed "s/$old_escaped/$new/g" "$FILE" > "$TMP_FILE" && mv "$TMP_FILE" "$FILE"
        echo "Updated '$FILE': Replaced '$old' with '$new'"
    done
done < "$CSV_FILE"

echo "Replacements completed."
