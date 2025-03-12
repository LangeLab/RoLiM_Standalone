#!/bin/bash

# MySQL credentials
MYSQL_USER="zigzag" # Change this to the username used for db
MYSQL_PASSWORD="MyStrongP@ss123" # Change this to your MySQL password
MAIN_FOLDER="tables" # Path to your main folder.
MYSQL_PATH="/run/mysqld/mysqld.sock" # Path to your MySQL socket

# Loop through the subfolders
for database_folder in "$MAIN_FOLDER"/merops "$MAIN_FOLDER"/patterndetection; do
  database_name=$(basename "$database_folder")

  # Loop through SQL files in the subfolder
  for dump_file in "$database_folder"/*.sql; do
    # Check if the dump file exists
    if [ -f "$dump_file" ]; then
      # Populate the database
      mysql -u "$MYSQL_USER" -p"$MYSQL_PASSWORD" --socket="$MYSQL_PATH" "$database_name" < "$dump_file"

      # Check the exit status of the mysql command
      if [ $? -eq 0 ]; then
        echo "Successfully populated database: $database_name from file: $(basename "$dump_file")"
      else
        echo "Error populating database: $database_name from file: $(basename "$dump_file")"
      fi
    else
      echo "Dump file not found: $dump_file"
    fi
  done
done

echo "Database population completed."