# Setting up Databases

## Establishing MySQL Databases

After installing, start a MySQL session by running ```"sudo mysql"```.

Next, create a database called "patterndetection" by running the following command:

```CREATE DATABASE IF NOT EXISTS patterndetection;```

Create a new user named "patterndetection" by running the following command:

```CREATE USER "zigzag"@"localhost" identified by "MyStrongP@ss123";```

Now grant this user full privileges on the new patterndetection database by running the following command:

```GRANT ALL PRIVILEGES ON patterndetection.* TO 'zigzag'@'localhost';```

Create a second database called "merops" by running the following command:

```CREATE DATABASE IF NOT EXISTS merops;```

Then give privileges to the zigzag user on the merops database by running the following command:

```GRANT ALL PRIVILEGES ON merops.* TO 'zigzag'@'localhost';```

Finally, make sure that the new user privileges take effect by running the following command:

```FLUSH PRIVILEGES;```

## Populating the MySQL Databases

To populate both databases, run the ```populate_dbs.sh``` script located in the ```/data``` directory. This script will use the .sql files to create the necessary tables and populate them with data. Need to ensure the MySQL server is running and the credentials in the script are correct.

To run the script, navigate to the ```/data``` directory and run the following command:

```./populate_dbs.sh```

### Warnings

- [ ] Ensure that the MySQL server is running before running the script.
- [ ] Ensure that the credentials in the script are correct.
- [ ] Ensure that the script has the necessary permissions to run.
- [ ] Ensure the socket path is correct in the script.
- [ ] Ensure the script is executable.
- [ ] Ensure the script is in the correct directory.
