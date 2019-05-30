# University of Delaware VIP HPC Research

## Ranking Supercomputers using SPEC Benchmarks

Gather real-world application benchmark suite results from SPEC HPG 2020 using a web-scraper. Store these results into our database on university server (yoda) so website can extract this information and have a sorting/filter system for users to view metrics.

### Web-Scraper

The process by which data will be gathered(currently):

1. Run ListLinks program to create notepad file with download links to all .csv for a given test suite.
2. Run wget -i on the file with all of the links in the directory where we want every raw .csv to be place
3. Run FoldAllCreator program in order to create a .sh file that is just “java -jar csvscraper1.4.jar” followed by all of the file names separated by spaces.
4. Run the foldAll.sh file that we created, which runs the "FileFolder" program on all raw input.
5. Now every .csv file is scraped, with a new output, and is ready to be fed to the database.

NOTE: In the future, these programs will exist in a cron job, and will have directory organization built in, so that they are able to constantly pull information and feed it into the database. Currently all the scripts are manually executed.

### Database

MongoDB database installed on yoda. Run python scripts to move scraped .csv files into the database.

1. Pull the vip github since the webscraped data is stored on github.
2. get the path to the new ouput csv files, change the path in each script based on the benchmarks, located in the importDB folder.
3. run each script using "python <ScriptName.py>" once
4. Now everthing is in the database.

NOTE: In the future, these programs will exist in a cron job, and will have directory organization built in, so that they are able to constantly updating database. Currently all the scripts are manually executed.
### User-Interface

*subject to change* 

Currently Yoda is housing the Frontend and the Backend. The Frontend is roughly defined as the React Website and the Backend as the Nodemon (looking to change this) website server on port 3000.

### Contributors

**Advisors:** Sunita Chandrasekaran, Rudolf Eigenmann, Mayara Gimenes,

**Web-Scraping Team:** Derek Baum, Ryan Emenheiser, Matt Benvenuto

**Database Team:** Max Luu, Jake Wise

**User-Interface Team:** Matthew Stack, Chris Munley
