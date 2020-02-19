# This code sets up a Cromwell workflow at FredHutch using the SciComp infrastructure


# 1. Set up MariaDB on gizmo: request a MariaDB instance at https://mydb.fredhutch.org/login
# note the settings generated, eg. `mysql --host mydb --port 32217 --user pvan --password`
#
# 2. Clone the Cromwell server skeleton: 
# `git clone https://github.com/FredHutch/diy-cromwell-server cromwell`
#
# 3. edit cromwellParams.sh as appropriate, specifically file paths
# 
# 4. Check the status of your Cromwell job(s), either by `squeue -u mylogin` or
# browsing `http://gizmo[nodenumber]:2020`, where 2020 is the port specified below

library(fh.wdlR)

pworkDir <- "/fh/fast/gottardo_r/pvan_working/"

if (Sys.info()["sysname"] == "Darwin") {
  pworkDir <- "~/rhino"
}

cromwellCreate(FredHutchId = "pvan", 
               port = "2020",
               pathToServerLogs = paste0(pworkDir, "/cromwell/cromwell-serverlogs/%A.txt"),
               pathToServerScript = paste0(pworkDir, "/cromwell/config/cromServer.sh"),
               pathToParams = paste0(pworkDir, "/cromwell/config/cromwellParams.sh")
               )