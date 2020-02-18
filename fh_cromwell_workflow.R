# 1. Set up MariaDB on gizmo: request a MariaDB instance at https://mydb.fredhutch.org/login
# note the settings generated, eg. `mysql --host mydb --port 32217 --user pvan --password`
#
# 2. Clone the Cromwell server skeleton: 
# `git clone https://github.com/FredHutch/diy-cromwell-server cromwell`
# edit cromwellParams.sh as appropriate

library(fh.wdlR)
cromwellCreate(FredHutchId = "pvan", port = "2020",
               pathToServerLogs = "/fh/fast/gottardo_r/pvan_working/cromwell/cromwell-serverlogs/%A.txt",
               pathToScript = "/fh/fast/gottardo_r/pvan_working/cromwell/config/cromServer.sh",
               pathToParams = "/fh/fast/gottardo_r/pvan_working/cromwell/config/cromwellParams.sh")